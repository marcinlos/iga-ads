#ifndef PROBLEMS_STOKES_STOKES_DG_HPP_
#define PROBLEMS_STOKES_STOKES_DG_HPP_

#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "problems/stokes/space_set.hpp"
#include "mumps.hpp"



namespace ads {

class stokes_dg : public simulation_2d {
private:
    using Base = simulation_2d;

    struct value_pair {
        value_type v1;
        value_type v2;
        bool boundary;
    };

    galois_executor executor{8};

    space_set trial, test;

    double h;

    double Cpen;// = 5 * (3 + 1);
    double hF;// = 1 / 20.;

    double eta;

    mumps::solver solver;

    Galois::StatTimer solver_timer{"solver"};
    Galois::StatTimer assembly_timer{"assembly"};
    Galois::StatTimer rhs_timer{"rhs"};

    output_manager<2> outputU1, outputU2, outputP;

public:
    stokes_dg(space_set trial_, space_set test_, const timesteps_config& steps)
    : Base{ test_.Px, test_.Py, steps }
    , trial{ std::move(trial_) }
    , test{ std::move(test_) }
    , h{ element_diam(trial.Px, trial.Py) }
    , outputU1{ trial.U1x.B, trial.U1y.B, 500 }
    , outputU2{ trial.U2x.B, trial.U2y.B, 500 }
    , outputP{ trial.Px.B, trial.Py.B, 500 }
    {
        // 5(p + 1)
        Cpen = 5 * trial.U1x.B.degree;
        hF = h;//1. / trial.Px.B.elements();
        auto p = trial.Px.B.degree;
        // eta = 3 * (p + 1) * (p + 2);
        eta = 10 * (p + 1) * (p + 2);
    }

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    void before() override {
        trial.U1x.factorize_matrix();
        trial.U1y.factorize_matrix();
        trial.U2x.factorize_matrix();
        trial.U2y.factorize_matrix();
        trial.Px.factorize_matrix();
        trial.Py.factorize_matrix();

        test.U1x.factorize_matrix();
        test.U1y.factorize_matrix();
        test.U2x.factorize_matrix();
        test.U2y.factorize_matrix();
        test.Px.factorize_matrix();
        test.Py.factorize_matrix();

        // output_exact();
    }

    value_type exact_p(point_type p) const {
        // auto x = p[0];
        // return {x * (1 - x) - 1./6, 1 - 2 * x, 0.0};

        // non-polynomial
        using std::exp;
        using std::pow;
        auto x = p[0], y = p[1];
        auto xx = x * x, yy = y * y;
        auto ex = exp(x);
        auto e = exp(1);

        return {
            -424 + 156*e + (yy - y) * (-456 + ex * (456 + xx * (228 - 5 *(yy - y)) + 2*x*(-228 + (yy - y)) + 2*pow(x,3) * (-36 + (yy - y)) + pow(x,4) * (12 + (yy - y)))),
            ex * (y - 1) * y * (pow(x,4) * (yy - y + 12) + 6 * pow(x,3) * (yy - y - 4) + xx * (yy - y + 12) - 8 * x * (y - 1) * y + 2 * (y - 1) * y),
            2 * (2 * y - 1) * (ex * (pow(x,4) * (yy - y + 6) + 2 * pow(x,3) * (yy - y - 18) + xx * (-5 * yy + 5 * y + 114) + 2 * x * (yy - y - 114) + 228) - 228)
        };
    }

    std::array<value_type, 2> exact_v(point_type p) const {
        // auto f = [](double x, double y) {
        //     return x*x * (1 - x) * (1 - x) * (2 * y - 6 * y*y + 4 * y*y*y);
        // };

        // auto dfx = [](double x, double y) {
        //     return (4 * x*x*x - 6 * x*x + 2 * x) * (2 * y - 6 * y*y + 4 * y*y*y);
        // };

        // auto dfy = [](double x, double y) {
        //     return x*x * (1 - x) * (1 - x) * (2 - 12 * y + 12 * y*y);
        // };

        // double x = p[0], y = p[1];
        // value_type vx = {f(x, y), dfx(x, y), dfy(x, y)};
        // value_type vy = {-f(y, x), -dfy(y, x), -dfx(y, x)};

        // return { vx ,vy };

        // non-polynomial
        using std::exp;
        using std::pow;
        auto x = p[0], y = p[1];
        auto ex = exp(x);

        auto ux = value_type{
            2*ex * pow(-1 + x, 2) * x * x * (y * y - y) * (-1 + 2*y),
            2*ex * x * (pow(x,3) + 2 * x * x - 5 * x + 2) * y*  (2 * y * y - 3 * y + 1),
            2*ex * pow(x - 1, 2) * x * x * (6 *  y * y - 6 * y + 1)
        };

        auto uy = value_type{
            -ex * (-1 + x) * x* (-2 + x * (3 + x)) * pow(-1 + y, 2) * y * y,
            -ex * (pow(x,4) + 6 * pow(x,3) + x * x - 8 * x + 2) * pow(y - 1, 2) * y * y,
            -2 * ex * x * (pow(x,3) + 2 * x * x - 5 *  x + 2) * y * (2 * y * y - 3 * y + 1)
        };

        return {ux, uy};
    }

    value_type exact_div(point_type p) const {
        return {0, 0, 0};
    }

    void output_exact() {
        auto p = [this](point_type x) { return exact_p(x).val; };
        auto vx = [this](point_type x) { return exact_v(x)[0].val; };
        auto vy = [this](point_type x) { return exact_v(x)[1].val; };

        auto project = [&](auto& x, auto& y, auto fun) {
            vector_type rhs{{ x.dofs(), y.dofs() }};
            vector_type buffer{{ x.dofs(), y.dofs() }};
            compute_projection(rhs, x.basis, y.basis, [&](double x, double y) { return fun({x, y}); });
            ads_solve(rhs, buffer, x.data(), y.data());
            return rhs;
        };

        outputP.to_file(project(trial.Px, trial.Py, p), "pressure_ref.data");
        outputU1.to_file(project(trial.U1x, trial.U1y, vx), "vx_ref.data");
        outputU2.to_file(project(trial.U2x, trial.U2y, vy), "vy_ref.data");
    }

    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p,
                     const vector_view& Rvx, const vector_view& Rvy, const vector_view& Rp) const {
        auto e_vx = [this](point_type x) { return exact_v(x)[0]; };
        auto e_vy = [this](point_type x) { return exact_v(x)[1]; };
        auto e_p = [this](point_type x) { return exact_p(x); };
        auto div = [this](point_type x) { return exact_div(x); };

        double vxL2 = errorL2(vx, trial.U1x, trial.U1y, e_vx) / normL2(trial.U1x, trial.U1y, e_vx) * 100;
        double vxH1 = errorH1(vx, trial.U1x, trial.U1y, e_vx) / normH1(trial.U1x, trial.U1y, e_vx) * 100;

        double vyL2 = errorL2(vy, trial.U2x, trial.U2y, e_vy) / normL2(trial.U2x, trial.U2y, e_vy) * 100;
        double vyH1 = errorH1(vy, trial.U2x, trial.U2y, e_vy) / normH1(trial.U2x, trial.U2y, e_vy) * 100;

        double pL2 = errorL2(p, trial.Px, trial.Py, e_p) / normL2(trial.Px, trial.Py, e_p) * 100;
        double pH1 = errorH1(p, trial.Px, trial.Py, e_p) / normH1(trial.Px, trial.Py, e_p) * 100;

        double divL2 = div_errorL2(vx, vy, trial, div) * 100;
        double divH1 = div_errorH1(vx, vy, trial, div) * 100;

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
        std::cout << "p   : L2 = " << pL2    << "%, H1 = " << pH1    << "%" << std::endl;
        std::cout << "div : L2 = " << divL2  << ", H1 = " << divH1  << std::endl;

        auto H10 = [](value_type a) { return a.dx * a.dx + a.dy * a.dy; };

        double err_vx = error(vx, trial.U1x, trial.U1y, H10, e_vx);
        double err_vy = error(vy, trial.U2x, trial.U2y, H10, e_vy);
        double err_p  = errorL2(p, trial.Px, trial.Py, e_p);
        double error_Vh = std::sqrt(err_vx * err_vx + err_vy * err_vy + err_p * err_p);

        std::cout << "Error in Vh:   " << error_Vh << std::endl;

        double r_norm = normV(Rvx, Rvy, Rp);
        std::cout << "Residual norm: " << r_norm << std::endl;

    }

    point_type forcing(point_type p) const {
        // double x = p[0], y = p[1];

        // auto fx =
        //     (12 - 24 * y) * x*x*x*x +
        //     (-24 + 48 * y) * x*x*x +
        //     (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
        //     (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
        //     1 - 4 * y + 12 * y*y - 8 * y*y*y;

        // auto fy =
        //     (8 - 48 * y + 48 * y*y) * x*x*x +
        //     (-12 + 72 * y - 72 * y*y) * x*x +
        //     (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x -
        //     12 * y*y + 24 * y*y*y - 12 * y*y*y*y;

        // return { fx, fy };

        // cavity flow
        // return { 0, 0 };

        // non-polynomial
        using std::exp;
        using std::pow;
        auto x = p[0], y = p[1];
        auto xx = x * x, yy = y * y;
        auto ex = exp(x);

        auto px = ex * (y - 1) * y * (pow(x,4) * (yy - y + 12) + 6 * pow(x,3) * (yy - y - 4) + xx * (yy - y + 12) - 8 * x * (y - 1) * y + 2 * (y - 1) * y);
        auto py = 2 * (2 * y - 1) * (ex * (pow(x,4) * (yy - y + 6) + 2 * pow(x,3) * (yy - y - 18) + xx * (-5 * yy + 5 * y + 114) + 2 * x * (yy - y - 114) + 228) - 228);

        auto Lux = 2 * ex * (pow(x,4) * (2 * pow(y,3) - 3 * yy + 13 * y - 6) + 6 * pow(x,3) * (2 * pow(y,3) - 3 * yy - 3 * y + 2) +
                             xx * (2 * pow(y,3) - 3 * yy + 13 * y - 6) - 8 * x * y * (2 * yy - 3 * y + 1) + 2 * y * (2 * yy - 3 * y + 1));
        auto Luy = -ex * (pow(x,4) * (pow(y,4) - 2 * pow(y,3) + 13 * yy - 12 * y + 2) + 2 * pow(x,3) * (5 * pow(y,4) - 10 * pow(y,3) + 17 * yy - 12 * y + 2) +
                          xx * (19 * pow(y,4) - 38 * pow(y,3) - 41 * yy + 60 * y - 10) + x * (-6 * pow(y,4) + 12 * pow(y,3) + 18 * yy - 24 * y + 4) - 6 * pow(y - 1, 2) * yy);

        return {-Lux + px, -Luy + py};
    }

    auto shifted(int n, int k, mumps::problem& problem) const {
        return [&problem,n,k](int i, int j, double val) {
            problem.add(n + i, k + j, val);
        };
    }

    bool overlap(int a, const dimension& U, int b, const dimension& V) const {
        auto ar = U.basis.element_ranges[a];
        auto br = V.basis.element_ranges[b];
        return (ar.first >= br.first && ar.first <= br.second) || (br.first >= ar.first && br.first <= ar.second);
    }

    bool overlap(index_type a, const dimension& Ux, const dimension& Uy,
                 index_type b, const dimension& Vx, const dimension& Vy) const {
        return overlap(a[0], Ux, b[0], Vx) && overlap(a[1], Uy, b[1], Vy);
    }

    bool dofs_touch(int a, const dimension& U, int b, const dimension& V) const {
        auto ar = U.basis.element_ranges[a];
        auto br = V.basis.element_ranges[b];
        return (ar.first >= br.first && ar.first <= br.second + 1) || (br.first >= ar.first && br.first <= ar.second + 1);
    }

    bool dofs_touch(index_type a, const dimension& Ux, const dimension& Uy,
                    index_type b, const dimension& Vx, const dimension& Vy) const {
        return dofs_touch(a[0], Ux, b[0], Vx) && dofs_touch(a[1], Uy, b[1], Vy);
    }

    bool is_pressure_fixed(index_type dof) const {
        return dof[0] == 0 && dof[1] == 0;
        // return false;
    }

    template <typename Form>
    double integrate(index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                     const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;

        for (auto e : elements_supporting_dof(i, Ux, Uy)) {
            if (! supported_in(j, e, Vx, Vy)) continue;

            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q);
                value_type ww = eval_basis(e, q, i, Ux, Uy);
                value_type uu = eval_basis(e, q, j, Vx, Vy);

                double fuw = form(ww, uu);
                val += fuw * w * J;
            }
        }
        return val;
    }

    auto touching_dofs(int dof, int begin, int end, const dimension& x) const {
        using std::min;
        using std::max;

        auto elems = x.basis.element_ranges[dof];
        auto first = x.basis.first_dof(max(elems.first - 1, 0));
        auto last = x.basis.last_dof(min(elems.second + 1, x.elements - 1));

        auto minx = max(begin, first);
        auto maxx = min(end, last);

        return boost::counting_range(minx, maxx + 1);
    }

    index_range touching_dofs(index_type dof, const dimension& x, const dimension& y) const {
        auto rx = touching_dofs(dof[0], 0, x.dofs(), x);
        auto ry = touching_dofs(dof[1], 0, y.dofs(), y);
        return util::product_range<index_type>(rx, ry);
    }

    void assemble_matrix(mumps::problem& problem) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dP = trial.Px.dofs() * trial.Py.dofs();

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();

        auto D = DU1 + DU2 + DP;

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DU1, DU1, problem);
        auto test_p = shifted(DU1 + DU2, DU1 + DU2, problem);

        // auto test_vxy = shifted(0, DU1, problem);
        // auto test_vyx = shifted(DU1, 0, problem);

        auto trial_vx = shifted(D, D, problem);
        auto trial_vy = shifted(D + dU1, D + dU1, problem);
        auto trial_p = shifted(D + dU1 + dU2, D + dU1 + dU2, problem);

        auto N = D + dU1 + dU2 + dP;
        bool bc = false;
        bool fix_p = false;

        // auto hh = h * h;

        // Gram matrix
        // G(w, u)
        // w = (tx, ty, w)
        // u = (vx, vy, p)

        // w, p -> (w, p) +  hh ([w], [p])
        for (auto i : dofs(test.Px, test.Py)) {
            // for (auto j : overlapping_dofs(i, test.Px, test.Py)) {
            for (auto j : touching_dofs(i, test.Px, test.Py)) {
                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, test.Px, test.Py) + 1;

                if (!fix_p || (! is_pressure_fixed(i) && ! is_pressure_fixed(j))) {
                    auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };
                    double product = eval([](auto w, auto p) { return w.val * p.val; });

                    auto form = [this](auto w, auto p, auto, auto) {
                        return h * jump(w).val * jump(p).val;
                    };
                    double boundary = integrate_over_internal_skeleton(i, j, test.Px, test.Py, test.Px, test.Py, form);
                    test_p(ii, jj, product + boundary);
                }
            }
        }

        // tx, vx -> (\/tx, \/vx) + 1/h ([tx], [vx])
        for (auto i : dofs(test.U1x, test.U1y)) {
            // for (auto j : overlapping_dofs(i, test.U1x, test.U1y)) {
            for (auto j : touching_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                // Weak BC
                // if (! is_boundary(i[0], test.U1x) && ! is_boundary(j[0], test.U1x)) {
                // Strong BC
                if (!bc || (! is_boundary(i, test.U1x, test.U1y) && ! is_boundary(j, test.U1x, test.U1y))) {
                    double val = eval([this](auto tx, auto vx) { return grad_dot(tx, vx); });
                    // skeleton
                    auto form = [this](auto tx, auto vx, auto, auto) {
                        return eta/h * jump(tx).val * jump(vx).val;
                    };
                    val += integrate_over_skeleton(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form);

                    test_vx(ii, jj, val);
                }
            }
        }

        // ty, vy -> (\/ty, \/vy) + 1/h ([ty], [vy])
        for (auto i : dofs(test.U2x, test.U2y)) {
            // for (auto j : overlapping_dofs(i, test.U2x, test.U2y)) {
            for (auto j : touching_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                // Weak BC
                // if (! is_boundary(i[1], test.U2y) && ! is_boundary(j[1], test.U2y)) {
                // Strong BC
                if (!bc || (! is_boundary(i, test.U2x, test.U2y) && ! is_boundary(j, test.U2x, test.U2y))) {
                    double val = eval([this](auto ty, auto vy) { return grad_dot(ty, vy); });
                    // skeleton
                    auto form = [this](auto ty, auto vy, auto, auto) {
                        return eta/h * jump(ty).val * jump(vy).val;
                    };
                    val += integrate_over_skeleton(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form);

                    test_vy(ii, jj, val);
                }
            }
        }

        // Dirichlet BC - test space
        // Weak BC
        // for (auto iy = 0; iy < test.U1y.dofs(); ++ iy) {
        //     int i0 = linear_index({0, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i0, i0, 1);
        //     int i1 = linear_index({test.U1x.dofs() - 1, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i1, i1, 1);
        // }
        // for (auto ix = 0; ix < test.U2x.dofs(); ++ ix) {
        //     int i0 = linear_index({ix, 0}, test.U2x, test.U2y) + 1;
        //     test_vy(i0, i0, 1);
        //     int i1 = linear_index({ix, test.U2y.dofs() - 1}, test.U2x, test.U2y) + 1;
        //     test_vy(i1, i1, 1);
        // }
        // Strong BC
        if (bc) {
            for_boundary_dofs(test.U1x, test.U1y, [&](index_type dof) {
                int i = linear_index(dof, test.U1x, test.U1y) + 1;
                test_vx(i, i, 1);
            });
            for_boundary_dofs(test.U2x, test.U2y, [&](index_type dof) {
                int i = linear_index(dof, test.U2x, test.U2y) + 1;
                test_vy(i, i, 1);
            });
        }
        if (fix_p) {
            int i = linear_index({0, 0}, test.Px, test.Py) + 1;
            test_p(i, i, 1.0);
        }


        // B, B^t
        auto put = [&](int i, int j, int si, int sj, double val, bool fixed_i, bool fixed_j) {
            int ii = i + si;
            int jj = j + sj;

            if (!fixed_i) {
                problem.add(ii, D + jj, -val);
            }
            if (!fixed_i && !fixed_j) {
                problem.add(D + jj, ii, val);
            }
        };

        // B(u, v)
        // u = (ux, uy, p)  trial
        // v = (vx, vy, q)  test

        // auto pen_term = [this](auto w, auto u) { return - 2 * Cpen / hF * w.val * u.val; };

        // ux, vx -> (\/ux, \/vx) - ({\/ux}, [vx])
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! dofs_touch(i, test.U1x, test.U1y, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[0], test.U1x);
                // bool bd_j = is_boundary(j[0], trial.U1x);
                // Strong BC
                bool bd_i = is_boundary(i, test.U1x, test.U1y) && bc;
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y) && bc;

                double value = eval([this](auto vx, auto ux) { return grad_dot(vx, ux); });

                // skeleton
                auto form = [this](auto vx, auto ux, auto, auto n) {
                    return eta / hF * jump(ux).val * jump(vx).val
                         - dot(average(ux), n) * jump(vx).val
                         - dot(average(vx), n) * jump(ux).val;
                };
                value += integrate_over_skeleton(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form);

                // boundary terms
                // top/bottom
                // if (is_boundary(i[1], test.U1y) || is_boundary(j[1], trial.U1y)) {
                //     auto side = i[1] == 0 ? boundary::bottom : boundary::top;
                //     // int n = side == boundary::top ? 1 : -1;
                //     value -= integrate_boundary(side, i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, [&](auto w, auto u, auto) {
                //         return pen_term(w, u);
                //     });
                // }
                put(ii, jj, 0, 0, value, bd_i, bd_j);
            }
        }

        // uy, vy -> (\/uy, \/vy) - ({\/uy}, [vy])
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! dofs_touch(i, test.U2x, test.U2y, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[1], test.U2y);
                // bool bd_j = is_boundary(j[1], trial.U2y);
                // Strong BC
                bool bd_i = is_boundary(i, test.U2x, test.U2y) && bc;
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y) && bc;

                double value = eval([this](auto vy, auto uy) { return grad_dot(vy, uy); });

                // skeleton
                auto form = [this](auto vy, auto uy, auto, auto n) {
                    return eta / hF * jump(uy).val * jump(vy).val
                         - dot(average(uy), n) * jump(vy).val
                         - dot(average(vy), n) * jump(uy).val;
                };
                value += integrate_over_skeleton(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form);

                // boundary terms
                // left/right
                // if (is_boundary(i[0], test.U2x) || is_boundary(j[0], trial.U2x)) {
                //     auto side = i[0] == 0 ? boundary::left : boundary::right;
                //     // int n = side == boundary::right ? 1 : -1;
                //     value -= integrate_boundary(side, i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, [&](auto w, auto u, auto) {
                //         return pen_term(w, u);
                //     });
                // }
                put(ii, jj, DU1, dU1, value, bd_i, bd_j);
            }
        }

        // ux, q -> (ux,x, q)- [ux] n1 {q}
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U1x, trial.U1y, form); };

                bool fixed_i = is_pressure_fixed(i) && fix_p;
                // Weak BC
                // bool bd_j = is_boundary(j[0], trial.U1x);
                // Strong BC
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y) && bc;

                double val = eval([](auto q, auto ux) { return q.val * ux.dx; });

                // skeleton
                auto form = [this](auto q, auto ux, auto, auto n) {
                    return - average(q).val * n[0] * jump(ux).val;
                };
                val += integrate_over_skeleton(i, j, test.Px, test.Py, trial.U1x, trial.U1y, form);

                put(ii, jj, DU1 + DU2, 0, val, fixed_i, bd_j);
            }
        }

        // uy, q -> (uy,y, q) - [uy] n2 {q}
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U2x, trial.U2y, form); };

                bool fixed_i = is_pressure_fixed(i) && fix_p;
                // Weak BC
                // bool bd_j = is_boundary(j[1], trial.U2y);
                // Strong BC
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y) && bc;

                double val = eval([](auto q, auto uy) { return q.val * uy.dy; });

                // skeleton
                auto form = [this](auto q, auto uy, auto, auto n) {
                    return - average(q).val * n[1] * jump(uy).val;
                };
                val += integrate_over_skeleton(i, j, test.Px, test.Py, trial.U2x, trial.U2y, form);

                put(ii, jj, DU1 + DU2, dU1, val, fixed_i, bd_j);
            }
        }

        // p, vx ->  - (p, vx,x) + hF [vx] n1 {p}
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.U1x, test.U1y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.Px, trial.Py, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[0], test.U1x);
                // Strong BC
                bool bd_i = is_boundary(i, test.U1x, test.U1y) && bc;

                bool fixed_j = is_pressure_fixed(j) && fix_p;

                double val = eval([](auto vx, auto p) { return -vx.dx * p.val; });
                // skeleton
                auto form = [this](auto vx, auto p, auto, auto n) {
                    return average(p).val * n[0] * jump(vx).val;
                };
                val += integrate_over_skeleton(i, j, test.U1x, test.U1y, trial.Px, trial.Py, form);

                put(ii, jj, 0, dU1 + dU2, val, bd_i, fixed_j);
            }
        }

        // p, vy ->  - (p, vy,y) + hF [vy] n2 {p}
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.U2x, test.U2y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.Px, trial.Py, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[1], test.U2y);
                // Strong BC
                bool bd_i = is_boundary(i, test.U2x, test.U2y) && bc;

                bool fixed_j = is_pressure_fixed(j) && fix_p;

                double val = eval([](auto vy, auto p) { return -vy.dy * p.val; });
                // skeleton
                auto form = [this](auto vy, auto p, auto, auto n) {
                    return average(p).val * n[1] * jump(vy).val;
                };
                val += integrate_over_skeleton(i, j, test.U2x, test.U2y, trial.Px, trial.Py, form);

                put(ii, jj, DU1, dU1 + dU2, val, bd_i, fixed_j);
            }
        }

        // Stabilization term
        // S(p, q) = hF [p][q]
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                // auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.Px, trial.Py, form); };

                bool fixed_i = is_pressure_fixed(i) && fix_p;
                bool fixed_j = is_pressure_fixed(j) && fix_p;

                // skeleton
                auto form = [this](auto q, auto p, auto, auto) {
                    return h * jump(p).val * jump(q).val;
                };
                double val = integrate_over_internal_skeleton(i, j, test.Px, test.Py, trial.Px, trial.Py, form);
                // val += h*h*eval([this](auto q, auto p) { return grad_dot(p, q); }); // Minev

                put(ii, jj, DU1 + DU2, dU1 + dU2, val, fixed_i, fixed_j);
            }
        }

        // Lagrange multiplier
        // for (auto i : dofs(trial.Px, trial.Py)) {
        //     int ii = linear_index(i, trial.Px, trial.Py) + 1;

        //     auto eval = [&](auto form) { return integrate(i, i, trial.Px, trial.Py, trial.Px, trial.Py, form); };
        //     auto val = eval([this](auto q, auto p) { return p.val; });

        //     problem.add(N + 1, D + dU1 + dU2 + ii, val);
        //     problem.add(D + dU1 + dU2 + ii, N + 1, val);
        // }

        // Dirichlet BC - trial space
        // Weak BC
        // for (auto iy = 0; iy < trial.U1y.dofs(); ++ iy) {
        //     int i0 = linear_index({0, iy}, trial.U1x, trial.U1y) + 1;
        //     trial_vx(i0, i0, 1);
        //     int i1 = linear_index({trial.U1x.dofs() - 1, iy}, trial.U1x, trial.U1y) + 1;
        //     trial_vx(i1, i1, 1);
        // }
        // for (auto ix = 0; ix < trial.U2x.dofs(); ++ ix) {
        //     int i0 = linear_index({ix, 0}, trial.U2x, trial.U2y) + 1;
        //     trial_vy(i0, i0, 1);
        //     int i1 = linear_index({ix, trial.U2y.dofs() - 1}, trial.U2x, trial.U2y) + 1;
        //     trial_vy(i1, i1, 1);
        // }
        // Strong BC
        if (bc) {
            for_boundary_dofs(trial.U1x, trial.U1y, [&](index_type dof) {
                int i = linear_index(dof, trial.U1x, trial.U1y) + 1;
                trial_vx(i, i, 1);
            });
            for_boundary_dofs(trial.U2x, trial.U2y, [&](index_type dof) {
                int i = linear_index(dof, trial.U2x, trial.U2y) + 1;
                trial_vy(i, i, 1);
            });
        }
        if (fix_p) {
            int ii = linear_index({0, 0}, trial.Px, trial.Py) + 1;
            trial_p(ii, ii, 1.0);
        }
    }



    void assemble_gram_matrix_no_jumps(mumps::problem& problem) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dP = trial.Px.dofs() * trial.Py.dofs();

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();

        auto D = DU1 + DU2 + DP;

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DU1, DU1, problem);
        auto test_p = shifted(DU1 + DU2, DU1 + DU2, problem);

        auto trial_vx = shifted(D, D, problem);
        auto trial_vy = shifted(D + dU1, D + dU1, problem);
        auto trial_p = shifted(D + dU1 + dU2, D + dU1 + dU2, problem);

        auto N = D + dU1 + dU2 + dP;
        bool bc = false;
        bool fix_p = false;

        // Gram matrix
        // G(w, u)
        // w = (tx, ty, w)
        // u = (vx, vy, p)

        // w, p -> (w, p)
        for (auto i : dofs(test.Px, test.Py)) {
            // for (auto j : overlapping_dofs(i, test.Px, test.Py)) {
            for (auto j : touching_dofs(i, test.Px, test.Py)) {
                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, test.Px, test.Py) + 1;

                if (!fix_p || (! is_pressure_fixed(i) && ! is_pressure_fixed(j))) {
                    auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };
                    double product = eval([](auto w, auto p) { return w.val * p.val; });
                    test_p(ii, jj, product);
                }
            }
        }

        // tx, vx -> (\/tx, \/vx)
        for (auto i : dofs(test.U1x, test.U1y)) {
            // for (auto j : overlapping_dofs(i, test.U1x, test.U1y)) {
            for (auto j : touching_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                // Weak BC
                // if (! is_boundary(i[0], test.U1x) && ! is_boundary(j[0], test.U1x)) {
                // Strong BC
                if (!bc || (! is_boundary(i, test.U1x, test.U1y) && ! is_boundary(j, test.U1x, test.U1y))) {
                    double val = eval([this](auto tx, auto vx) { return grad_dot(tx, vx); });
                    test_vx(ii, jj, val);
                }
            }
        }

        // ty, vy -> (\/ty, \/vy)
        for (auto i : dofs(test.U2x, test.U2y)) {
            // for (auto j : overlapping_dofs(i, test.U2x, test.U2y)) {
            for (auto j : touching_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                // Weak BC
                // if (! is_boundary(i[1], test.U2y) && ! is_boundary(j[1], test.U2y)) {
                // Strong BC
                if (!bc || (! is_boundary(i, test.U2x, test.U2y) && ! is_boundary(j, test.U2x, test.U2y))) {
                    double val = eval([this](auto ty, auto vy) { return grad_dot(ty, vy); });
                    test_vy(ii, jj, val);
                }
            }
        }

        // Dirichlet BC - test space
        // Weak BC
        // for (auto iy = 0; iy < test.U1y.dofs(); ++ iy) {
        //     int i0 = linear_index({0, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i0, i0, 1);
        //     int i1 = linear_index({test.U1x.dofs() - 1, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i1, i1, 1);
        // }
        // for (auto ix = 0; ix < test.U2x.dofs(); ++ ix) {
        //     int i0 = linear_index({ix, 0}, test.U2x, test.U2y) + 1;
        //     test_vy(i0, i0, 1);
        //     int i1 = linear_index({ix, test.U2y.dofs() - 1}, test.U2x, test.U2y) + 1;
        //     test_vy(i1, i1, 1);
        // }
        // Strong BC
        if (bc) {
            for_boundary_dofs(test.U1x, test.U1y, [&](index_type dof) {
                int i = linear_index(dof, test.U1x, test.U1y) + 1;
                test_vx(i, i, 1);
            });
            for_boundary_dofs(test.U2x, test.U2y, [&](index_type dof) {
                int i = linear_index(dof, test.U2x, test.U2y) + 1;
                test_vy(i, i, 1);
            });
        }
        if (fix_p) {
            int i = linear_index({0, 0}, test.Px, test.Py) + 1;
            test_p(i, i, 1.0);
        }
    }


    void compute_rhs(vector_view& vx, vector_view& vy, vector_view& p) const {
        using shape = std::array<std::size_t, 2>;
        auto u1_shape = shape{ test.U1x.basis.dofs_per_element(), test.U1y.basis.dofs_per_element() };
        auto u2_shape = shape{ test.U2x.basis.dofs_per_element(), test.U2y.basis.dofs_per_element() };


        executor.for_each(elements(test.Px, test.Py), [&](index_type e) {
            auto vx_loc = vector_type{ u1_shape };
            auto vy_loc = vector_type{ u2_shape };

            double J = jacobian(e);
            for (auto q : quad_points(test.Px, test.Py)) {
                double W = weigth(q);
                auto x = point(e, q);
                auto F = forcing(x);
                for (auto a : dofs_on_element(e, test.U1x, test.U1y)) {
                    auto aa = dof_global_to_local(e, a, test.U1x, test.U1y);
                    value_type v = eval_basis(e, q, a, test.U1x, test.U1y);

                    double Lvx = F[0] * v.val;
                    vx_loc(aa[0], aa[1]) -= Lvx * W * J;
                }
                for (auto a : dofs_on_element(e, test.U2x, test.U2y)) {
                    auto aa = dof_global_to_local(e, a, test.U2x, test.U2y);
                    value_type v = eval_basis(e, q, a, test.U2x, test.U2y);

                    double Lvy = F[1] * v.val;
                    vy_loc(aa[0], aa[1]) -= Lvy * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(vx, vx_loc, e, test.U1x, test.U1y);
                update_global_rhs(vy, vy_loc, e, test.U2x, test.U2y);
            });
        });

        // // cavity flow
        // auto side = boundary::top;

        // // vx
        // for (auto i : dofs(test.U1x, test.U1y)) {
        //     double val = 0;

        //     auto form = [&](auto v, auto x, auto n) {
        //         auto g = 1;
        //         return - dot(v, n) * g + eta / hF * v.val * g;
        //     };
        //     if (touches(i, side, test.U1x, test.U1y)) {
        //         val += integrate_boundary(side, i, test.U1x, test.U1y, form);
        //     }
        //     vx(i[0], i[1]) -= val;
        // }
        // // vy - 0

        // // p
        // for (auto i : dofs(test.Px, test.Py)) {
        //     double val = 0;

        //     auto form = [&](auto q, auto x, auto n) {
        //         auto g = 1;
        //         return q.val * n[0];
        //     };
        //     if (touches(i, side, test.Px, test.Py)) {
        //         val += integrate_boundary(side, i, test.Px, test.Py, form);
        //     }
        //     p(i[0], i[1]) -= val;
        // }
    }


    void apply_bc(vector_view& vx, vector_view& vy, vector_view& p) {
        // Strong BC
        // zero_bc(vx, trial.U1x, trial.U1y);
        // zero_bc(vy, trial.U2x, trial.U2y);

        // Cavity flow
        // constexpr double eps = 1e-4;
        // auto drop = [](double t) { return t < 1 - eps ? 0 : 1 - (1 - t) / eps; };
        // dirichlet_bc(vx, boundary::left,   trial.U1x, trial.U1y, drop);
        // dirichlet_bc(vx, boundary::right,  trial.U1x, trial.U1y, drop);
        // dirichlet_bc(vx, boundary::top,    trial.U1x, trial.U1y, 1.0);
        // dirichlet_bc(vx, boundary::bottom, trial.U1x, trial.U1y, 0);
        // zero_bc(vy, trial.U2x, trial.U2y);


        // Weak BC
        // vx = 0 on left/right edge
        // for (auto i = 0; i < trial.U1y.dofs(); ++ i) {
        //     vx(0, i) = 0;
        //     vx(trial.U1x.dofs() - 1, i) = 0;
        // }

        // vy = 0 on top/bottom edge
        // for (auto i = 0; i < trial.U2x.dofs(); ++ i) {
        //     vy(i, 0) = 0;
        //     vy(i, trial.U2y.dofs() - 1) = 0;
        // }

        // int i = linear_index({0, 0}, trial.Px, trial.Py);
        // p(i, i) = 0; // fix pressure at a point
    }

    void apply_bc_test(vector_view& Rvx, vector_view& Rvy, vector_view& Rp) {
        // zero_bc(Rvx, test.U1x, test.U1y);
        // zero_bc(Rvy, test.U2y, test.U2y);
    }

    void step(int /*iter*/, double /*t*/) override {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dP = trial.Px.dofs() * trial.Py.dofs();
        auto dim_trial = dU1 + dU2 + dP;

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();
        auto dim_test = DU1 + DU2 + DP;

        std::vector<double> rhs(dim_test + dim_trial);

        vector_view Rvx{rhs.data(), {test.U1x.dofs(), test.U1y.dofs()}};
        vector_view Rvy{Rvx.data() + DU1, {test.U2x.dofs(), test.U2y.dofs()}};
        vector_view Rp{Rvy.data() + DU2, {test.Px.dofs(), test.Py.dofs()}};

        vector_view vx{Rp.data() + DP, {trial.U1x.dofs(), trial.U1y.dofs()}};
        vector_view vy{vx.data() + dU1, {trial.U2x.dofs(), trial.U2y.dofs()}};
        vector_view p{vy.data() + dU2, {trial.Px.dofs(), trial.Py.dofs()}};


        mumps::problem problem(rhs.data(), rhs.size());

        std::cout << "Assembling matrix" << std::endl;
        assembly_timer.start();
        assemble_matrix(problem);
        assembly_timer.stop();

        std::cout << "Assembling no jumps matrix" << std::endl;
        mumps::problem problem2(nullptr, dim_test);
        assemble_gram_matrix_no_jumps(problem2);
        solver.save_to_file(problem2, "no-jumps");

        std::cout << "Computing RHS" << std::endl;
        rhs_timer.start();
        compute_rhs(Rvx, Rvy, Rp);
        apply_bc(vx, vy, p);
        apply_bc_test(Rvx, Rvy, Rp);
        rhs_timer.stop();

        std::cout << "Solving" << std::endl;
        solver_timer.start();
        solver.solve(problem, "problem");
        solver_timer.stop();

        std::cout << "  matrix time:       " << static_cast<double>(assembly_timer.get()) << " ms" << std::endl;
        std::cout << "  RHS time:          " << static_cast<double>(rhs_timer.get()) << " ms" << std::endl;
        std::cout << "  solver time:       " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
        std::cout << "  assembly    FLOPS: " << solver.flops_assembly() << std::endl;
        std::cout << "  elimination FLOPS: " << solver.flops_elimination() << std::endl;

        // std::cout << "Reading the solution" << std::endl;
        // read_vector(vx.data(), dim_trial, "solution.mtx");

        auto p_avg = correct_pressure(p);
        std::cout << "Avg pressure (pre-correction): " << p_avg << std::endl;

        std::cout << "Error:" << std::endl;
        print_error(vx, vy, p, Rvx, Rvy, Rp);

        std::cout << "Outputting" << std::endl;
        outputP.to_file(p, "pressure.data");
        outputU1.to_file(vx, "vx.data");
        outputU2.to_file(vy, "vy.data");

        // save_vector(Rvx.data(), dim_test, "residuum.rhs");
        // save_vector(vx.data(), dim_trial, "solution.rhs");
    }

    void save_vector(const double* data, int size, const std::string& path) const {
        std::ofstream os{path};

        os << "%%MatrixMarket matrix array real general" << std::endl;
        os << size << " " << 1 << std::endl;

        for (int i = 0; i < size; ++ i) {
            os << data[i] << std::endl;
        }
    }

    void read_vector(double* data, int size, const std::string& path) const {
        std::ifstream is{path};

        // ignore mtx header
        for (int i = 0; i < 3; ++ i) {
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        for (int i = 0; i < size; ++ i) {
            is >> data[i];
        }
    }

    template <typename Sol>
    double correct_pressure(Sol& pressure) const {
        auto p_avg = average_value(pressure, trial.Px, trial.Py);
        for (auto i : dofs(trial.Px, trial.Py)) {
            pressure(i[0], i[1]) -= p_avg;
        }
        return p_avg;
    }

    template <typename Sol>
    double average_value(const Sol& u, const dimension& Ux, const dimension& Uy) const {
        double val = 0;
        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q, Ux, Uy);
                value_type uu = eval(u, e, q, Ux, Uy);
                val += uu.val * w * J;
            }
        }
        return val;
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) const {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

    template <typename Sol, typename Fun>
    double div_errorL2(const Sol& u, const Sol& v, const space_set& space, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return div_error(u, v, space, L2, fun);
    }

    template <typename Sol, typename Fun>
    double div_errorH1(const Sol& u, const Sol& v, const space_set& space, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy; };
        return div_error(u, v, space, H1, fun);
    }

    template <typename Sol, typename Fun, typename Norm>
    double div_error(const Sol& u, const Sol& v, const space_set& space, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(space.Px, space.Px)) {
            double J = jacobian(e, space.Px, space.Py);
            for (auto q : quad_points(space.Px, space.Py)) {
                double w = weigth(q, space.Px, space.Py);
                auto x = point(e, q, space.Px, space.Py);
                value_type div = divergence(u, v, e, q, space);

                auto d = div - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol>
    value_type divergence(const Sol& u, const Sol& v, index_type e, index_type q, const space_set& space) const {
        value_type div{};
        for (auto b : dofs_on_element(e, space.U1x, space.U1y)) {
            double c = u(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, space.U1x, space.U1y);

            const auto& bx = space.U1x.basis;
            const auto& by = space.U1y.basis;

            double B2  = by.b[e[1]][q[1]][0][loc[1]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];

            double dx = dB1 *  B2;
            double dxx = ddB1 * B2;
            double dxy = dB1 * dB2;

            div.val += c * dx;
            div.dx += c * dxx;
            div.dy += c * dxy;
        }
        for (auto b : dofs_on_element(e, space.U2x, space.U2y)) {
            double d = v(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, space.U2x, space.U2y);

            const auto& bx = space.U2x.basis;
            const auto& by = space.U2y.basis;

            double B1  = bx.b[e[0]][q[0]][0][loc[0]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB2 = by.b[e[1]][q[1]][2][loc[1]];

            double dy =  B1 * dB2;
            double dyy = B1 * ddB2;
            double dxy = dB1 * dB2;

            div.val += d * dy;
            div.dx += d * dxy;
            div.dy += d * dyy;
        }
        return div;
    }

    double dot(value_type a, point_type b) const {
        return a.dx * b[0] + a.dy * b[1];
    }

    value_type average(const value_pair& v) const {
        return 0.5 * (v.v1 + v.v2);
    }

    value_type jump(const value_pair& v) const {
        if (! v.boundary) {
            return v.v1 - v.v2;
        } else {
            return v.v1;
        }
    }

    point_type normal(boundary side) const {
        switch (side) {
        case boundary::left:   return {-1,  0};
        case boundary::right:  return { 1,  0};
        case boundary::bottom: return { 0, -1};
        case boundary::top:    return { 0,  1};
        default: return {0, 0};
        }
    }

    value_type eval_basis_at(point_type p, index_type span, index_type dof,
                             const dimension& x, const dimension& y) const {
        bspline::eval_ders_ctx cx{x.p, 1};
        bspline::eval_ders_ctx cy{y.p, 1};

        double** bvx = cx.basis_vals();
        double** bvy = cy.basis_vals();

        eval_basis_with_derivatives(span[0], p[0], x.B, bvx, 1, cx);
        eval_basis_with_derivatives(span[1], p[1], y.B, bvy, 1, cy);

        int offsetx = span[0] - x.p;
        int offsety = span[1] - y.p;

        int ix = dof[0] - offsetx;
        int iy = dof[1] - offsety;

        if (ix < 0 || ix > x.p || iy < 0 || iy > y.p) return {};

        auto value = bvx[0][ix] * bvy[0][iy];
        auto dx    = bvx[1][ix] * bvy[0][iy];
        auto dy    = bvx[0][ix] * bvy[1][iy];

        return { value, dx, dy };
    }

    value_type eval_basis_at(point_type p, index_type dof, const dimension& x, const dimension& y) const {
        int spanx = bspline::find_span(p[0], x.B);
        int spany = bspline::find_span(p[1], y.B);

        return eval_basis_at(p, {spanx, spany}, dof, x, y);
    }

    value_pair eval_basis_at_edge(point_type p, boundary orientation, index_type dof,
                                  const dimension& x, const dimension& y) const {
        int spanx = bspline::find_span(p[0], x.B);
        int spany = bspline::find_span(p[1], y.B);
        auto span = index_type{spanx, spany};
        auto val1 = eval_basis_at(p, span, dof, x, y);

        if (orientation == boundary::vertical) {
            int spanx0 = bspline::find_span(p[0] - 1e-10, x.B);
            auto val0 = eval_basis_at(p, index_type{spanx0, spany}, dof, x, y);
            return {val0, val1, spanx == spanx0};
        } else {
            int spany0 = bspline::find_span(p[1] - 1e-10, y.B);
            auto val0 = eval_basis_at(p, index_type{spanx, spany0}, dof, x, y);
            return {val0, val1, spany == spany0};
        }
    }

    bool supported_in_1d(int dof, int e, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return e >= xrange.first && e <= xrange.second;
    }

    bool touches_point(int dof, int sx, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return sx >= xrange.first && sx <= xrange.second + 1;
    }

    template <typename Form>
    double integrate_vface(int sx, int ey, index_type i, index_type j,
                           const dimension& Ux, const dimension& Uy,
                           const dimension& Vx, const dimension& Vy,
                           Form&& form) const {

        if (! touches_point(i[0], sx, Ux)   || ! touches_point(j[0], sx, Vx))   return 0;
        if (! supported_in_1d(i[1], ey, Uy) || ! supported_in_1d(j[1], ey, Vy)) return 0;

        double x0 = Ux.basis.points[sx];
        auto n = point_type{1, 0};

        double val = 0;
        for (int q = 0; q < Uy.basis.quad_order; ++ q) {
            double w = Uy.basis.w[q];
            double J = Uy.basis.J[ey];

            auto x = point_type{x0, Uy.basis.x[ey][q]};
            auto uu = eval_basis_at_edge(x, boundary::vertical, i, Ux, Uy);
            auto vv = eval_basis_at_edge(x, boundary::vertical, j, Vx, Vy);

            auto fuv = form(uu, vv, x, n);
            val += fuv * w * J;
        }
        return val;
    }

    template <typename Form>
    double integrate_hface(int sy, int ex, index_type i, index_type j,
                           const dimension& Ux, const dimension& Uy,
                           const dimension& Vx, const dimension& Vy,
                           Form&& form) const {

        if (! supported_in_1d(i[0], ex, Ux) || ! supported_in_1d(j[0], ex, Vx)) return 0;
        if (! touches_point(i[1], sy, Uy)   || ! touches_point(j[1], sy, Vy))   return 0;

        double y0 = Uy.basis.points[sy];
        auto n = point_type{0, 1};

        double val = 0;
        for (int q = 0; q < Ux.basis.quad_order; ++ q) {
            double w = Ux.basis.w[q];
            double J = Ux.basis.J[ex];

            auto x = point_type{Ux.basis.x[ex][q], y0};
            auto uu = eval_basis_at_edge(x, boundary::horizontal, i, Ux, Uy);
            auto vv = eval_basis_at_edge(x, boundary::horizontal, j, Vx, Vy);

            auto fuv = form(uu, vv, x, n);
            val += fuv * w * J;
        }
        return val;
    }

    template <typename Form>
    double integrate_over_skeleton(index_type i, index_type j,
                                   const dimension& Ux, const dimension& Uy,
                                   const dimension& Vx, const dimension& Vy,
                                   Form&& form) const {
        using std::min;
        using std::max;

        auto rUx = Ux.basis.element_ranges[i[0]];
        auto rVx = Vx.basis.element_ranges[j[0]];
        auto rUy = Uy.basis.element_ranges[i[1]];
        auto rVy = Vy.basis.element_ranges[j[1]];
        auto ex0 = min(rUx.first, rVx.first);
        auto ex1 = max(rUx.second, rVx.second);
        auto ey0 = min(rUy.first, rVy.first);
        auto ey1 = max(rUy.second, rVy.second);

        double val = 0;
        for (int ix = ex0; ix <= ex1 + 1; ++ ix) {
            for (auto ey = ey0; ey <= ey1; ++ ey) {
                val += integrate_vface(ix, ey, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        for (auto ex = ex0; ex <= ex1; ++ ex) {
            for (int iy = ey0; iy <= ey1 + 1; ++ iy) {
                val += integrate_hface(iy, ex, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        return val;
    }

    template <typename Form>
    double integrate_over_internal_skeleton(index_type i, index_type j,
                                            const dimension& Ux, const dimension& Uy,
                                            const dimension& Vx, const dimension& Vy,
                                            Form&& form) const {
        using std::min;
        using std::max;

        auto rUx = Ux.basis.element_ranges[i[0]];
        auto rVx = Vx.basis.element_ranges[j[0]];
        auto rUy = Uy.basis.element_ranges[i[1]];
        auto rVy = Vy.basis.element_ranges[j[1]];
        auto ex0 = min(rUx.first, rVx.first);
        auto ex1 = max(rUx.second, rVx.second);
        auto ey0 = min(rUy.first, rVy.first);
        auto ey1 = max(rUy.second, rVy.second);

        double val = 0;
        for (int ix = max(1, ex0); ix < min(Ux.elements, ex1 + 2); ++ ix) {
            for (auto ey = ey0; ey <= ey1; ++ ey) {
                val += integrate_vface(ix, ey, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        for (auto ex = ex0; ex <= ex1; ++ ex) {
            for (int iy = max(1, ey0); iy < min(Uy.elements, ey1 + 2); ++ iy) {
                val += integrate_hface(iy, ex, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        return val;
    }


    template <typename Form>
    double integrate_boundary(boundary side, index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                              const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Uy.elements - 1;
            if (! supported_in_1d(j[1], ey, Vy) || ! supported_in_1d(i[1], ey, Uy)) return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                if (! supported_in_1d(j[0], e, Vx)) continue;

                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++ q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type ww = eval_basis_at(x, {e, ey}, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, {e, ey}, j, Vx, Vy);
                    double fuw = form(ww, uu, x);
                    val += fuw * w * J;
                }
            }
        } else {
            int ex = side == boundary::left ? 0 : Ux.elements - 1;
            if (! supported_in_1d(j[0], ex, Vx) || ! supported_in_1d(i[0], ex, Ux)) return 0;

            auto x0 = side == boundary::left ? Ux.a : Ux.b;

            for (auto e : Uy.basis.element_range(i[1])) {
                if (! supported_in_1d(j[1], e, Vy)) continue;

                double J = Uy.basis.J[e];

                for (int q = 0; q < Uy.basis.quad_order; ++ q) {
                    double w = Uy.basis.w[q];
                    point_type x{x0, Uy.basis.x[e][q]};
                    value_type ww = eval_basis_at(x, {ex, e}, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, {ex, e}, j, Vx, Vy);
                    double fuw = form(ww, uu, x);
                    val += fuw * w * J;
                }
            }
        }
        return val;
    }

    template <typename Form>
    double integrate_boundary(boundary side, index_type i, const dimension& Ux, const dimension& Uy, Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;
        auto nv = normal(side);

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Uy.elements - 1;
            if (! supported_in_1d(i[1], ey, Uy)) return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++ q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    double fuw = form(ww, x, nv);
                    val += fuw * w * J;
                }
            }
        } else {
            int ex = side == boundary::left ? 0 : Ux.elements - 1;
            if (! supported_in_1d(i[0], ex, Ux)) return 0;

            auto x0 = side == boundary::left ? Ux.a : Ux.b;

            for (auto e : Uy.basis.element_range(i[1])) {
                double J = Uy.basis.J[e];

                for (int q = 0; q < Uy.basis.quad_order; ++ q) {
                    double w = Uy.basis.w[q];
                    point_type x{x0, Uy.basis.x[e][q]};
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    double fuw = form(ww, x, nv);
                    val += fuw * w * J;
                }
            }
        }
        return val;
    }


    bool touches(index_type dof, boundary side, const dimension& x, const dimension& y) const {
        if (side == boundary::left || side == boundary::right) {
            auto e = side == boundary::left ? 0 : x.elements - 1;
            return supported_in_1d(dof[0], e, x);
        } else {
            auto e = side == boundary::bottom ? 0 : y.elements - 1;
            return supported_in_1d(dof[1], e, y);
        }
    }


    template <typename Sol>
    double normV(const Sol& vx, const Sol& vy, const Sol& p) const {
        double norm = 0;

        // w, p -> (w, p) +  hh ([w], [p])
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : touching_dofs(i, test.Px, test.Py)) {
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };
                double val = eval([](auto w, auto p) { return w.val * p.val; });

                auto form = [this](auto w, auto p, auto, auto) {
                    return h * jump(w).val * jump(p).val;
                };
                val += integrate_over_internal_skeleton(i, j, test.Px, test.Py, test.Px, test.Py, form);
                norm += val * p(i[0], i[1]) * p(j[0], j[1]);
            }
        }

        // tx, vx -> (\/tx, \/vx) + 1/h ([tx], [vx])
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : touching_dofs(i, test.U1x, test.U1y)) {
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                double val = eval([this](auto tx, auto vx) { return grad_dot(tx, vx); });
                // skeleton
                auto form = [this](auto tx, auto vx, auto, auto) {
                    return eta/h * jump(tx).val * jump(vx).val;
                };
                val += integrate_over_skeleton(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form);
                norm += val * vx(i[0], i[1]) * vx(j[0], j[1]);
            }
        }

        // ty, vy -> (\/ty, \/vy) + 1/h ([ty], [vy])
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : touching_dofs(i, test.U2x, test.U2y)) {
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                double val = eval([this](auto ty, auto vy) { return grad_dot(ty, vy); });
                // skeleton
                auto form = [this](auto ty, auto vy, auto, auto) {
                    return eta/h * jump(ty).val * jump(vy).val;
                };
                val += integrate_over_skeleton(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form);
                norm += val * vy(i[0], i[1]) * vy(j[0], j[1]);
            }
        }

        return std::sqrt(norm);
    }
};

}

#endif // PROBLEMS_STOKES_STOKES_HPP_
