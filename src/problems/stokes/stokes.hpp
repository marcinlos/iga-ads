#ifndef PROBLEMS_STOKES_STOKES_HPP_
#define PROBLEMS_STOKES_STOKES_HPP_

#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "problems/stokes/space_set.hpp"
#include "mumps.hpp"

#include <galois/Timer.h>


namespace ads {

class stokes : public simulation_2d {
private:
    using Base = simulation_2d;

    galois_executor executor{8};

    space_set trial, test, ref;

    double h;

    mumps::solver solver;
    galois::StatTimer solver_timer{"solver"};
    output_manager<2> outputU1, outputU2, outputP;

public:
    stokes(space_set trial, space_set test, space_set ref, const timesteps_config& steps)
    : Base{ trial.U1x, trial.U1y, steps }
    , trial{ std::move(trial) }
    , test{ std::move(test) }
    , ref{ std::move(ref) }
    , h{ element_diam(this->trial.Px, this->trial.Py) }
    , outputU1{ this->trial.U1x.B, this->trial.U1y.B, 500 }
    , outputU2{ this->trial.U2x.B, this->trial.U2y.B, 500 }
    , outputP{ this->trial.Px.B, this->trial.Py.B, 500 }
    { }

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

        output_exact();
    }

    value_type exact_p(point_type p) const {
        auto x = p[0];
        return {x * (1 - x), 1 - 2 * x, 0.0};
    }

    std::array<value_type, 2> exact_v(point_type p) const {
        auto f = [](double x, double y) {
            return x*x * (1 - x) * (1 - x) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        auto dfx = [](double x, double y) {
            return (4 * x*x*x - 6 * x*x + 2 * x) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        auto dfy = [](double x, double y) {
            return x*x * (1 - x) * (1 - x) * (2 - 12 * y + 12 * y*y);
        };

        double x = p[0], y = p[1];
        value_type vx = {f(x, y), dfx(x, y), dfy(x, y)};
        value_type vy = {-f(y, x), -dfy(y, x), -dfx(y, x)};

        return { vx ,vy };
    }

    value_type exact_div(point_type p) const {
        auto v = exact_v(p);
        auto div = v[0].dx + v[1].dy;

        auto dfxy = [](double x, double y) {
            return (4 * x*x*x - 6 * x*x + 2 * x) * (2 - 12 * y + 12 * y*y);
        };

        auto dfxx = [](double x, double y) {
            return (12 * x*x - 12 * x + 2) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        double x = p[0], y = p[1];
        double dx = dfxx(x, y) - dfxy(y, x);
        double dy = dfxy(x, y) - dfxx(y, x);

        return { div, dx, dy };
    }

    void output_exact() {
        auto p = [this](point_type x) { return exact_p(x).val; };
        auto vx = [this](point_type x) { return exact_v(x)[0].val; };
        auto vy = [this](point_type x) { return exact_v(x)[1].val; };

        auto Ux = trial.U1x;
        auto Uy = trial.U1y;
        auto buffer = vector_type{{ Ux.dofs(), Uy.dofs() }};

        auto project = [&](auto fun) {
            vector_type rhs{{ Ux.dofs(), Uy.dofs() }};
            compute_projection(rhs, Ux.basis, Uy.basis, [&](double x, double y) { return fun({x, y}); });
            ads_solve(rhs, buffer, Ux.data(), Uy.data());
            return rhs;
        };

        outputU1.to_file(project(p), "pressure_ref.data");
        outputU2.to_file(project(vx), "vx_ref.data");
        outputP.to_file(project(vy), "vy_ref.data");
    }

    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p) const {
        auto e_vx = [this](point_type x) { return exact_v(x)[0]; };
        auto e_vy = [this](point_type x) { return exact_v(x)[1]; };
        auto e_p = [this](point_type x) { return exact_p(x); };
        auto div = [this](point_type x) { return exact_div(x); };

        double vxL2 = error_relative_L2(vx, trial.U1x, trial.U1y, e_vx) * 100;
        double vxH1 = error_relative_H1(vx, trial.U1x, trial.U1y, e_vx) * 100;

        double vyL2 = error_relative_L2(vy, trial.U2x, trial.U2y, e_vy) * 100;
        double vyH1 = error_relative_H1(vy, trial.U2x, trial.U2y, e_vy) * 100;

        double vL2 = std::hypot(vxL2, vyL2) / std::sqrt(2);
        double vH1 = std::hypot(vxH1, vyH1) / std::sqrt(2);

        double pL2 = error_relative_L2(p, trial.Px, trial.Py, e_p) * 100;
        double pH1 = error_relative_H1(p, trial.Px, trial.Py, e_p) * 100;

        double divL2 = div_errorL2(vx, vy, trial.U1x, trial.U1y, div) * 100;
        double divH1 = div_errorH1(vx, vy, trial.U1x, trial.U1y, div) * 100;

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
        std::cout << "v   : L2 = " << vL2    << "%, H1 = " << vH1    << "%" << std::endl;
        std::cout << "p   : L2 = " << pL2    << "%, H1 = " << pH1    << "%" << std::endl;
        std::cout << "div : L2 = " << divL2  << ", H1 = " << divH1  << std::endl;
    }

    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p,
                     const vector_type& ref_vx, const vector_type& ref_vy, const vector_type& ref_p) const {

        auto cU1x = bspline::eval_ders_ctx{ref.U1x.p, 1};
        auto cU1y = bspline::eval_ders_ctx{ref.U1y.p, 1};
        auto cU2x = bspline::eval_ders_ctx{ref.U2x.p, 1};
        auto cU2y = bspline::eval_ders_ctx{ref.U2y.p, 1};
        auto cPx = bspline::eval_ders_ctx{ref.Px.p, 1};
        auto cPy = bspline::eval_ders_ctx{ref.Py.p, 1};

        auto e_vx = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_vx, ref.U1x.B, ref.U1y.B, cU1x, cU1y); };
        auto e_vy = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_vy, ref.U2x.B, ref.U2y.B, cU2x, cU2y); };
        auto e_p = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_p, ref.Px.B, ref.Py.B, cPx, cPy); };
        auto div = [this](point_type x) { return value_type{}; };

        double vxL2 = error_relative_L2(vx, trial.U1x, trial.U1y, e_vx) * 100;
        double vxH1 = error_relative_H1(vx, trial.U1x, trial.U1y, e_vx) * 100;

        double vyL2 = error_relative_L2(vy, trial.U2x, trial.U2y, e_vy) * 100;
        double vyH1 = error_relative_H1(vy, trial.U2x, trial.U2y, e_vy) * 100;

        double vL2 = std::hypot(vxL2, vyL2) / std::sqrt(2);
        double vH1 = std::hypot(vxH1, vyH1) / std::sqrt(2);

        double pL2 = error_relative_L2(p, trial.Px, trial.Py, e_p) * 100;
        double pH1 = error_relative_H1(p, trial.Px, trial.Py, e_p) * 100;

        double divL2 = div_errorL2(vx, vy, trial.Px, trial.Py, div) * 100;
        double divH1 = div_errorH1(vx, vy, trial.Px, trial.Py, div) * 100;

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
        std::cout << "v   : L2 = " << vL2    << "%, H1 = " << vH1    << "%" << std::endl;
        std::cout << "p   : L2 = " << pL2    << "%, H1 = " << pH1    << "%" << std::endl;
        std::cout << "div : L2 = " << divL2  << ", H1 = " << divH1  << std::endl;
    }



    point_type forcing(point_type p) const {
        double x = p[0], y = p[1];

        auto fx =
            (12 - 24 * y) * x*x*x*x +
            (-24 + 48 * y) * x*x*x +
            (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
            (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
            1 - 4 * y + 12 * y*y - 8 * y*y*y;

        auto fy =
            (8 - 48 * y + 48 * y*y) * x*x*x +
            (-12 + 72 * y - 72 * y*y) * x*x +
            (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x -
            12 * y*y + 24 * y*y*y - 12 * y*y*y*y;

        return { fx, fy };
        // return { 0, 0 };
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

    bool is_pressure_fixed(index_type dof) const {
        return dof[0] == 0 && dof[1] == 0;
    }

    void assemble_matrix(mumps::problem& problem) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();

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

        auto hh = h * h;

        // Gram matrix
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : overlapping_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                if (! is_boundary(i, test.U1x, test.U1y) && ! is_boundary(j, test.U1x, test.U1y)) {
                    double value = eval([](auto u, auto v) { return u.dx * v.dx + u.dy * v.dy; });
                    test_vx(ii, jj, value);
                }
            }
        }

        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : overlapping_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                if (! is_boundary(i, test.U2x, test.U2y) && ! is_boundary(j, test.U2x, test.U2y)) {
                    double value = eval([](auto u, auto v) { return u.dx * v.dx + u.dy * v.dy; });
                    test_vy(ii, jj, value);
                }
            }
        }

        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : overlapping_dofs(i, test.Px, test.Py)) {
                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, test.Px, test.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };

                if (! is_pressure_fixed(i) && ! is_pressure_fixed(j)) {
                    double value = eval([](auto p, auto q) { return p.val * q.val; });
                    test_p(ii, jj, value);
                }
            }
        }

        // Dirichlet BC - test space
        for_boundary_dofs(test.U1x, test.U1y, [&](index_type dof) {
            int i = linear_index(dof, test.U1x, test.U1y) + 1;
            test_vx(i, i, 1);
        });
        for_boundary_dofs(test.U2x, test.U2y, [&](index_type dof) {
            int i = linear_index(dof, test.U2x, test.U2y) + 1;
            test_vy(i, i, 1);
        });
        test_p(1, 1, 1.0);

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

        // B(v, u)
        // u = (ux, uy, p)
        // v = (vx, vy, q)

        // vx, ux -> (\/vx, \/ux)
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! overlap(i, test.U1x, test.U1y, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;

                bool bd_i = is_boundary(i, test.U1x, test.U1y);
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y);

                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form); };
                double value = eval([](auto v, auto u) { return v.dx * u.dx + v.dy * u.dy; });

                put(ii, jj, 0, 0, value, bd_i, bd_j);
            }
        }
        // vy, uy -> (\/vy, \/uy)
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! overlap(i, test.U2x, test.U2y, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;

                bool bd_i = is_boundary(i, test.U2x, test.U2y);
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y);

                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form); };
                double value = eval([](auto v, auto u) { return v.dx * u.dx + v.dy * u.dy; });

                put(ii, jj, DU1, dU1, value, bd_i, bd_j);
            }
        }
        // q, ux -> (q, ux,x)
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! overlap(i, test.Px, test.Py, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;

                bool fixed_i = is_pressure_fixed(i);
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y);

                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U1x, trial.U1y, form); };
                double value = eval([](auto q, auto u) { return q.val * u.dx; });
                put(ii, jj, DU1 + DU2, 0, value, fixed_i, bd_j);
            }
        }
        // q, uy -> (q, uy,y)
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! overlap(i, test.Px, test.Py, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;

                bool fixed_i = is_pressure_fixed(i);
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y);

                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U2x, trial.U2y, form); };
                double value = eval([](auto q, auto u) { return q.val * u.dy; });
                put(ii, jj, DU1 + DU2, dU1, value, fixed_i, bd_j);
            }
        }
        // vx, p -> - (vx,x, p)
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! overlap(i, test.U1x, test.U1y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;

                bool bd_i = is_boundary(i, test.U1x, test.U1y);
                bool fixed_j = is_pressure_fixed(j);

                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.Px, trial.Py, form); };
                double value = eval([](auto v, auto p) { return -v.dx * p.val; });
                put(ii, jj, 0, dU1 + dU2, value, bd_i, fixed_j);
            }
        }
        // vy, p ->  - (vy,y, p)
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! overlap(i, test.U2x, test.U2y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;

                bool bd_i = is_boundary(i, test.U2x, test.U2y);
                bool fixed_j = is_pressure_fixed(j);

                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.Px, trial.Py, form); };
                double value = eval([](auto v, auto p) { return -v.dy * p.val; });
                put(ii, jj, DU1, dU1 + dU2, value, bd_i, fixed_j);
            }
        }

        // Dirichlet BC - trial space
        for_boundary_dofs(trial.U1x, trial.U1y, [&](index_type dof) {
            int i = linear_index(dof, trial.U1x, trial.U1y) + 1;
            trial_vx(i, i, 1);
        });
        for_boundary_dofs(trial.U2x, trial.U2y, [&](index_type dof) {
            int i = linear_index(dof, trial.U2x, trial.U2y) + 1;
            trial_vy(i, i, 1);
        });
        trial_p(1, 1, 1.0);
    }

    void compute_rhs(vector_view& vx, vector_view& vy, vector_view& /*p*/) const {
        using shape = std::array<std::size_t, 2>;
        auto u1_shape = shape{ test.U1x.basis.dofs_per_element(), test.U1y.basis.dofs_per_element() };
        auto u2_shape = shape{ test.U2x.basis.dofs_per_element(), test.U2y.basis.dofs_per_element() };


        executor.for_each(elements(trial.Px, trial.Py), [&](index_type e) {
            auto vx_loc = vector_type{ u1_shape };
            auto vy_loc = vector_type{ u2_shape };

            double J = jacobian(e);
            for (auto q : quad_points(trial.Px, trial.Py)) {
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
    }


    void apply_bc(vector_view& Rvx, vector_view& Rvy, vector_view& Rp) {
        zero_bc(Rvx, trial.U1x, trial.U1y);

        // constexpr double eps = 1e-4;
        // auto drop = [](double t) { return t < 1 - eps ? 0 : 1 - (1 - t) / eps; };
        // dirichlet_bc(Rvx, boundary::left,   trial.U1x, trial.U1y, drop);
        // dirichlet_bc(Rvx, boundary::right,  trial.U1x, trial.U1y, drop);
        // dirichlet_bc(Rvx, boundary::top,    trial.U1x, trial.U1y, 1);
        // dirichlet_bc(Rvx, boundary::bottom, trial.U1x, trial.U1y, 0);

        zero_bc(Rvy, trial.U2x, trial.U2y);

        Rp(0, 0) = 0; // fix pressure at a point
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
        assemble_matrix(problem);

        std::cout << "Computing RHS" << std::endl;
        compute_rhs(Rvx, Rvy, Rp);
        apply_bc(vx, vy, p);


        std::cout << "Solving" << std::endl;
        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        std::cout << "  solver time:       " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
        std::cout << "  assembly    FLOPS: " << solver.flops_assembly() << std::endl;
        std::cout << "  elimination FLOPS: " << solver.flops_elimination() << std::endl;

        std::cout << "Error:" << std::endl;
        print_error(vx, vy, p);

        std::cout << "Outputting" << std::endl;
        outputU1.to_file(p, "pressure.data");
        outputU2.to_file(vx, "vx.data");
        outputP.to_file(vy, "vy.data");

        // print_solution("vx.sol", vx, Ux, Uy);
        // print_solution("vy.sol", vy, Ux, Uy);
        // print_solution("p.sol", p, Ux, Uy);

        // auto ref_vx = read_solution("vx.sol", ref_x, ref_y);
        // auto ref_vy = read_solution("vy.sol", ref_x, ref_y);
        // auto ref_p = read_solution("p.sol", ref_x, ref_y);

        // print_error(vx, vy, p, ref_vx, ref_vy, ref_p);
    }

    vector_type read_solution(const char* filename, const dimension& Ux, const dimension& Uy) const {
        auto sol = vector_type{{ Ux.dofs(), Uy.dofs() }};
        std::ifstream in{filename};
        while (in) {
            index_type dof;
            in >> dof[0] >> dof[1];
            in >> sol(dof[0], dof[1]);
        }
        return sol;
    }

    template <typename V>
    void print_solution(const char* filename, const V& u, const dimension& Ux, const dimension& Uy) const {
        std::ofstream out{filename};
        for (auto dof : dofs(Ux, Uy)) {
            out << dof[0] << " " << dof[1] << " " << std::setprecision(16) << u(dof[0], dof[1]) << std::endl;
        }
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) const {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

    template <typename Sol, typename Fun>
    double div_errorL2(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return div_error(u, v, Ux, Uy, L2, fun);
    }

    template <typename Sol, typename Fun>
    double div_errorH1(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy; };
        return div_error(u, v, Ux, Uy, H1, fun);
    }

    template <typename Sol, typename Fun, typename Norm>
    double div_error(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q, Ux, Uy);
                auto x = point(e, q, Ux, Uy);
                value_type div = divergence(u, v, e, q, Ux, Uy);

                auto d = div - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol>
    value_type divergence(const Sol& u, const Sol& v, index_type e, index_type q,
                          const dimension& x, const dimension& y) const {
        value_type div{};
        for (auto b : dofs_on_element(e, x, y)) {
            double c = u(b[0], b[1]);
            double d = v(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, x, y);

            const auto& bx = x.basis;
            const auto& by = y.basis;

            double B1  = bx.b[e[0]][q[0]][0][loc[0]];
            double B2  = by.b[e[1]][q[1]][0][loc[1]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];
            double ddB2 = by.b[e[1]][q[1]][2][loc[1]];

            double dx = dB1 *  B2;
            double dy =  B1 * dB2;
            double dxx = ddB1 * B2;
            double dyy = B1 * ddB2;
            double dxy = dB1 * dB2;

            div.val += c * dx + d * dy;
            div.dx += c * dxx + d * dxy;
            div.dy += c * dxy + d * dyy;
        }
        return div;
    }

    template <typename Form>
    double integrate(index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                     const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;

        for (auto e : elements_supporting_dof(i, Ux, Uy)) {
            if (! supported_in(j, e, Vx, Vy)) continue;

            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q, Ux, Uy);
                value_type ww = eval_basis(e, q, i, Ux, Uy);
                value_type uu = eval_basis(e, q, j, Vx, Vy);

                double fuw = form(ww, uu);
                val += fuw * w * J;
            }
        }
        return val;
    }

};

}

#endif // PROBLEMS_STOKES_STOKES_HPP_
