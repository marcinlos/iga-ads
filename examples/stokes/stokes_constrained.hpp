// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef STOKES_STOKES_CONSTRAINED_HPP
#define STOKES_STOKES_CONSTRAINED_HPP

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/solver/mumps.hpp"
#include "space_set.hpp"


namespace ads {

class stokes_constrained : public simulation_2d {
private:
    using Base = simulation_2d;

    galois_executor executor{8};

    space_set trial, test;

    double h;

    double Cpen;// = 10;//5 * (3 + 1);
    double hF;// = 1 / 20.;

    mumps::solver solver;
    output_manager<2> outputU1, outputU2, outputP;

public:
    stokes_constrained(space_set trial_, space_set test_, const timesteps_config& steps)
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
        hF = 1. / trial.Px.dofs();
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

    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p) const {
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

        double divL2 = div_errorL2(vx, vy, trial, div);
        double divH1 = div_errorH1(vx, vy, trial, div);

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
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

    template <typename Form>
    double integrate(index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                     const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;

        for (auto e : elements_supporting_dof(i, Ux, Uy)) {
            if (! supported_in(j, e, Vx, Vy)) continue;

            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q);
                value_type ww = eval_basis(e, q, i, Ux, Uy);
                value_type uu = eval_basis(e, q, j, Vx, Vy);

                double fuw = form(ww, uu);
                val += fuw * w * J;
            }
        }
        return val;
    }

    value_type eval_basis_at(point_type p, index_type e, index_type dof, const dimension& x, const dimension& y) const {
        int spanx = bspline::find_span(p[0], x.B);
        int spany = bspline::find_span(p[1], y.B);

        bspline::eval_ders_ctx cx{x.p, 1};
        bspline::eval_ders_ctx cy{y.p, 1};

        double** bvx = cx.basis_vals();
        double** bvy = cy.basis_vals();

        eval_basis_with_derivatives(spanx, p[0], x.B, bvx, 1, cx);
        eval_basis_with_derivatives(spany, p[1], y.B, bvy, 1, cy);

        int offsetx = spanx - x.p;
        int offsety = spany - y.p;

        int ix = dof[0] - offsetx;
        int iy = dof[1] - offsety;

        auto value = bvx[0][ix] * bvy[0][iy];
        auto dx    = bvx[1][ix] * bvy[0][iy];
        auto dy    = bvx[0][ix] * bvy[1][iy];

        return { value, dx, dy };
    }

    bool supported_in_1d(int dof, int e, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return e >= xrange.first && e <= xrange.second;
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

    void assemble_matrix(mumps::problem& problem) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();

        auto D = DU1 + DU2;

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DU1, DU1, problem);

        auto trial_vx = shifted(D, D, problem);
        auto trial_vy = shifted(D + dU1, D + dU1, problem);
        auto trial_p = shifted(D + dU1 + dU2, D + dU1 + dU2, problem);

        auto hh = h * h;

        // Gram matrix
        // G(w, u)
        // w = (tx, ty)
        // u = (vx, vy)

        // tx, vx -> (tx, vx) + hh (tx,x, vx,x)     ---- old
        // tx, vx -> (tx, vx) + (\/tx,x, \/vx)
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : overlapping_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                if (! is_boundary(i[0], test.U1x) && ! is_boundary(j[0], test.U1x)) {
                    test_vx(ii, jj, eval([hh](auto w, auto u) {
                        return w.val * u.val + hh * (w.dx * u.dx + w.dy * u.dy);
                    }));
                }
            }
        }

        // ty, vy -> (ty, vy) + hh (ty,y, vy,y)     ---- old
        // ty, vy -> (ty, vy) + (\/ty, \/vy)
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : overlapping_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                if (! is_boundary(i[1], test.U2y) && ! is_boundary(j[1], test.U2y)) {
                    test_vy(ii, jj, eval([hh](auto w, auto u) {
                        return w.val * u.val + hh * (w.dy * u.dy + w.dx + u.dx);
                    }));
                }
            }
        }

        // Dirichlet BC - test space
        for (auto iy = 0; iy < test.U1y.dofs(); ++ iy) {
            int i0 = linear_index({0, iy}, test.U1x, test.U1y) + 1;
            test_vx(i0, i0, 1);
            int i1 = linear_index({test.U1x.dofs() - 1, iy}, test.U1x, test.U1y) + 1;
            test_vx(i1, i1, 1);
        }
        for (auto ix = 0; ix < test.U2x.dofs(); ++ ix) {
            int i0 = linear_index({ix, 0}, test.U2x, test.U2y) + 1;
            test_vy(i0, i0, 1);
            int i1 = linear_index({ix, test.U2y.dofs() - 1}, test.U2x, test.U2y) + 1;
            test_vy(i1, i1, 1);
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        // B(w, u)
        // w = (tx, ty, w)
        // u = (vx, vy, p)

        auto pen_term = [this](auto w, auto u) { return - 2 * Cpen / hF * w.val * u.val; };

        // tx, vx -> (\/tx, \/vx)
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! overlap(i, test.U1x, test.U1y, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form); };

                bool bd_i = is_boundary(i[0], test.U1x);
                bool bd_j = is_boundary(j[0], trial.U1x);

                double value = eval([](auto w, auto u) { return w.dx * u.dx + w.dy * u.dy; });

                // boundary terms
                // top/bottom
                if (is_boundary(i[1], test.U1y) || is_boundary(j[1], trial.U1y)) {
                    auto side = i[1] == 0 ? boundary::bottom : boundary::top;
                    int n = side == boundary::top ? 1 : -1;
                    value -= integrate_boundary(side, i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, [&](auto w, auto u, auto) {
                        return n * w.dy * u.val + pen_term(w, u);
                    });
                }
                put(ii, jj, 0, 0, value, bd_i, bd_j);
            }
        }

        // ty, vy -> (\/ty, \/vy)
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! overlap(i, test.U2x, test.U2y, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form); };

                bool bd_i = is_boundary(i[1], test.U2y);
                bool bd_j = is_boundary(j[1], trial.U2y);

                double value = eval([](auto w, auto u) { return w.dx * u.dx + w.dy * u.dy; });

                // boundary terms
                // left/right
                if (is_boundary(i[0], test.U2x) || is_boundary(j[0], trial.U2x)) {
                    auto side = i[0] == 0 ? boundary::left : boundary::right;
                    int n = side == boundary::right ? 1 : -1;
                    value -= integrate_boundary(side, i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, [&](auto w, auto u, auto) {
                        return n * w.dx * u.val + pen_term(w, u);
                    });
                }
                put(ii, jj, DU1, dU1, value, bd_i, bd_j);
            }
        }

        // tx, p -> (tx, p,x) = - (tx,x, p)
        for (auto i : dofs(trial.U1x, trial.U1y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! overlap(i, trial.U1x, trial.U1y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, trial.U1x, trial.U1y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, trial.U1x, trial.U1y, trial.Px, trial.Py, form); };

                bool bd_i = is_boundary(i[0], trial.U1x);
                bool fixed_j = is_pressure_fixed(j);

                put(ii, jj, DU1 + DU2, dU1 + dU2, eval([](auto w, auto u) { return w.dx * u.val; }), bd_i, fixed_j);
            }
        }

        // ty, p -> (ty, p,y) = - (ty,y, p)
        for (auto i : dofs(trial.U2x, trial.U2y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! overlap(i, trial.U2x, trial.U2y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, trial.U2x, trial.U2y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, trial.U2x, trial.U2y, trial.Px, trial.Py, form); };

                bool bd_i = is_boundary(i[1], trial.U2y);
                bool fixed_j = is_pressure_fixed(j);

                put(ii, jj, DU1 + DU2 + dU1, dU1 + dU2, eval([](auto w, auto u) { return w.dy * u.val; }), bd_i, fixed_j);
            }
        }

        // Dirichlet BC - trial space
        for (auto iy = 0; iy < trial.U1y.dofs(); ++ iy) {
            int i0 = linear_index({0, iy}, trial.U1x, trial.U1y) + 1;
            trial_vx(i0, i0, 1);
            int i1 = linear_index({trial.U1x.dofs() - 1, iy}, trial.U1x, trial.U1y) + 1;
            trial_vx(i1, i1, 1);
        }
        for (auto ix = 0; ix < trial.U2x.dofs(); ++ ix) {
            int i0 = linear_index({ix, 0}, trial.U2x, trial.U2y) + 1;
            trial_vy(i0, i0, 1);
            int i1 = linear_index({ix, trial.U2y.dofs() - 1}, trial.U2x, trial.U2y) + 1;
            trial_vy(i1, i1, 1);
        }

        // for_boundary_dofs(trial.U1x, trial.U1y, [&](index_type dof) {
        //     int i = linear_index(dof, trial.U1x, trial.U1y) + 1;
        //     trial_vx(i, i, 1);
        // });
        // for_boundary_dofs(trial.U2x, trial.U2y, [&](index_type dof) {
        //     int i = linear_index(dof, trial.U2x, trial.U2y) + 1;
        //     trial_vy(i, i, 1);
        // });
        int ii = linear_index({0, 0}, trial.Px, trial.Py) + 1;
        trial_p(ii, ii, 1.0);

        // for (auto j : dofs(trial.Px, trial.Py)) {
        //     int jj = linear_index(j, trial.Px, trial.Py) + 1;
        //     trial_p(ii, jj, integrate({0,0}, j, trial.Px, trial.Py, trial.Px, trial.Py, [](auto w, auto u) { return u.val; }));
        // }
    }

    void compute_rhs(vector_view& vx, vector_view& vy) const {
        using shape = std::array<std::size_t, 2>;
        auto u1_shape = shape{ test.U1x.basis.dofs_per_element(), test.U1y.basis.dofs_per_element() };
        auto u2_shape = shape{ test.U2x.basis.dofs_per_element(), test.U2y.basis.dofs_per_element() };


        executor.for_each(elements(test.Px, test.Py), [&](index_type e) {
            auto vx_loc = vector_type{ u1_shape };
            auto vy_loc = vector_type{ u2_shape };

            double J = jacobian(e);
            for (auto q : quad_points(test.Px, test.Py)) {
                double W = weight(q);
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


    void apply_bc(vector_view& Rvx, vector_view& Rvy) {
        // zero_bc(Rvx, test.U1x, test.U1y);
        // zero_bc(Rvy, test.U2x, test.U2y);

        // vx = 0 on left/right edge
        for (auto i = 0; i < test.U1y.dofs(); ++ i) {
            Rvx(0, i) = 0;
            Rvx(test.U1x.dofs() - 1, i) = 0;
        }

        // vy = 0 on top/bottom edge
        for (auto i = 0; i < test.U2x.dofs(); ++ i) {
            Rvy(i, 0) = 0;
            Rvy(i, test.U2y.dofs() - 1) = 0;
        }
    }

    void step(int /*iter*/, double /*t*/) override {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dP = trial.Px.dofs() * trial.Py.dofs();
        auto dim_trial = dU1 + dU2 + dP;

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto dim_test = DU1 + DU2;

        std::vector<double> rhs(dim_test + dim_trial);

        vector_view Rvx{rhs.data(), {test.U1x.dofs(), test.U1y.dofs()}};
        vector_view Rvy{Rvx.data() + DU1, {test.U2x.dofs(), test.U2y.dofs()}};

        vector_view vx{rhs.data() + dim_test, {trial.U1x.dofs(), trial.U1y.dofs()}};
        vector_view vy{vx.data() + dU1, {trial.U2x.dofs(), trial.U2y.dofs()}};
        vector_view p{vy.data() + dU2, {trial.Px.dofs(), trial.Py.dofs()}};


        mumps::problem problem(rhs.data(), rhs.size());

        std::cout << "Assembling matrix" << std::endl;
        assemble_matrix(problem);

        std::cout << "Computing RHS" << std::endl;
        compute_rhs(Rvx, Rvy);
        apply_bc(Rvx, Rvy);
        int i = linear_index({0, 0}, trial.Px, trial.Py);
        p(i, i) = 0; // fix pressure at a point

        std::cout << "Solving" << std::endl;
        solver.solve(problem);

        std::cout << "Error:" << std::endl;
        print_error(vx, vy, p);

        std::cout << "Pressure integral: " << total_pressure(p, trial.Px, trial.Py) << std::endl;

        std::cout << "Outputting" << std::endl;
        outputP.to_file(p, "pressure.data");
        outputU1.to_file(vx, "vx.data");
        outputU2.to_file(vy, "vy.data");
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
                double w = weight(q, space.Px, space.Py);
                auto x = point(e, q, space.Px, space.Py);
                value_type div = divergence(u, v, e, q, space);

                auto d = div - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol>
    double total_pressure(const Sol& p, const dimension& Px, const dimension& Py) const {
        double pressure = 0;
        for (auto e : elements(Px, Py)) {
            double J = jacobian(e, Px, Py);
            for (auto q : quad_points(Px, Py)) {
                double w = weight(q, Px, Py);
                value_type v = eval(p, e, q, Px, Py);
                pressure += v.val * w * J;
            }
        }
        return pressure;
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

};

}

#endif // STOKES_STOKES_CONSTRAINED_HPP
