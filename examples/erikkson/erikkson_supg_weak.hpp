// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ERIKKSON_ERIKKSON_SUPG_WEAK_HPP
#define ERIKKSON_ERIKKSON_SUPG_WEAK_HPP

#include <galois/Timer.h>

#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/output_manager.hpp"
#include "ads/solver/mumps.hpp"
#include "erikkson_base.hpp"

namespace ads {

class erikkson_supg_weak : public erikkson_base {
private:
    using Base = erikkson_base;

    galois_executor executor{8};

    // New stuff
    lin::band_matrix Mx, My;
    lin::band_matrix Kx, Ky;
    lin::band_matrix Ax, Ay;

    vector_type u, rhs;

    int save_every = 1;

    double peclet = 1e6;
    double epsilon = 1 / peclet;

    // double C1 = 4, C2 = 2;

    // point_type c_diff{{ epsilon, epsilon }};

    // double angle = 0;
    // double angle = M_PI / 6;

    // double len = 1;

    double hmax;

    // point_type beta{{ len * cos(angle), len * sin(angle) }};
    point_type beta{{1, 0}};

    mumps::solver solver;

    output_manager<2> output;

    galois::StatTimer integration_timer{"integration"};
    galois::StatTimer solver_timer{"solver"};

public:
    erikkson_supg_weak(const dimension& trial_x, const dimension& trial_y,
                       const timesteps_config& steps)
    : Base{trial_x, trial_y, steps}
    , Mx{x.p, x.p, x.dofs(), x.dofs(), 0}
    , My{y.p, y.p, y.dofs(), y.dofs(), 0}
    , Kx{x.p, x.p, x.dofs(), x.dofs(), 0}
    , Ky{y.p, y.p, y.dofs(), y.dofs(), 0}
    , Ax{x.p, x.p, x.dofs(), x.dofs(), 0}
    , Ay{y.p, y.p, y.dofs(), y.dofs(), 0}
    , u{{x.dofs(), y.dofs()}}
    , rhs{{x.dofs(), y.dofs()}}
    , output{x.B, y.B, 500} {
        hmax = element_diam(x, y);
    }

private:
    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    double elem_diam(index_type e) const {
        double hx = 2 * x.basis.J[e[0]];
        double hy = 2 * y.basis.J[e[1]];
        // double h = hx;
        return std::sqrt(hx * hy);
    }

    double dot(point_type a, point_type b) const { return a[0] * b[0] + a[1] * b[1]; }

    double dot(value_type a, point_type b) const { return a.dx * b[0] + a.dy * b[1]; }

    double dot(point_type a, value_type b) const { return a[0] * b.dx + a[1] * b.dy; }

    double diffusion(double /*x*/, double /*y*/) const {
        return 1 / peclet;
        // constexpr double eta = 1e6;
        // bool left = x < 0.5, right = !left;
        // bool bottom = y < 0.5, top = !bottom;

        // if ((bottom && left) || (top && right)) {
        //     return eta;
        // } else {
        //     return 1;
        // }
    }

    bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
        // return dof[0] == 0 || dof[0] == x.dofs() - 1 || dof[1] == 0 || dof[1] == y.dofs() - 1;
        return dof[0] == 0;
    }

    value_type eval_basis_at(point_type p, index_type e, index_type dof) const {
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

        if (ix < 0 || ix > x.p || iy < 0 || iy > y.p) {
            std::cout << "dof = (" << dof[0] << ", " << dof[1] << ") "
                      << "e = (" << e[0] << ", " << e[1] << ") "
                      << "(x, y) = (" << p[0] << ", " << p[1] << ") "
                      << "span = (" << spanx << ", " << spany << ") "
                      << "offset = (" << offsetx << ", " << offsety << ") "
                      << "i = (" << ix << ", " << iy << ")" << std::endl;
        }

        auto value = bvx[0][ix] * bvy[0][iy];
        auto dx = bvx[1][ix] * bvy[0][iy];
        auto dy = bvx[0][ix] * bvy[1][iy];

        return {value, dx, dy};
    }

    bool supported_in_1d(int dof, int e, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return e >= xrange.first && e <= xrange.second;
    }

    template <typename Form>
    double integrate_boundary(boundary side, index_type i, index_type j, const dimension& Ux,
                              const dimension& Uy, const dimension& Vx, const dimension& Vy,
                              Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Uy.elements - 1;
            if (!supported_in_1d(j[1], ey, Vy) || !supported_in_1d(i[1], ey, Uy))
                return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                if (!supported_in_1d(j[0], e, Vx))
                    continue;

                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type ww = eval_basis_at(x, {e, ey}, i);
                    value_type uu = eval_basis_at(x, {e, ey}, j);
                    double fuw = form(ww, uu);
                    val += fuw * w * J;
                }
            }
        } else {
            int ex = side == boundary::left ? 0 : Ux.elements - 1;
            if (!supported_in_1d(j[0], ex, Vx) || !supported_in_1d(i[0], ex, Ux))
                return 0;

            auto x0 = side == boundary::left ? Ux.a : Ux.b;

            for (auto e : Uy.basis.element_range(i[1])) {
                if (!supported_in_1d(j[1], e, Vy))
                    continue;

                double J = Uy.basis.J[e];

                for (int q = 0; q < Uy.basis.quad_order; ++q) {
                    double w = Uy.basis.w[q];
                    point_type x{x0, Uy.basis.x[e][q]};
                    value_type ww = eval_basis_at(x, {ex, e}, i);
                    value_type uu = eval_basis_at(x, {ex, e}, j);
                    double fuw = form(ww, uu);
                    val += fuw * w * J;
                }
            }
        }
        return val;
    }

    point_type normal(boundary side) const {
        switch (side) {
        case boundary::left: return {-1, 0};
        case boundary::right: return {1, 0};
        case boundary::bottom: return {0, -1};
        case boundary::top: return {0, 1};
        default: return {0, 0};
        }
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

    void assemble_problem(mumps::problem& problem, double /*dt*/) {
        for (auto a : dofs(x, y)) {
            for (auto b : overlapping_dofs(a, x, y)) {
                if (is_fixed(a, x, y))
                    continue;

                double val = 0;
                for (auto e : elements_supporting_dof(a, x, y)) {
                    if (!supported_in(b, e, x, y))
                        continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(x, y)) {
                        double w = weight(q, x, y);
                        auto pt = point(e, q, x, y);
                        value_type ww = eval_basis(e, q, a, x, y);
                        value_type uu = eval_basis(e, q, b, x, y);

                        auto diff = diffusion(pt[0], pt[1]);
                        double bwu = diff * grad_dot(uu, ww) + dot(beta, uu) * ww.val;

                        // double h = elem_diam(e);
                        // double tau = 1 / (C1 * epsilon / (h * h) + C2 / h);
                        double hx = 2 * x.basis.J[e[0]];
                        double hy = 2 * y.basis.J[e[1]];
                        double hh = hx * hx + hy * hy;
                        // double tau = \tau ^(-2) = ( bx/hx+by/hy)^2+ 3 \mu * 1 / (hx^2+hy^2)
                        double tau = 1. / std::sqrt(1 / (hx * hx) + 3 * diff / hh);

                        double bDw = dot(beta, ww);

                        double lap = laplacian(e, q, b, x, y);
                        double res = -epsilon * lap + dot(beta, uu);
                        double v = bDw * res;
                        val += (bwu + tau * v) * w * J;
                    }
                }

                auto int_bd = [&](auto side, auto form) {
                    if (touches(a, side, x, y) && touches(b, side, x, y)) {
                        val += integrate_boundary(side, a, b, x, y, x, y, [&](auto w, auto u) {
                            return form(w, u, this->normal(side));
                        });
                    }
                };

                auto boundary_term = [&](auto form) {
                    // int_bd(boundary::left, form);
                    int_bd(boundary::right, form);
                    int_bd(boundary::top, form);
                    int_bd(boundary::bottom, form);
                };

                int p = x.basis.degree;
                double gamma = 3 * epsilon * p * p / hmax;

                // <eps \/u*n, v>
                boundary_term([&](auto w, auto u, auto n) { return -epsilon * w.val * dot(u, n); });
                // <u, eps \/v*n>
                boundary_term([&](auto w, auto u, auto n) { return -epsilon * u.val * dot(w, n); });
                // <u, v beta*n>
                boundary_term(
                    [&](auto w, auto u, auto n) { return -u.val * w.val * dot(beta, n); });
                // <u, gamma v>
                boundary_term([&](auto w, auto u, auto) { return -u.val * w.val * gamma; });

                if (val != 0) {
                    int i = linear_index(a, x, y) + 1;
                    int j = linear_index(b, x, y) + 1;
                    problem.add(i, j, val);
                }
            }
        }

        // 1's for Dirichlet BC
        for_boundary_dofs(x, y, [&](index_type dof) {
            if (is_fixed(dof, x, y)) {
                int i = linear_index(dof, x, y) + 1;
                problem.add(i, i, 1);
            }
        });
    }

    void prepare_matrices() {
        gram_matrix_1d(Mx, x.basis);
        gram_matrix_1d(My, y.basis);

        stiffness_matrix_1d(Kx, x.basis);
        stiffness_matrix_1d(Ky, y.basis);

        advection_matrix_1d(Ax, x.basis);
        advection_matrix_1d(Ay, y.basis);
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min((dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
        // return 0;
    };

    void before() override {
        prepare_matrices();
        x.factorize_matrix();
        y.factorize_matrix();

        skew_bc(u, x, y);
        output.to_file(u, "out_0.data");
    }

    void step(int /*iter*/, double /*t*/) override {
        zero(rhs);

        integration_timer.start();

        compute_rhs();

        // stationary_bc(rhs, x, y);
        dirichlet_bc(rhs, boundary::left, x, y, [](double t) { return std::sin(M_PI * t); });

        // skew_bc(rhs, x, y);
        // zero_bc(rhs, x, y);

        mumps::problem problem(rhs.data(), rhs.size());
        assemble_problem(problem, steps.dt);
        integration_timer.stop();

        solver_timer.start();
        solver.solve(problem, "matrix");
        solver_timer.stop();

        u = rhs;
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % save_every == 0) {
            // std::cout << "Step " << (iter + 1) << " : " << errorL2() << " " << errorH1() <<
            // std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data", u, x, y);
        // std::cout << "{ 'L2': '" << errorL2() << "', 'H1': '" << errorH1() << "'}" << std::endl;
        std::cout << "L2 abs: " << errorL2_abs() << std::endl;
        std::cout << "H1 abs: " << errorH1_abs() << std::endl;
        std::cout << "L2 rel: " << errorL2() << std::endl;
        std::cout << "H1 rel: " << errorH1() << std::endl;
        std::cout << "integration: " << static_cast<double>(integration_timer.get()) << std::endl;
        std::cout << "solver:      " << static_cast<double>(solver_timer.get()) << std::endl;

        print_solution("solution.data", u, x, y);
    }

    void compute_rhs() {
        zero(rhs);
        executor.for_each(elements(x, y), [&](index_type e) {
            auto U = vector_type{{x.basis.dofs_per_element(), y.basis.dofs_per_element()}};

            double J = jacobian(e, x, y);
            for (auto q : quad_points(x, y)) {
                double w = weight(q, x, y);
                // auto x = point(e, q);

                for (auto a : dofs_on_element(e, x, y)) {
                    auto aa = dof_global_to_local(e, a, x, y);
                    value_type v = eval_basis(e, q, a, x, y);

                    // double F = 1;
                    double F = 0;
                    double val = F * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() { update_global_rhs(rhs, U, e, x, y); });
        });
    }

    double errorL2_abs() const {
        auto sol = exact(epsilon);
        return Base::errorL2(u, x, y, sol);
    }

    double errorH1_abs() const {
        auto sol = exact(epsilon);
        return Base::errorH1(u, x, y, sol);
    }

    double errorL2() const {
        auto sol = exact(epsilon);
        return Base::errorL2(u, x, y, sol) / normL2(x, y, sol) * 100;
    }

    double errorH1() const {
        auto sol = exact(epsilon);
        return Base::errorH1(u, x, y, sol) / normH1(x, y, sol) * 100;
    }
};

}  // namespace ads

#endif  // ERIKKSON_ERIKKSON_SUPG_WEAK_HPP
