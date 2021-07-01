// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_ERIKKSON_ERIKKSON_CG_WEAK_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_CG_WEAK_HPP_

#include <galois/Timer.h>

#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/output_manager.hpp"
#include "ads/solver/mumps.hpp"
#include "erikkson_base.hpp"
#include "solution.hpp"


namespace ads {

class erikkson_CG_weak : public erikkson_base {
private:
    using Base = erikkson_base;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;


    // New stuff
    lin::band_matrix MVx, MVy;
    lin::band_matrix MUx, MUy;

    lin::band_matrix KVx, KVy;
    lin::band_matrix KUx, KUy;


    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;

    lin::dense_matrix MUUx, MUUy;
    lin::dense_matrix KUUx, KUUy;
    lin::dense_matrix AUUx, AUUy;

    struct residuum {
        vector_type data;
        const dimension* Vx;
        const dimension* Vy;
    };

    vector_type u, u_prev;
    residuum r;
    vector_type u_buffer;
    std::vector<double> full_rhs;

    int save_every = 1;

    double h;
    double gamma;

    double peclet = 1e30;
    double epsilon = 1 / peclet;

    point_type c_diff{{ epsilon, epsilon }};

    double b = 1;

    // double angle = 0;
    // double angle = M_PI / 6;

    // double len = 1;

    // point_type beta{{ len * cos(angle), len * sin(angle) }};
    // point_type beta{{ 1, 1 }};
    // point_type beta{{ 1, 0 }};

    mumps::solver solver;

    output_manager<2> output;

    galois::StatTimer integration_timer{"integration"};
    galois::StatTimer solver_timer{"solver"};

public:
    erikkson_CG_weak(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
    : Base{std::move(test_x), std::move(test_y), steps}
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , Vx{ x }
    , Vy{ y }
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , MUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , KUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , MUVx{ Vx.dofs(), Ux.dofs() }
    , MUVy{ Vy.dofs(), Uy.dofs() }
    , KUVx{ Vx.dofs(), Ux.dofs() }
    , KUVy{ Vy.dofs(), Uy.dofs() }
    , AUVx{ Vx.dofs(), Ux.dofs() }
    , AUVy{ Vy.dofs(), Uy.dofs() }
    , MUUx{ Ux.dofs(), Ux.dofs() }
    , MUUy{ Uy.dofs(), Uy.dofs() }
    , KUUx{ Ux.dofs(), Ux.dofs() }
    , KUUy{ Uy.dofs(), Uy.dofs() }
    , AUUx{ Ux.dofs(), Ux.dofs() }
    , AUUy{ Uy.dofs(), Uy.dofs() }
    , u{{ Ux.dofs(), Uy.dofs() }}
    , u_prev{{ Ux.dofs(), Uy.dofs() }}
    , r{ vector_type{{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , h{ element_diam(Ux, Uy) }
    , output{ Ux.B, Uy.B, 500 }
    {
        int p = Vx.basis.degree;
        gamma = 3 * epsilon * p * p / h;
    }

private:

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    struct matrix_set {
        using band_matrix_ref = lin::band_matrix&;
        using dense_matrix_ref = lin::dense_matrix&;

        band_matrix_ref MVx, MVy, KVx, KVy;
        dense_matrix_ref MUVx, MUVy, KUVx, KUVy, AUVx, AUVy;
    };

    bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
        return dof[0] == 0;
        // return false;
    }

    double diffusion(point_type /*x*/) const {
        return epsilon;

        // bool left = x < 1 - 1/peclet, right = !left;
        // bool bottom = y < 1 - 1/peclet, top = !bottom;

        // if ((bottom && left) || (top && right)) {
            // return peclet;
        // } else {
            // return 1 / peclet;
        // }
    }

    point_type beta(point_type x) const {
        // return {1, 0};
        double r = b / std::sqrt(x[0] * x[0] + x[1] * x[1]);
        return { - r * x[1], r * x[0] };
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

        if (ix < 0 || ix > x.p || iy < 0 || iy > y.p) {
            std::cout << "dof = (" << dof[0] << ", " << dof[1] << ") "
                      << "e = (" << e[0] << ", " << e[1] << ") "
                      << "(x, y) = (" << p[0] << ", " << p[1] << ") "
                      << "span = (" << spanx << ", " << spany << ") "
                      << "offset = (" << offsetx << ", " << offsety << ") "
                      << "i = (" << ix << ", " << iy << ")" << std::endl;
        }

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

    template <typename Fun, typename Form>
    double integrate_boundary(boundary side, index_type i, const dimension& Ux, const dimension& Uy,
                              Fun&& g, Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Uy.elements - 1;
            if (! supported_in_1d(i[1], ey, Uy)) return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++ q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type ww = eval_basis_at(x, {e, ey}, i, Ux, Uy);
                    auto gg = g(x);
                    double fuw = form(ww, gg);
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
                    value_type ww = eval_basis_at(x, {ex, e}, i, Ux, Uy);
                    auto gg = g(x);
                    double fuw = form(ww, gg);
                    val += fuw * w * J;
                }
            }
        }
        return val;
    }

    double dot(point_type a, point_type b) const {
        return a[0] * b[0] + a[1] * b[1];
    }

    double dot(value_type a, point_type b) const {
        return a.dx * b[0] + a.dy * b[1];
    }

    double dot(point_type a, value_type b) const {
        return a[0] * b.dx + a[1] * b.dy;
    }

    point_type normal(boundary side) const {
        switch (side) {
        case boundary::left:   return {-1,  0};
        case boundary::right:  return { 1,  0};
        case boundary::bottom: return { 0, -1};
        case boundary::top:    return { 0,  1};
        default:               return {0, 0};
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

    void assemble_problem(mumps::problem& problem, const dimension& Vx, const dimension& Vy, const matrix_set& /*M*/) {
        auto N = Vx.dofs() * Vy.dofs();
        auto hh = h * h;

        // Gram matrix - upper left
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : overlapping_dofs(i, Vx, Vy)) {
                if (is_fixed(i, Vx, Vy) || is_fixed(j, Vx, Vy)) continue;

                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double val = kron(MVx, MVy, i, j) + hh * (kron(KVx, MVy, i, j) + kron(MVx, KVy, i, j));
                // val += hh * hh * kron(M.KVx, M.KVy, i, j);
                problem.add(ii, jj, val);
            }
        }

        // B, B^T
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : dofs(Ux, Uy)) {
                if (is_fixed(i, Vx, Vy) || is_fixed(j, Ux, Uy)) continue;

                double val = 0;
                for (auto e : elements_supporting_dof(i, Vx, Vy)) {
                    if (! supported_in(j, e, Ux, Uy)) continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(Vx, Vy)) {
                        double w = weigth(q);
                        auto x = point(e, q);
                        value_type ww = eval_basis(e, q, i, Vx, Vy);
                        value_type uu = eval_basis(e, q, j, Ux, Uy);

                        auto diff = diffusion(x);
                        double bwu = diff * grad_dot(uu, ww) + dot(beta(x), uu) * ww.val;
                        val += bwu * w * J;
                    }
                }

                auto int_bd = [&](auto side, auto form) {
                    if (touches(i, side, Vx, Vy) && touches(j, side, Ux, Uy)) {
                        val += integrate_boundary(side, i, j, Vx, Vy, Ux, Uy, [&](auto w, auto u, auto x) {
                            return form(w, u, x, this->normal(side));
                        });
                    }
                };

                auto boundary_term = [&](auto form) {
                    int_bd(boundary::left, form);
                    int_bd(boundary::right, form);
                    int_bd(boundary::top, form);
                    int_bd(boundary::bottom, form);
                };

                // <eps \/u*n, v>
                boundary_term([&](auto w, auto u, auto  , auto n) { return - epsilon * w.val * dot(u, n); });
                // <u, eps \/v*n>
                boundary_term([&](auto w, auto u, auto  , auto n) { return - epsilon * u.val * dot(w, n); });
                // <u, v beta*n>
                boundary_term([&](auto w, auto u, point_type x, auto n) { return - u.val * w.val * dot(beta(x), n); });
                // <u, gamma v>
                boundary_term([&](auto w, auto u, auto  , auto  ) { return - u.val * w.val * gamma; });

                if (val != 0) {
                    int ii = linear_index(i, Vx, Vy) + 1;
                    int jj = linear_index(j, Ux, Uy) + 1;

                    problem.add(ii, N + jj, -val);
                    problem.add(N + jj, ii, val);
                }
            }
        }

        // Dirichlet BC - upper left
        for_boundary_dofs(Vx, Vy, [&](index_type dof) {
            if (is_fixed(dof, Vx, Vy)) {
                int i = linear_index(dof, Vx, Vy) + 1;
                problem.add(i, i, 1);
            }
        });

        // Dirichlet BC - lower right
        for_boundary_dofs(Ux, Uy, [&](index_type dof) {
            if (is_fixed(dof, Ux, Uy)) {
                int i = linear_index(dof, Ux, Uy) + 1;
                problem.add(N + i, N + i, 1);
            }
        });
    }


    matrix_set matrices(bool x_refined, bool y_refined) {
        if (x_refined && y_refined) {
            return { MVx, MVy, KVx, KVy, MUVx, MUVy, KUVx, KUVy, AUVx, AUVy };
        } else if (x_refined) {
            return { MVx, MUy, KVx, KUy, MUVx, MUUy, KUVx, KUUy, AUVx, AUUy };
        } else if (y_refined) {
            return { MUx, MVy, KUx, KVy, MUUx, MUVy, KUUx, KUVy, AUUx, AUVy };
        } else {
            return { MUx, MUy, KUx, KUy, MUUx, MUUy, KUUx, KUUy, AUUx, AUUy };
        }
    }

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);

        gram_matrix_1d(MUx, Ux.basis);
        gram_matrix_1d(MUy, Uy.basis);

        gram_matrix_1d(MUUx, Ux.basis, Ux.basis);
        gram_matrix_1d(MUUy, Uy.basis, Uy.basis);

        gram_matrix_1d(MUVx, Ux.basis, Vx.basis);
        gram_matrix_1d(MUVy, Uy.basis, Vy.basis);


        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);

        stiffness_matrix_1d(KUx, Ux.basis);
        stiffness_matrix_1d(KUy, Uy.basis);

        stiffness_matrix_1d(KUVx, Ux.basis, Vx.basis);
        stiffness_matrix_1d(KUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KUUx, Ux.basis, Ux.basis);
        stiffness_matrix_1d(KUUy, Uy.basis, Uy.basis);


        advection_matrix_1d(AUVx, Ux.basis, Vx.basis);
        advection_matrix_1d(AUVy, Uy.basis, Vy.basis);

        advection_matrix_1d(AUUx, Ux.basis, Ux.basis);
        advection_matrix_1d(AUUy, Uy.basis, Uy.basis);
    }

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();
        Vx.factorize_matrix();
        Vy.factorize_matrix();

        // auto init = [this](double x, double y) { return init_state(x, y); };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        zero(r.data);
        zero(u);

        // stationary_bc(u, Ux, Uy);
        // skew_bc(u, Ux, Uy);
        // zero_bc(u, Ux, Uy);

        dirichlet_bc(u, boundary::left, Ux, Uy, [&](double t) {
            double tt = std::abs(t);
            return 0.5 * (std::tanh(b / epsilon * (tt < 0.5 ? tt - 0.35 : 0.65 - tt)) + 1);
            // return std::sin(M_PI * t);
        });


        output.to_file(u, "out_0.data");
    }

    void add_solution(const vector_view& u_rhs, const vector_view& r_rhs, const dimension& Vx, const dimension& Vy) {
        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) += u_rhs(i[0], i[1]);
        }
        for (auto i : dofs(Vx, Vy)) {
            r.data(i[0], i[1]) += r_rhs(i[0], i[1]);
        }
    }

    double norm(const vector_view& u) const {
        double norm = 0;
        for (auto i : dofs(Ux, Uy)) {
            auto a = u(i[0], i[1]);
            norm += a * a;
        }
        norm /= (Ux.dofs() * Uy.dofs());
        return std::sqrt(norm);
    }

    double substep(bool x_refined, bool y_refined, double /*t*/) {
        dimension& Vx = x_refined ? this->Vx : Ux;
        dimension& Vy = y_refined ? this->Vy : Uy;

        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        integration_timer.start();
        // compute_rhs_nonstationary(Vx, Vy, r_rhs, u_rhs, t);
        compute_rhs(Vx, Vy, r_rhs, u_rhs);

        // BC
        dirichlet_bc(u_rhs, boundary::left, Ux, Uy, 0);
        // dirichlet_bc(u_rhs, boundary::bottom, Ux, Uy, 0);
        // dirichlet_bc(u, boundary::left, Ux, Uy, [](double t) { return std::sin(M_PI * t); });


        dirichlet_bc(r_rhs, boundary::left, Vx, Vy, 0);
        // dirichlet_bc(r_rhs, boundary::bottom, Vx, Vy, 0);

        // zero_bc(r_rhs, Vx, Vy);
        // zero_bc(u_rhs, Ux, Uy);

        int size = Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs();
        mumps::problem problem(full_rhs.data(), size);
        assemble_problem(problem, Vx, Vy, matrices(x_refined, y_refined));
        integration_timer.stop();

        solver_timer.start();
        solver.solve(problem, "matrix");
        solver_timer.stop();

        add_solution(u_rhs, r_rhs, Vx, Vy);
        return norm(u_rhs);
    }

    void step(int /*iter*/, double t) override {
        // bool xrhs[] = { false, true };

        // bool x_rhs = xrhs[iter % sizeof(xrhs)];
        // bool y_rhs = !x_rhs;
        // std::cout << iter << ": x " << (x_rhs ? "rhs" : "lhs") << ", y " << (y_rhs ? "rhs" : "lhs") << std::endl;
        // substep(x_rhs, y_rhs, true, true);
        // substep(x_rhs, y_rhs, !x_rhs, !y_rhs);

        // using std::swap;
        // swap(u, u_prev);
        // zero(u);

        // std::cout << "Step " << (iter + 1) << std::endl;
        constexpr auto max_iters = 1;
        for (int i = 1; ; ++ i) {
            auto norm = substep(true, true, t);
            // std::cout << "  substep " << i << ": |eta| = " << norm << std::endl;
            if (norm < 1e-7 || i >= max_iters) {
                break;
            }
        }
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % save_every == 0) {
            // std::cout << "Step " << (iter + 1) << " : " << errorL2(t) << " " << errorH1(t) << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data", u, Ux, Uy);
        // plot_horizontal("horizontal.data", 0, u, Ux, Uy);
        // plot_vertical("vertical.data", 0.2, u, Ux, Uy);

        double T = steps.dt * steps.step_count;
        std::cout << "L2 abs: " << errorL2_abs(T) << std::endl;
        std::cout << "H1 abs: " << errorH1_abs(T) << std::endl;
        std::cout << "L2 rel: " << errorL2(T) << std::endl;
        std::cout << "H1 rel: " << errorH1(T) << std::endl;
        std::cout << "integration: " << static_cast<double>(integration_timer.get()) << std::endl;
        std::cout << "solver:      " << static_cast<double>(solver_timer.get()) << std::endl;

        print_solution("solution.data", u, Ux, Uy);
    }

    void compute_rhs(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs) {
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weigth(q);
                double WJ = W * J;
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);
                value_type rr = eval(r.data, e, q, *r.Vx, *r.Vy);
                auto diff = diffusion(x);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    // double F = 1;
                    double F = 0;
                    // double F = erikkson2_forcing(x[0], x[1], epsilon);
                    double lv = F * v.val;

                    double val = -lv;

                    // Bu
                    val += diff * grad_dot(uu, v) + dot(beta(x), uu) * v.val;
                    // -Aw
                    val -= (rr.val * v.val + h * h * grad_dot(rr, v));

                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);
                    double val = 0;

                    // -B'w
                    val -= diff * grad_dot(w, rr) + dot(beta(x), w) * rr.val;
                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
        // // Boundary terms
        // for (auto i : dofs(Vx, Vy)) {
        //     double val = 0;
        //     auto g = [](auto x) { return x[0] == 0 ? std::sin(M_PI * x[1]) : 0; };

        //     auto int_bd = [&](auto side, auto form) {
        //         if (touches(i, side, Vx, Vy)) {
        //             val += integrate_boundary(side, i, Vx, Vy, g, [&](auto w, auto u) {
        //                 return form(w, u, this->normal(side));
        //             });
        //         }
        //     };

        //     auto boundary_term = [&](auto form) {
        //         int_bd(boundary::left, form);
        //         int_bd(boundary::right, form);
        //         int_bd(boundary::top, form);
        //         int_bd(boundary::bottom, form);
        //     };
        //     // <g, eps \/v*n>
        //     boundary_term([&](auto w, auto g, auto n) { return - epsilon * g * dot(w, n); });
        //     // <g, v beta*n>
        //     boundary_term([&](auto w, auto g, auto n) { return - g * w.val * dot(beta, n); });
        //     // <u, gamma v>
        //     boundary_term([&](auto w, auto g, auto  ) { return - g * w.val * gamma; });

        //     r_rhs(i[0], i[1]) += val;
        // }
    }

    void compute_rhs_nonstationary(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs,
                                   double t) {
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weigth(q);
                double WJ = W * J;

                auto x = point(e, q);

                value_type uu = eval(u, e, q, Ux, Uy);
                value_type uu_prev = eval(u_prev, e, q, Ux, Uy);
                value_type rr = eval(r.data, e, q, *r.Vx, *r.Vy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    double F = erikkson_forcing(x[0], x[1], epsilon, t + steps.dt);
                    double lv = (uu_prev.val + steps.dt * F) * v.val;

                    double val = -lv;

                    // Bu
                    val += uu.val * v.val;
                    val += steps.dt * (c_diff[0] * uu.dx * v.dx + beta(x)[0] * uu.dx * v.val);
                    val += steps.dt * (c_diff[1] * uu.dy * v.dy + beta(x)[1] * uu.dy * v.val);
                    // -Aw
                    val -= (rr.val * v.val + h * h * (rr.dx * v.dx + rr.dy * v.dy));

                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);
                    double val = 0;

                    // -B'w
                    val -= w.val * rr.val;
                    val -= steps.dt * (c_diff[0] * w.dx * rr.dx + beta(x)[0] * w.dx * rr.val);
                    val -= steps.dt * (c_diff[1] * w.dy * rr.dy + beta(x)[1] * w.dy * rr.val);

                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
    }

    double errorL2_abs(double /*t*/) const {
        auto sol = exact(epsilon);
        return Base::errorL2(u, Ux, Uy, sol);
    }

    double errorH1_abs(double /*t*/) const {
        auto sol = exact(epsilon);
        return Base::errorH1(u, Ux, Uy, sol);
    }

    double errorL2(double /*t*/) const {
        auto sol = exact(epsilon);
        // auto sol = [&](point_type x) { return erikkson2_exact(x[0], x[1], epsilon); };
        // auto sol = [&](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };

        return Base::errorL2(u, Ux, Uy, sol) / normL2(Ux, Uy, sol) * 100;
    }

    double errorH1(double /*t*/) const {
        auto sol = exact(epsilon);
        // auto sol = [&](point_type x) { return erikkson2_exact(x[0], x[1], epsilon); };
        // auto sol = [&](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };

        return Base::errorH1(u, Ux, Uy, sol) / normH1(Ux, Uy, sol) * 100;
    }
};

}



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_CG_WEAK_HPP */
