// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CG_ADVECTION_HPP
#define CG_ADVECTION_HPP

#include <galois/Timer.h>

#include "../erikkson/erikkson_base.hpp"
#include "../erikkson/solution.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/output_manager.hpp"
#include "ads/solver/mumps.hpp"


namespace ads {

struct advection_config {
    double tol_outer = 1e-7;
    double tol_inner = 1e-7;
    int max_outer_iters = 100;
    int max_inner_iters = 250;

    bool use_cg = true;
    bool weak_bc = false;

    bool print_inner = false;
    bool print_inner_count = true;
    bool print_outer = true;
    bool print_outer_count = true;
    bool print_inner_total = true;

    bool print_times = false;
    bool print_errors = true;
    bool plot = true;

    int threads = 8;
};


template <typename Problem>
class advection : public erikkson_base {
public:


private:
    using Base = erikkson_base;
    using Base::errorL2;
    using Base::errorH1;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;

    lin::band_matrix MVx, MVy;
    lin::band_matrix KVx, KVy;

    lin::band_matrix Ax, Ay;
    lin::solver_ctx Ax_ctx, Ay_ctx;

    vector_type u, r;
    vector_type buffer;

    std::vector<double> full_rhs;

    double h;
    double eta;
    double gamma;

    double peclet;
    double epsilon = 1 / peclet;

    advection_config cfg;
    Problem problem;

    mumps::solver solver;

    output_manager<2> output;

    galois::StatTimer integration_timer{"integration"};
    galois::StatTimer solver_timer{"solver"};
    galois::StatTimer total_timer{"total"};

    int total_CG_iters = 0;

    boundary dirichlet;

public:
    advection(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, double peclet, double eta, const advection_config& cfg, Problem problem)
    : Base{std::move(test_x), std::move(test_y), timesteps_config(1, 0.0)}
    , executor{ cfg.threads }
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , Vx{ x }
    , Vy{ y }
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , Ax{Vx.p, Vx.p, Vx.dofs()}
    , Ay{Vy.p, Vy.p, Vy.dofs()}
    , Ax_ctx{ Ax }
    , Ay_ctx{ Ay }
    , u{{ Ux.dofs(), Uy.dofs() }}
    , r{{ Vx.dofs(), Vy.dofs() }}
    , buffer{{ Vx.dofs(), Vy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , h{ element_diam(Ux, Uy) }
    , eta{ eta }
    , peclet{ peclet }
    , cfg{ cfg }
    , problem{ problem }
    , output{ Ux.B, Uy.B, 500 }
    , dirichlet{ cfg.weak_bc ? boundary::none : boundary::full }
    {
        int p = Vx.basis.degree;
        gamma = 3 * epsilon * p * p / h;
    }

private:

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    value_type eval_basis_at(point_type p, index_type dof, const dimension& x, const dimension& y) const {
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

    template <typename Sol>
    value_type eval_at(point_type p, const Sol& v, const dimension& x, const dimension& y) const {
        bspline::eval_ders_ctx cx{x.p, 1};
        bspline::eval_ders_ctx cy{y.p, 1};
        return bspline::eval_ders(p[0], p[1], v, x.B, y.B, cx, cy);
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
        auto nv = normal(side);

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
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, j, Vx, Vy);
                    double fuw = form(ww, uu, x, nv);
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
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, j, Vx, Vy);
                    double fuw = form(ww, uu, x, nv);
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

    template <typename Form, typename Sol>
    double integrate_boundary(boundary side, index_type i, const dimension& Ux, const dimension& Uy,
                              const Sol& u, const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;
        auto nv = normal(side);

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Vy.elements - 1;
            if (! supported_in_1d(i[1], ey, Uy)) return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++ q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type uu = eval_at(x, u, Vx, Vy);
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    double fuw = form(uu, ww, x, nv);
                    val += fuw * w * J;
                }
            }
        } else {
            int ex = side == boundary::left ? 0 : Vx.elements - 1;
            if (! supported_in_1d(i[0], ex, Ux)) return 0;

            auto x0 = side == boundary::left ? Ux.a : Ux.b;

            for (auto e : Uy.basis.element_range(i[1])) {
                double J = Uy.basis.J[e];

                for (int q = 0; q < Uy.basis.quad_order; ++ q) {
                    double w = Uy.basis.w[q];
                    point_type x{x0, Uy.basis.x[e][q]};
                    value_type uu = eval_at(x, u, Vx, Vy);
                    value_type ww = eval_basis_at(x, i, Ux, Uy);
                    double fuw = form(uu, ww, x, nv);
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

    void assemble_problem(mumps::problem& problem, const dimension& Vx, const dimension& Vy) {
        auto N = Vx.dofs() * Vy.dofs();
        // auto eta = h * h;

        // Gram matrix - upper left
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : overlapping_dofs(i, Vx, Vy)) {
                if (is_fixed(i, Vx, Vy) || is_fixed(j, Vx, Vy)) continue;

                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double val = kron(MVx, MVy, i, j) + eta * (kron(KVx, MVy, i, j) + kron(MVx, KVy, i, j));
                val += eta * eta * kron(KVx, KVy, i, j);
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
                        double w = weight(q);
                        auto x = point(e, q);
                        value_type ww = eval_basis(e, q, i, Vx, Vy);
                        value_type uu = eval_basis(e, q, j, Ux, Uy);
                        val += B(uu, ww, x) * w * J;
                    }
                }

                auto form = [&](auto u, auto w, auto x, auto n) { return this->bdB(u, w, x, n); };
                for_sides(~dirichlet, [&](auto side) {
                    if (touches(i, side, Vx, Vy) && touches(j, side, Ux, Uy)) {
                        val += integrate_boundary(side, j, i, Ux, Uy, Vx, Vy, form);
                    }
                });

                if (val != 0) {
                    int ii = linear_index(i, Vx, Vy) + 1;
                    int jj = linear_index(j, Ux, Uy) + 1;

                    problem.add(ii, N + jj, val);
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

    void fix_dof(int k, const dimension& dim, lin::band_matrix& K) {
        int last = dim.dofs() - 1;
        for (int i = clamp(k - dim.p, 0, last); i <= clamp(k + dim.p, 0, last); ++ i) {
            K(k, i) = 0;
        }
        K(k, k) = 1;
    }

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);
        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);

        // double eta = h * h;

        for (int i = 0; i < Vx.dofs(); ++ i) {
            for (int j : overlapping_dofs(i, 0, Vx.dofs(), Vx)) {
                Ax(i, j) = MVx(i, j) + eta * KVx(i, j);
            }
        }
        for (int i = 0; i < Vy.dofs(); ++ i) {
            for (int j : overlapping_dofs(i, 0, Vy.dofs(), Vy)) {
                Ay(i, j) = MVy(i, j) + eta * KVy(i, j);
            }
        }

        if (contains(dirichlet, boundary::left))   fix_dof(0, Vx, Ax);
        if (contains(dirichlet, boundary::right))  fix_dof(Vx.dofs() - 1, Vx, Ax);
        if (contains(dirichlet, boundary::bottom)) fix_dof(0, Vy, Ay);
        if (contains(dirichlet, boundary::top))    fix_dof(Vy.dofs() - 1, Vy, Ay);
    }

    void before() override {
        prepare_matrices();

        Ux.factorize_matrix();
        Uy.factorize_matrix();
        Vx.factorize_matrix();
        Vy.factorize_matrix();

        lin::factorize(Ax, Ax_ctx);
        lin::factorize(Ay, Ay_ctx);

        zero(r);
        zero(u);

        // auto init = [this](double x, double y) { return 2*x*x + 2*y*y; };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // vector_type u_buffer{{ Ux.dofs(), Uy.dofs() }};
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        apply_bc(u, Ux, Uy);
    }

    template <typename U, typename V>
    void update(const U& c, const V& d) {
        for (auto i : dofs(Vx, Vy)) {
            r(i[0], i[1]) += d(i[0], i[1]);
        }
        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) += c(i[0], i[1]);
        }
    }

    template <typename U>
    void update_solution(const U& c) {
        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) += c(i[0], i[1]);
        }
    }

    template <typename U, typename V>
    void update_residual(const U& dd, const V& Bc) {
        for (auto i : dofs(Vx, Vy)) {
            r(i[0], i[1]) = dd(i[0], i[1]) - Bc(i[0], i[1]);
        }
    }

    template <typename U, typename V>
    double dot(const U& u, const V& v, const dimension& x, const dimension& y) const {
        double val = 0;
        for (auto i : dofs(x, y)) {
            val += u(i[0], i[1]) * v(i[0], i[1]);
        }
        return val;
    }

    template <typename V>
    double norm_sq(const V& u, const dimension& x, const dimension& y) const {
        return dot(u, u, x, y);
    }

    template <typename V>
    double norm(const V& u, const dimension& x, const dimension& y) const {
        return std::sqrt(norm_sq(u, x, y));
    }

    vector_view substep_mumps() {
        vector_view dr{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view du{full_rhs.data() + dr.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        integration_timer.start();
        compute_rhs(Vx, Vy, dr, du);

        // BC
        for_sides(dirichlet, [&](auto side) {
            dirichlet_bc(du, side, Ux, Uy, 0);
            dirichlet_bc(dr, side, Vx, Vy, 0);
        });

        mumps::problem problem(full_rhs.data(), full_rhs.size());
        assemble_problem(problem, Vx, Vy);
        integration_timer.stop();

        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        update_solution(du);
        return du;
    }

    vector_view substep_CG(const vector_type& dc) {
        vector_type u_prev = u;

        auto p = dc, q = dc;
        auto theta = vector_type{{ Vx.dofs(), Vy.dofs() }};
        auto delta = theta;

        auto Mp = vector_type{{ Ux.dofs(), Uy.dofs() }};

        int i = 0;
        for (; i < cfg.max_inner_iters; ++ i) {
            // theta = Bp
            apply_B(p, theta);
            zero_bc(theta, Vx, Vy);

            // delta = A~ \ theta
            delta = theta;
            solve_A(delta);

            // alpha = (p, q) / (p, Mp)
            double alpha = dot(p, q, Ux, Uy) / dot(theta, delta, Vx, Vy);

            // Mp = B' delta
            apply_Bt(delta, Mp);
            zero_bc(Mp, Ux, Uy);

            double qnorm2_prev = norm_sq(q, Ux, Uy);
            // u := u + alpha p
            // q := q - alpha Mp
            for (auto i : dofs(Ux, Uy)) {
                u(i[0], i[1]) += alpha * p(i[0], i[1]);
                q(i[0], i[1]) -= alpha * Mp(i[0], i[1]);
            }

            // beta = |q_prev|^2 / |q|^2
            double qnorm2 = norm_sq(q, Ux, Uy);
            double beta = qnorm2 / qnorm2_prev;

            // p := q + beta * p
            for (auto i : dofs(Ux, Uy)) {
                p(i[0], i[1]) = q(i[0], i[1]) + beta * p(i[0], i[1]);
            }

            // convergence
            auto dimU = Ux.dofs() * Uy.dofs();
            double residuum = std::sqrt(qnorm2) / dimU;

            if (cfg.print_inner)
                std::cout << "     inner " << (i + 1) << ": |q| = " << residuum << std::endl;

            if (residuum < cfg.tol_inner) {
                ++ i;
                break;
            }
        }
        total_CG_iters += i;
        if (cfg.print_inner_count)
            std::cout << "  CG iters: " << i << std::endl;


        vector_view du{full_rhs.data(), {Ux.dofs(), Uy.dofs()}};
        for (auto i : dofs(Ux, Uy)) {
            du(i[0], i[1]) = u(i[0], i[1]) - u_prev(i[0], i[1]);
        }
        return du;
    }

    void solve_A(vector_type& v) {
        ads_solve(v, buffer, dim_data{Ax, Ax_ctx}, dim_data{Ay, Ay_ctx});
    }

    void step(int /*iter*/, double /*t*/) override {
        auto dd = vector_type{{ Vx.dofs(), Vy.dofs() }};
        auto dc = vector_type{{ Ux.dofs(), Uy.dofs() }};
        auto Bc = vector_type{{ Vx.dofs(), Vy.dofs() }};

        total_timer.start();
        int i = 0;
        for (; i < cfg.max_outer_iters ; ++ i) {
            // dd = A~ \ (F + Kr - Bu)
            compute_dd(Vx, Vy, dd);
            zero_bc(dd, Vx, Vy);
            solve_A(dd);

            // dc = B' dd
            apply_Bt(dd, dc);
            zero_bc(dc, Ux, Uy);

            auto c = cfg.use_cg ? substep_CG(dc) : substep_mumps();

            // Bc = A~ \ B c
            apply_B(c, Bc);
            zero_bc(Bc, Vx, Vy);
            solve_A(Bc);

            // r + d = dd - A~ \ B c
            update_residual(dd, Bc);

            auto dimU = Ux.dofs() * Uy.dofs();
            auto cc = norm(c, Ux, Uy) / dimU;

            if (cfg.print_outer)
                std::cout << "  outer " << (i + 1) << ": |c| = " << cc << std::endl;

            if (cc < cfg.tol_outer) {
                ++ i;
                break;
            }
        }
        total_timer.stop();

        if (cfg.print_outer_count)
            std::cout << "outer iters: " << i << std::endl;

        if (cfg.print_inner_total)
            std::cout << "total CG iters: " << total_CG_iters << std::endl;
    }

    void after() override {
        if (cfg.plot) {
            output.to_file(u, "result.data");
            print_solution("coefficients.data", u, Ux, Uy);
            plot_middle("xsection.data", u, Ux, Uy);
        }

        if (cfg.print_errors) {
            std::cout << "L2 abs: " << errorL2_abs() << std::endl;
            std::cout << "H1 abs: " << errorH1_abs() << std::endl;
            std::cout << "L2 rel: " << errorL2() << "%" << std::endl;
            std::cout << "H1 rel: " << errorH1() << "%" << std::endl;
        }

        if (cfg.print_times) {
            std::cout << "integration: " << static_cast<double>(integration_timer.get()) << " ms" << std::endl;
            std::cout << "solver:      " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
            std::cout << "total:       " << static_cast<double>(total_timer.get()) << " ms" << std::endl;
        }
    }

    void compute_rhs(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs) {
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);
                value_type rr = eval(r, e, q, Vx, Vy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    double Lv = F(x) * v.val;

                    // F - Ar - Bu
                    double val = Lv - A(rr, v) - B(uu, v, x);
                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);

                    // -B'r
                    double val = -B(w, rr, x);
                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });

        // Boundary terms of -Bu
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto form = [&](auto v, auto u, auto x, auto n) { return this->bdB(u, v, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, u, Ux, Uy, form);
                }
            });
            r_rhs(i[0], i[1]) -= val;
        }

        // Boundary terms of -B'r
        for (auto i : dofs(Ux, Uy)) {
            double val = 0;

            auto form = [&](auto u, auto v, auto x, auto n) { return this->bdB(u, v, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Ux, Uy)) {
                    val += integrate_boundary(side, i, Ux, Uy, r, Vx, Vy, form);
                }
            });
            u_rhs(i[0], i[1]) -= val;
        }

        // Boundary terms of RHS
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto form = [&](auto w, auto x, auto n) { return this->bdL(w, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, form);
                }
            });
            r_rhs(i[0], i[1]) += val;
        }
    }

    double eval_basis_dxy(index_type e, index_type q, index_type dof, const dimension& x, const dimension& y) const  {
        auto loc = dof_global_to_local(e, dof);
        const auto& bx = x.basis;
        const auto& by = y.basis;
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];
        return dB1 * dB2;
    }

    double eval_fun_dxy(const vector_type& v, index_type e, index_type q, const dimension& x, const dimension& y) const {
        double u = 0;
        for (auto b : dofs_on_element(e)) {
            double c = v(b[0], b[1]);
            double B = eval_basis_dxy(e, q, b, x, y);
            u += c * B;
        }
        return u;
    }

    void compute_dd(const dimension& Vx, const dimension& Vy, vector_type& dd) {
        zero(dd);
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);
                double rxy = eval_fun_dxy(r, e, q, Vx, Vy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    double vxy = eval_basis_dxy(e, q, a, Vx, Vy);
                    // double eta = h * h;
                    double Lv = F(x) * v.val;

                    // F + Kr - Bu
                    double Kr = eta * eta * rxy * vxy;
                    double val = Lv + Kr - B(uu, v, x);
                    R(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(dd, R, e, Vx, Vy);
            });
        });

        // Boundary terms of -Bu
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto form = [&](auto v, auto u, auto x, auto n) { return this->bdB(u, v, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, u, Ux, Uy, form);
                }
            });
            dd(i[0], i[1]) -= val;
        }

        // Boundary terms of RHS
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto form = [&](auto w, auto x, auto n) { return this->bdL(w, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, form);
                }
            });
            dd(i[0], i[1]) += val;
        }
    }

    template <typename U, typename Res>
    void apply_B(const U& u, Res& result) {
        zero(result);
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto rhs = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    // Bu
                    double val = B(uu, v, x);
                    rhs(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(result, rhs, e, Vx, Vy);
            });
        });

        // Boundary terms of Bu
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto form = [&](auto v, auto u, auto x, auto n) { return this->bdB(u, v, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, u, Ux, Uy, form);
                }
            });
            result(i[0], i[1]) += val;
        }
    }

    template <typename U, typename Res>
    void apply_Bt(const U& r, Res& result) {
        zero(result);
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto rhs = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                auto x = point(e, q);
                value_type rr = eval(r, e, q, Vx, Vy);

                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);

                    // B' r
                    double val = B(w, rr, x);
                    rhs(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(result, rhs, e, Ux, Uy);
            });
        });

        // Boundary terms of B'r
        for (auto i : dofs(Ux, Uy)) {
            double val = 0;

            auto form = [&](auto u, auto v, auto x, auto n) { return this->bdB(u, v, x, n); };
            for_sides(~dirichlet, [&](auto side) {
                if (this->touches(i, side, Ux, Uy)) {
                    val += integrate_boundary(side, i, Ux, Uy, r, Vx, Vy, form);
                }
            });
            result(i[0], i[1]) += val;
        }
    }

    auto exact() const {
        return [&](point_type x) { return solution(x); };
    }

    double errorL2_abs() const {
        auto sol = exact();
        return errorL2(u, Ux, Uy, sol);
    }

    double errorH1_abs() const {
        auto sol = exact();
        return errorH1(u, Ux, Uy, sol);
    }

    double errorL2() const {
        auto sol = exact();
        return errorL2(u, Ux, Uy, sol) / normL2(Ux, Uy, sol) * 100;
    }

    double errorH1() const {
        auto sol = exact();
        return errorH1(u, Ux, Uy, sol) / normH1(Ux, Uy, sol) * 100;
    }

    bool is_fixed(index_type dof, const dimension& x, const dimension& y) const {
        if (dof[0] == 0            && contains(dirichlet, boundary::left))   return true;
        if (dof[0] == x.dofs() - 1 && contains(dirichlet, boundary::right))  return true;
        if (dof[1] == 0            && contains(dirichlet, boundary::bottom)) return true;
        if (dof[1] == y.dofs() - 1 && contains(dirichlet, boundary::top))    return true;

        return false;
    }

    template <typename Fun>
    void for_sides(boundary sides, Fun&& f) {
        auto g = [&](auto s) { if (contains(sides, s)) f(s); };
        g(boundary::left);
        g(boundary::right);
        g(boundary::bottom);
        g(boundary::top);
    }

    point_type boundary_point(double t, boundary side) {
        double x0 = Ux.a, x1 = Ux.b;
        double y0 = Uy.a, y1 = Uy.b;
        switch (side) {
        case boundary::left:   return {x0, t};
        case boundary::right:  return {x1, t};
        case boundary::bottom: return {t, y0};
        case boundary::top:    return {t, y1};
        default:               return {0, 0};
        }
    }

    void apply_bc(vector_type& u, dimension& x, dimension& y) {
        for_sides(dirichlet, [&](auto side) {
            this->dirichlet_bc(u, side, x, y, [&](double t) {
                auto x = boundary_point(t, side);
                return this->g(x);
            });
        });
    }

    void zero_bc(vector_type& u, dimension& x, dimension& y) {
        for_sides(dirichlet, [&](auto side) { this->dirichlet_bc(u, side, x, y, 0); });
    }

    // Forms
    double B(value_type u, value_type v, point_type x) const {
        auto diff = diffusion(x);
        return diff * grad_dot(u, v) + dot(beta(x), u) * v.val;
    }

    double bdB(value_type u, value_type v, point_type x, point_type n) const {
        return
            - epsilon * v.val * dot(u, n)         // <eps \/u*n, v> -- consistency
            - epsilon * u.val * dot(v, n)         // <u, eps \/v*n> -- symmetrization
            - u.val * v.val * dot(beta(x), n)     // <u, v beta*n>  -- symmetrization
            - u.val * v.val * gamma;              // <u, gamma v>   -- penalty
    }

    double bdL(value_type v, point_type x, point_type n) const {
        return
            - epsilon * g(x) * dot(v, n)          // <g, eps \/v*n>
            - g(x) * v.val * dot(beta(x), n)      // <g, v beta*n>
            - g(x) * v.val * gamma;               // <u, gamma v> -- penalty
    }

    // Scalar product
    double A(value_type u, value_type v) const {
        return u.val * v.val + h * h * grad_dot(u, v);
    }

    // --------------------
    // Problem definition
    // --------------------

    double diffusion(point_type x) const {
        return problem.diffusion(x);
    }

    point_type beta(point_type x) const {
        return problem.beta(x);
    }

    // forcing
    double F(point_type x) const {
        return problem.forcing(x);
    }

    // Boundary conditions
    double g(point_type x) const {
        return problem.boundary(x);
    }

    // Exact solution
    value_type solution(point_type x) const {
        return problem.solution(x);
    }

};

}


#endif // CG_ADVECTION_HPP
