#ifndef PROBLEMS_CG_ADVECTION_HPP
#define PROBLEMS_CG_ADVECTION_HPP

#include "problems/erikkson/erikkson_base.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"
#include "problems/erikkson/solution.hpp"


#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class advection : public erikkson_base {
private:
    using Base = erikkson_base;
    using Base::errorL2;
    using Base::errorH1;

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

    vector_type u, r;
    std::vector<double> full_rhs;

    double h;
    double gamma;

    double peclet;
    double epsilon = 1 / peclet;

    double tol_outer = 1e-7;
    double tol_inner = 1e-7;
    int max_iters = 50;

    mumps::solver solver;

    output_manager<2> output;

    Galois::StatTimer integration_timer{"integration"};
    Galois::StatTimer solver_timer{"solver"};

public:
    advection(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, double peclet)
    : Base{std::move(test_x), std::move(test_y), timesteps_config(1, 0.0)}
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
    , r{{ Vx.dofs(), Vy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , h{ element_diam(Ux, Uy) }
    , peclet{ peclet }
    , output{ Ux.B, Uy.B, 500 }
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

    template <typename Form>
    double integrate_boundary(boundary side, index_type i, const dimension& Ux, const dimension& Uy,
                              const vector_type& u, const dimension& Vx, const dimension& Vy, Form&& form) const {
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
        auto eta = h * h;

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
                        double w = weigth(q);
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

        zero(r);
        zero(u);

        // auto init = [this](double x, double y) { return 2*x*x + 2*y*y; };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // vector_type u_buffer{{ Ux.dofs(), Uy.dofs() }};
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        apply_bc(u);
    }

    void update_solution(const vector_view& u_rhs, const vector_view& r_rhs, const dimension& Vx, const dimension& Vy) {
        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) += u_rhs(i[0], i[1]);
        }
        for (auto i : dofs(Vx, Vy)) {
            r(i[0], i[1]) += r_rhs(i[0], i[1]);
        }
    }

    double norm(const vector_view& u) const {
        double norm = 0;
        for (auto i : dofs(Ux, Uy)) {
            auto a = u(i[0], i[1]);
            norm += a * a;
        }
        auto n = Ux.dofs() * Uy.dofs();
        return std::sqrt(norm) / n;
    }

    double substep() {
        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        integration_timer.start();
        compute_rhs(Vx, Vy, r_rhs, u_rhs);

        // BC
        for_sides(dirichlet, [&](auto side) {
            dirichlet_bc(u_rhs, side, Ux, Uy, 0);
            dirichlet_bc(r_rhs, side, Vx, Vy, 0);
        });

        mumps::problem problem(full_rhs.data(), full_rhs.size());
        assemble_problem(problem, Vx, Vy);
        integration_timer.stop();

        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        update_solution(u_rhs, r_rhs, Vx, Vy);
        return norm(u_rhs);
    }

    void step(int /*iter*/, double /*t*/) override {
        for (int i = 1; ; ++ i) {
            auto norm = substep();
            std::cout << "  substep " << i << ": |eta| = " << norm << std::endl;

            if (norm < tol_outer || i >= max_iters) {
                break;
            }
        }
    }

    void after() override {
        output.to_file(u, "result.data");
        print_solution("coefficients.data", u, Ux, Uy);

        plot_middle("xsection.data", u, Ux, Uy);

        std::cout << "L2 abs: " << errorL2_abs() << std::endl;
        std::cout << "H1 abs: " << errorH1_abs() << std::endl;
        std::cout << "L2 rel: " << errorL2() << std::endl;
        std::cout << "H1 rel: " << errorH1() << std::endl;

        std::cout << "integration: " << static_cast<double>(integration_timer.get()) << " ms" << std::endl;
        std::cout << "solver:      " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
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

    double errorL2_abs() const {
        auto sol = exact(epsilon);
        return errorL2(u, Ux, Uy, sol);
    }

    double errorH1_abs() const {
        auto sol = exact(epsilon);
        return errorH1(u, Ux, Uy, sol);
    }

    double errorL2() const {
        auto sol = exact(epsilon);
        return errorL2(u, Ux, Uy, sol) / normL2(Ux, Uy, sol) * 100;
    }

    double errorH1() const {
        auto sol = exact(epsilon);
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

    void apply_bc(vector_type& u) {
        for_sides(dirichlet, [&](auto side) {
            dirichlet_bc(u, side, Ux, Uy, [&](double t) {
                auto x = boundary_point(t, side);
                return g(x);
            });
        });
    }



    // --------------------
    // Problem definition
    // --------------------

    double diffusion(point_type /*x*/) const {
        return epsilon;
    }

    point_type beta(point_type /*x*/) const {
        return {1, 0};
    }

    double F(point_type /*x*/) const {
        return 0;
    }

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
            - epsilon * g(x) * dot(v, n)       // <g, eps \/v*n>
            - g(x) * v.val * dot(beta(x), n)   // <g, v beta*n>
            - g(x) * v.val * gamma;            // <u, gamma v> -- penalty
    }

    // Scalar product

    double A(value_type u, value_type v) const {
        return u.val * v.val + h * h * grad_dot(u, v);
    }

    // Boundary conditions

    static constexpr boundary dirichlet = boundary::none;

    double g(point_type x) const {
        return x[0] == 0 ? std::sin(M_PI * x[1]) : 0;
    }



};

}

#endif // PROBLEMS_CG_ADVECTION_HPP
