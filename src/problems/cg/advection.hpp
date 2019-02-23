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

    double peclet;
    double epsilon = 1 / peclet;

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
    , u_prev{{ Ux.dofs(), Uy.dofs() }}
    , r{ {{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
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

    struct matrix_set {
        using band_matrix_ref = lin::band_matrix&;
        using dense_matrix_ref = lin::dense_matrix&;

        band_matrix_ref MVx, MVy, KVx, KVy;
        dense_matrix_ref MUVx, MUVy, KUVx, KUVy, AUVx, AUVy;
    };


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

    template <typename Form>
    double integrate_boundary(boundary side, index_type i, const dimension& Ux, const dimension& Uy, Form&& form) const {
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
                    double fuw = form(ww, x);
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
                    double fuw = form(ww, x);
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
        }
        return {0, 0};
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

    void assemble_problem(mumps::problem& problem, const dimension& Vx, const dimension& Vy, const matrix_set& M) {
        auto N = Vx.dofs() * Vy.dofs();
        auto eta = h * h;

        // Gram matrix - upper left
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : overlapping_dofs(i, Vx, Vy)) {
                if (is_fixed(i, Vx, Vy) || is_fixed(j, Ux, Uy)) continue;

                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double val = kron(MVx, MVy, i, j) + eta * (kron(KVx, MVy, i, j) + kron(MVx, KVy, i, j));
                // val += eta * eta * kron(M.KVx, M.KVy, i, j);
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

                auto int_bd = [&](auto side, auto form) {
                    if (touches(i, side, Vx, Vy) && touches(j, side, Ux, Uy)) {
                        val += integrate_boundary(side, i, j, Vx, Vy, Ux, Uy, [&](auto w, auto u, auto x) {
                            return form(w, u, x, this->normal(side));
                        });
                    }
                };

                auto boundary_term = [&](auto form) {
                    // int_bd(boundary::left, form);
                    int_bd(boundary::right, form);
                    int_bd(boundary::top, form);
                    int_bd(boundary::bottom, form);
                };

                // penalty & symmetrization
                boundary_term([&](auto w, auto u, point_type x, auto n) {
                    return
                        - epsilon * w.val * dot(u, n)         // <eps \/u*n, v> -- consistency
                        - epsilon * u.val * dot(w, n)         // <u, eps \/v*n>
                        - u.val * w.val * dot(beta(x), n)     // <u, v beta*n>
                        - u.val * w.val * gamma;              // <u, gamma v>   -- penalty
                });

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

        apply_bc(u);

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
        auto n = Ux.dofs() * Uy.dofs();
        return std::sqrt(norm) / n;
    }

    double substep(bool x_refined, bool y_refined, double t) {
        dimension& Vx = x_refined ? this->Vx : Ux;
        dimension& Vy = y_refined ? this->Vy : Uy;

        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        integration_timer.start();
        compute_rhs(Vx, Vy, r_rhs, u_rhs);

        // BC
        dirichlet_bc(u_rhs, boundary::left, Ux, Uy, 0);
        // dirichlet_bc(u_rhs, boundary::bottom, Ux, Uy, 0);
        // dirichlet_bc(u, boundary::left, Ux, Uy, [](double t) { return std::sin(M_PI * t); });


        dirichlet_bc(r_rhs, boundary::left, Vx, Vy, 0);
        // dirichlet_bc(r_rhs, boundary::bottom, Vx, Vy, 0);

        // zero_bc(r_rhs, Vx, Vy);
        // zero_bc(u_rhs, Ux, Uy);

        mumps::problem problem(full_rhs.data(), full_rhs.size());
        assemble_problem(problem, Vx, Vy, matrices(x_refined, y_refined));
        integration_timer.stop();

        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        add_solution(u_rhs, r_rhs, Vx, Vy);
        return norm(u_rhs);
    }

    void step(int iter, double t) override {
        constexpr auto max_iters = 1;

        for (int i = 1; ; ++ i) {
            auto norm = substep(true, true, t);
            std::cout << "  substep " << i << ": |eta| = " << norm << std::endl;

            if (norm < 1e-7 || i >= max_iters) {
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

        std::cout << "integration: " << static_cast<double>(integration_timer.get()) << std::endl;
        std::cout << "solver:      " << static_cast<double>(solver_timer.get()) << std::endl;
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

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    double lv = F(x) * v.val;

                    // -F + Bu - Ar
                    double val = -lv + B(uu, v, x) - A(rr, v);
                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);

                    // -B'w
                    double val = -B(w, rr, x);
                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
        // Boundary terms
        for (auto i : dofs(Vx, Vy)) {
            double val = 0;

            auto int_bd = [&](auto side, auto form) {
                if (touches(i, side, Vx, Vy)) {
                    val += integrate_boundary(side, i, Vx, Vy, [&](auto w, auto x) {
                        return form(w, x, this->normal(side));
                    });
                }
            };

            auto boundary_term = [&](auto form) {
                // int_bd(boundary::left, form);
                int_bd(boundary::right, form);
                int_bd(boundary::top, form);
                int_bd(boundary::bottom, form);
            };

            boundary_term([&](auto w, point_type x, auto n) {
                return
                    - epsilon * g(x) * dot(w, n)       // <g, eps \/v*n>
                    - g(x) * w.val * dot(beta(x), n)   // <g, v beta*n>
                    - g(x) * w.val * gamma;            // <u, gamma v> -- penalty
            });

            r_rhs(i[0], i[1]) += val;
        }
    }

    double errorL2_abs() const {
        auto sol = exact(epsilon);
        return Base::errorL2(u, Ux, Uy, sol);
    }

    double errorH1_abs() const {
        auto sol = exact(epsilon);
        return Base::errorH1(u, Ux, Uy, sol);
    }

    double errorL2() const {
        auto sol = exact(epsilon);
        return errorL2(u, Ux, Uy, sol) / normL2(Ux, Uy, sol) * 100;
    }

    double errorH1() const {
        auto sol = exact(epsilon);
        return errorH1(u, Ux, Uy, sol) / normH1(Ux, Uy, sol) * 100;
    }

    // --------------------
    // Problem definition
    // --------------------

    double diffusion(point_type x) const {
        return epsilon;
    }

    point_type beta(point_type x) const {
        return {1, 0};
    }

    double F(point_type x) const {
        return 0;
    }

    double B(value_type u, value_type v, point_type x) const {
        auto diff = diffusion(x);
        return diff * grad_dot(u, v) + dot(beta(x), u) * v.val;
    }

    double A(value_type u, value_type v) const {
        return u.val * v.val + h * h * grad_dot(u, v);
    }

    // Boundary conditions

    double g(point_type x) const {
        return x[0] == 0 ? std::sin(M_PI * x[1]) : 0;
    }

    bool is_fixed(index_type dof, const dimension& x, const dimension& y) const {
        return dof[0] == 0;
    }

    void apply_bc(vector_type& u) {
        // stationary_bc(u, Ux, Uy);
        // skew_bc(u, Ux, Uy);
        // zero_bc(u, Ux, Uy);

        dirichlet_bc(u, boundary::left, Ux, Uy, [&](double t) {
            // double tt = std::abs(t);
            // return 0.5 * (std::tanh(b / epsilon * (tt < 0.5 ? tt - 0.35 : 0.65 - tt)) + 1);
            return std::sin(M_PI * t);
        });
    }


};

}

#endif // PROBLEMS_CG_ADVECTION_HPP
