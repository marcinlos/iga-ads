#ifndef PROBLEMS_ERIKKSON_ERIKKSON_CG_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_CG_HPP_


#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"


#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class erikkson_CG : public simulation_2d {
private:
    using Base = simulation_2d;

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

    double peclet = 1e6;
    double epsilon = 1 / peclet;

    point_type c_diff{{ epsilon, epsilon }};

    double angle = 0;
    // double angle = M_PI / 6;

    double len = 1;

    point_type beta{{ len * cos(angle), len * sin(angle) }};

    mumps::solver solver;

    output_manager<2> output;

public:
    erikkson_CG(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
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
    , r{ {{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , h{ element_diam(Ux, Uy) }
    , output{ Ux.B, Uy.B, 500 }
    { }

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

    void assemble_problem(mumps::problem& problem, const dimension& Vx, const dimension& Vy, const matrix_set& M) {
        auto N = Vx.dofs() * Vy.dofs();
        auto hh = h * h;

        // Gram matrix - upper left
        for (auto i : internal_dofs(Vx, Vy)) {
            for (auto j : overlapping_internal_dofs(i, Vx, Vy)) {
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
                double val = 0;
                // val += kron(M.MUVx, M.MUVy, i, j);
                val += /*steps.dt */ (c_diff[0] * kron(M.KUVx, M.MUVy, i, j) + beta[0] * kron(M.AUVx, M.MUVy, i, j));
                val += /*steps.dt */ (c_diff[1] * kron(M.MUVx, M.KUVy, i, j) + beta[1] * kron(M.MUVx, M.AUVy, i, j));

                if (val != 0) {
                    if (! is_boundary(i, Vx, Vy) && ! is_boundary(j, Ux, Uy)) {
                        int ii = linear_index(i, Vx, Vy) + 1;
                        int jj = linear_index(j, Ux, Uy) + 1;

                        problem.add(ii, N + jj, -val);
                        problem.add(N + jj, ii, val);
                    }
                }
            }
        }

        // Dirichlet BC - upper left
        for_boundary_dofs(Vx, Vy, [&](index_type dof) {
            int i = linear_index(dof, Vx, Vy) + 1;
            problem.add(i, i, 1);
        });

        // Dirichlet BC - lower right
        for_boundary_dofs(Ux, Uy, [&](index_type dof) {
            int i = linear_index(dof, Ux, Uy) + 1;
            problem.add(N + i, N + i, 1);
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

        // auto init = [this](double x, double y) { return init_state(x, y); };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        zero(r.data);
        zero(u);
        apply_bc(u);

        output.to_file(u, "out_0.data");
    }

    void apply_bc(vector_type& u_rhs) {
        dirichlet_bc(u_rhs, boundary::left, Ux, Uy, [](double t) { return std::sin(M_PI * t); });
        dirichlet_bc(u_rhs, boundary::right, Ux, Uy, 0);
        dirichlet_bc(u_rhs, boundary::bottom, Ux, Uy, 0);
        dirichlet_bc(u_rhs, boundary::top, Ux, Uy, 0);
    }

    void zero_bc(vector_view& r_rhs, vector_view& u_rhs) {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u_rhs(i[0], i[1]) = 0; });
        for_boundary_dofs(Vx, Vy, [&](index_type i) { r_rhs(i[0], i[1]) = 0; });
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

    double substep(bool x_refined, bool y_refined, double t) {
        dimension& Vx = x_refined ? this->Vx : Ux;
        dimension& Vy = y_refined ? this->Vy : Uy;

        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        // compute_rhs_nonstationary(Vx, Vy, r_rhs, u_rhs, t);
        compute_rhs(Vx, Vy, r_rhs, u_rhs);

        zero_bc(r_rhs, u_rhs);

        int size = Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs();
        mumps::problem problem(full_rhs.data(), size);
        assemble_problem(problem, Vx, Vy, matrices(x_refined, y_refined));
        solver.solve(problem);

        add_solution(u_rhs, r_rhs, Vx, Vy);
        return norm(u_rhs);
    }

    void step(int iter, double t) override {
        // bool xrhs[] = { false, true };

        // bool x_rhs = xrhs[iter % sizeof(xrhs)];
        // bool y_rhs = !x_rhs;
        // std::cout << iter << ": x " << (x_rhs ? "rhs" : "lhs") << ", y " << (y_rhs ? "rhs" : "lhs") << std::endl;
        // substep(x_rhs, y_rhs, true, true);
        // substep(x_rhs, y_rhs, !x_rhs, !y_rhs);

        // using std::swap;
        // swap(u, u_prev);
        // zero(u);

        std::cout << "Step " << (iter + 1) << std::endl;
        constexpr auto max_iters = 50;
        for (int i = 1; ; ++ i) {
            auto norm = substep(true, true, t);
            std::cout << "  substep " << i << ": |eta| = " << norm << std::endl;
            if (norm < 1e-7 || i >= max_iters) {
                break;
            }
        }
    }

    void after_step(int iter, double t) override {
        if ((iter + 1) % save_every == 0) {
            std::cout << "Step " << (iter + 1) << " : " << errorL2(t) << " " << errorH1(t) << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data");
        double T = steps.dt * steps.step_count;
        std::cout << "{ 'L2': '" << errorL2(T) << "', 'H1': '" << errorH1(T) << "'}" << std::endl;

        std::ofstream sol("solution.data");
        for (int i = 0; i < Ux.dofs(); ++ i) {
            for (int j = 0; j < Uy.dofs(); ++ j) {
                sol << i << " " << j << " " << u(i, j) << std::endl;
            }
        }
    }

    void plot_middle(const char* filename) {
        std::ofstream out{filename};
        bspline::eval_ctx ctx_x{ Ux.B.degree }, ctx_y{ Uy.B.degree };

        auto print = [&](double xx) {
            auto val = bspline::eval(xx, 0.5, u, Ux.B, Uy.B, ctx_x, ctx_y);
            out << std::setprecision(16) << xx << " " << val << std::endl;
        };

        print(0);
        auto N = Ux.basis.quad_order;
        for (auto e : Ux.element_indices()) {
            std::vector<double> qs(Ux.basis.x[e], Ux.basis.x[e] + N);
            std::sort(begin(qs), end(qs));
            for (auto xx : qs) {
                print(xx);
            }
        }
        print(1);
    }

    void compute_rhs(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs) {
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weigth(q);
                double WJ = W * J;
                // auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);
                value_type rr = eval(r.data, e, q, *r.Vx, *r.Vy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    double lv = 0;//* v.val;

                    double val = -lv;

                    // Bu
                    val += c_diff[0] * uu.dx * v.dx + beta[0] * uu.dx * v.val;
                    val += c_diff[1] * uu.dy * v.dy + beta[1] * uu.dy * v.val;
                    // -Aw
                    val -= (rr.val * v.val + h * h * (rr.dx * v.dx + rr.dy * v.dy));

                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);
                    double val = 0;

                    // -B'w
                    val -= (c_diff[0] * w.dx * rr.dx + beta[0] * w.dx * rr.val);
                    val -= (c_diff[1] * w.dy * rr.dy + beta[1] * w.dy * rr.val);

                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
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

                    double F = forcing(x[0], x[1], t + steps.dt);
                    double lv = (uu_prev.val + steps.dt * F) * v.val;

                    double val = -lv;

                    // Bu
                    val += uu.val * v.val;
                    val += steps.dt * (c_diff[0] * uu.dx * v.dx + beta[0] * uu.dx * v.val);
                    val += steps.dt * (c_diff[1] * uu.dy * v.dy + beta[1] * uu.dy * v.val);
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
                    val -= steps.dt * (c_diff[0] * w.dx * rr.dx + beta[0] * w.dx * rr.val);
                    val -= steps.dt * (c_diff[1] * w.dy * rr.dy + beta[1] * w.dy * rr.val);

                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
    }

    double forcing(double x, double y, double t) {
        constexpr double pi = M_PI;
        auto s = [](double a) { return std::sin(pi * a); };
        auto c = [](double a) { return std::cos(pi * a); };
        return pi*s(x)*s(y)*c(t) + 2*pi*pi*epsilon*s(x)*s(y)*s(t) + pi*c(x)*s(y)*s(t);
    }

    value_type exact_nonstationary(double x, double y, double t) const {
        constexpr double pi = M_PI;
        auto s = [](double a) { return std::sin(pi * a); };
        auto c = [](double a) { return std::cos(pi * a); };

        return value_type{ s(t)*s(x)*s(y), pi*s(t)*c(x)*s(y), pi*s(t)*s(x)*c(y) };
    }

    value_type exact(double x, double y, double eps) const {
        using std::exp;
        auto lambda = M_PI * eps;
        auto del = std::sqrt(1 + 4 * lambda * lambda);
        auto r1 = (1 + del) / (2 * eps);
        auto r2 = (1 - del) / (2 * eps);

        auto norm = exp(-r1) - exp(-r2);
        auto alpha = (exp(r1 * (x - 1)) - exp(r2 * (x - 1))) / norm;
        double val = alpha * std::sin(M_PI * y);

        return {
            val,
            (r1 * exp(r1 * (x - 1)) - r2 * exp(r2 * (x - 1))) / norm * std::sin(M_PI * y),
            alpha * M_PI * std::cos(M_PI * y)
        };
    }

    double errorL2(double t) const {
        double error = 0;

        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e);
            for (auto q : quad_points(Ux, Ux)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - exact(x[0], x[1], epsilon);
                // auto d = uu - exact_nonstationary(x[0], x[1], t);
                error += d.val * d.val * w * J;
            }
        }
        return std::sqrt(error);
    }

    double errorH1(double t) const {
        double error = 0;
        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e);
            for (auto q : quad_points(Ux, Ux)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - exact(x[0], x[1], epsilon);
                // auto d = uu - exact_nonstationary(x[0], x[1], t);
                error += (d.val * d.val + d.dx * d.dx + d.dy * d.dy) * w * J;
            }
        }
        return std::sqrt(error);
    }
};

}



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_CG_HPP */
