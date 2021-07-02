// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ERIKKSON_POLLUTION_ROTATION_HPP
#define ERIKKSON_POLLUTION_ROTATION_HPP

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

class pollution_rotation : public erikkson_base {
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

    double peclet = 1e6;
    double epsilon = 1 / peclet;

    mumps::solver solver;

    output_manager<2> output;

    galois::StatTimer solver_timer{"solver"};

public:
    pollution_rotation(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
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

    double diffusion(double /*x*/, double /*y*/) const {
        return epsilon;
        // const double eta = epsilon;
        // return 1 / eta;

        // bool left = x < 0.5, right = !left;
        // bool bottom = y < 0.5, top = !bottom;

        // if ((bottom && left) || (top && right)) {
            // return 1 / eta;
        // } else {
            // return eta;
        // }
    }

    void assemble_problem2(mumps::problem& problem, const dimension& Vx, const dimension& Vy, const matrix_set& /*M*/) {
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
        for (auto i : internal_dofs(Vx, Vy)) {
            for (auto j : internal_dofs(Ux, Uy)) {

                double val = 0;
                for (auto e : elements_supporting_dof(i, Vx, Vy)) {
                    if (! supported_in(j, e, Ux, Uy)) continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(Vx, Vy)) {
                        double w = weight(q);
                        auto x = point(e, q);
                        value_type ww = eval_basis(e, q, i, Vx, Vy);
                        value_type uu = eval_basis(e, q, j, Ux, Uy);

                        auto diff = diffusion(x[0], x[1]);
                        auto bx = -x[1];
                        auto by =  x[0];

                        double bwu = uu.val * ww.val +
                            steps.dt * (diff * grad_dot(uu, ww) + bx * uu.dx * ww.val + by * uu.dy * ww.val);

                        val += bwu * w * J;
                    }
                }

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

        // double F = 0;//source(x[0], x[1]);

        // auto init = [this](double x, double y) { return source(x, y); };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());
        u(Ux.dofs() / 2 + Ux.dofs() / 4, Uy.dofs() / 2) = 1;

        zero(r.data);
        // zero(u);

        // stationary_bc(u, Ux, Uy);
        // skew_bc(u, Ux, Uy);
        zero_bc(u, Ux, Uy);

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

    double substep(bool x_refined, bool y_refined, double t) {
        dimension& Vx = x_refined ? this->Vx : Ux;
        dimension& Vy = y_refined ? this->Vy : Uy;

        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);

        compute_rhs_nonstationary(Vx, Vy, r_rhs, u_rhs, t);

        zero_bc(r_rhs, Vx, Vy);
        zero_bc(u_rhs, Ux, Uy);

        int size = Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs();
        mumps::problem problem(full_rhs.data(), size);
        assemble_problem2(problem, Vx, Vy, matrices(x_refined, y_refined));

        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        add_solution(u_rhs, r_rhs, Vx, Vy);
        return norm(u_rhs);
    }

    void step(int iter, double t) override {
        using std::swap;
        swap(u, u_prev);
        zero(u);

        std::cout << "Step " << (iter + 1) << std::endl;
        constexpr auto max_iters = 1;//30;
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
        plot_middle("final.data", u, Ux, Uy);
        double T = steps.dt * steps.step_count;
        std::cout << "{ 'L2': '" << errorL2(T)
                  << "', 'H1': '" << errorH1(T)
                  << "', 'solver': " << static_cast<double>(solver_timer.get()) / steps.step_count
                  << ", 'dofs': " << Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs()
                  << "}" << std::endl;
        print_solution("solution.data", u, Ux, Uy);
    }

    double source(double x, double y) const {
        double dx = x - 0.5;
        double dy = y;
        double d = std::sqrt(dx * dx + dy * dy);
        double t = d * 20;
        return t > 1 ? 0 : (d - 1) * (d - 1) * (d + 1) * (d + 1);
    }

    void compute_rhs_nonstationary(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs,
                                   double /*t*/) {
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;

                auto x = point(e, q);

                value_type uu = eval(u, e, q, Ux, Uy);
                value_type uu_prev = eval(u_prev, e, q, Ux, Uy);
                value_type rr = eval(r.data, e, q, *r.Vx, *r.Vy);

                auto diff = diffusion(x[0], x[1]);
                auto bx = -x[1];
                auto by =  x[0];

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    // double F = erikkson_forcing(x[0], x[1], epsilon, t + steps.dt);
                    double F = 0;//source(x[0], x[1]);
                    double lv = (uu_prev.val + steps.dt * F) * v.val;

                    double val = -lv;

                    // Bu
                    val += uu.val * v.val + steps.dt * (diff * grad_dot(uu, v) + (bx * uu.dx + by * uu.dy) * v.val);
                    // -Aw
                    val -= (rr.val * v.val + h * h * (rr.dx * v.dx + rr.dy * v.dy));

                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);
                    double val = 0;

                    // -B'w
                    val -= w.val * rr.val + steps.dt * (diff * grad_dot(w, rr) + (bx * w.dx + by * w.dy) * rr.val);

                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
    }

    double errorL2(double /*t*/) const {
        auto sol = [](auto) { return value_type{0, 0, 0}; };
        return Base::errorL2(u, Ux, Uy, sol);

    }

    double errorH1(double /*t*/) const {
        auto sol = [](auto) { return value_type{0, 0, 0}; };
        return Base::errorH1(u, Ux, Uy, sol);

    }
};

}



#endif // ERIKKSON_POLLUTION_ROTATION_HPP
