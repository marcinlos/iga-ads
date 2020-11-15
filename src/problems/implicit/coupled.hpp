#ifndef PROBLEMS_IMPLICIT_HPP_
#define PROBLEMS_IMPLICIT_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class coupled : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;

    std::vector<double> sol, buf; // Buffers for keeping the combined RHS vectors

    // Two sets of vectors for two coupled equations
    vector_type u, u_prev;
    vector_type u2, u2_prev;


    output_manager<2> output;
    galois_executor executor{4};

    // Large matrices - dense matrices + data for solver
    lin::dense_matrix Ax, Ay;
    lin::solver_ctx Ax_ctx, Ay_ctx;

    // NEW: Mass matrices for the periodic basis
    lin::dense_matrix Mx, My;
    lin::solver_ctx Mx_ctx, My_ctx;

    double s = 40;

public:
    coupled(const config_2d& config)
    : Base{ config }
    , sol(2 * (x.dofs() - 1) * (y.dofs() - 1))
    , buf(2 * (x.dofs() - 1) * (y.dofs() - 1))
    , u{ shape() }
    , u_prev{ shape() }
    , u2{ shape() }
    , u2_prev{ shape() }
    , output{ x.B, y.B, 200 }
    , Ax{ 2 * (x.dofs() - 1), 2 * (x.dofs() - 1) } // 2x2 matrix for x direction
    , Ay{ 2 * (y.dofs() - 1), 2 * (y.dofs() - 1) } // 2x2 matrix for y direction
    , Ax_ctx{ Ax }
    , Ay_ctx{ Ay }
    , Mx{ x.dofs() - 1, x.dofs() - 1 } // NEW: initialization of mass matrices
    , My{ y.dofs() - 1, y.dofs() - 1 }
    , Mx_ctx{ Mx }
    , My_ctx{ My }
    {
        // Fill the large matrices
        matrix(Ax, x.basis, steps.dt);
        matrix(Ay, y.basis, steps.dt);

        // NEW: Fill the mass matrices
        mass_matrix(Mx, x.basis, steps.dt);
        mass_matrix(My, y.basis, steps.dt);
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(20 * (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

private:
    // NEW: Mass matrix no longer diagonal
    void mass_matrix(lin::dense_matrix& M, const basis_data& d, double /*h*/) {
        auto N = d.basis.dofs() - 1;
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
                        M(ia % N, ib % N) += va * vb * d.w[q] * d.J[e]; // upper left
                    }
                }
            }
        }
    }

    void matrix(lin::dense_matrix& K, const basis_data& d, double h) {
        auto N = d.basis.dofs() - 1;
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
                        auto da = d.b[e][q][1][a];
                        auto db = d.b[e][q][1][b];
                        K(ia % N, ib % N) += (va * vb + 0.5 * h * da * db - 0.5 * h * s * va * db) * d.w[q] * d.J[e]; // upper left
                        K(N + ia % N, N + ib % N) += (va * vb + 0.5 * h * da * db - 0.5 * h * s * va * db) * d.w[q] * d.J[e]; // lower right

                        // Mass-Mass-Stiffness matrix (what you need if I recall correctly):
                        // K(ia, ib) += va * vb * d.w[q] * d.J[e]; // upper left
                        // K(N + ia, N + ib) += va * vb * d.w[q] * d.J[e]; // lower right
                        // K(N + ia, ib) += da * db * d.w[q] * d.J[e]; // lower left
                    }
                }
            }
        }
    }

    void prepare_matrices() {
        Base::prepare_matrices();
        // Dense factorization (overloaded function)
        lin::factorize(Ax, Ax_ctx);
        lin::factorize(Ay, Ay_ctx);
        // NEW: factorize mass matrices
        lin::factorize(Mx, Mx_ctx);
        lin::factorize(My, My_ctx);
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(2*x - 1, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "out1_0.data");


        auto init2 = [this](double x, double y) { return init_state(x, y - 0.3); };
        projection(u2, init2);
        solve(u2);
        output.to_file(u2, "out2_0.data");

    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
        swap(u2, u2_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs_1();

        // NEW: different size and indices are taken mod Nx, Ny to ensure periodic BC
        int Nx = x.dofs() - 1;
        int Ny = y.dofs() - 1;
        // Copy data from separate RHS vectors to the combined one
        std::fill(begin(sol), end(sol), 0);
        vector_view view_x{ sol.data(), {2*Nx, Ny}};

        for (auto i = 0; i < x.dofs(); ++ i) {
            for (auto j = 0; j < y.dofs(); ++ j) {
                view_x(i % Nx, j % Ny) += u(i, j);
                view_x(Nx + i % Nx, j % Ny) += u2(i, j);
            }
        }

        // ads_solve(u, buffer, dim_data{Kx, x.ctx}, y.data());
        // ads_solve(u2, buffer, dim_data{Kx, x.ctx}, y.data());
        // Expanded ADS call:
        lin::solve_with_factorized(Ax, view_x, Ax_ctx);
        auto F = lin::cyclic_transpose(view_x, buf.data());
        // NEW: periodic mass matrix instead of a standard one
        // lin::solve_with_factorized(y.data().M, F, y.data().ctx);
        lin::solve_with_factorized(My, F, My_ctx);
        lin::cyclic_transpose(F, view_x);

        // ... and copy the solutions back to the separate vectors
        for (auto i = 0; i < x.dofs(); ++ i) {
            for (auto j = 0; j < y.dofs(); ++ j) {
                u(i, j) = view_x(i % Nx, j % Ny);
                u2(i, j) = view_x(Nx + i % Nx, j % Ny);
            }
        }

        using std::swap;
        swap(u, u_prev);
        swap(u2, u2_prev);

        compute_rhs_2();

        // Copy data from separate RHS vectors to the combined one
        std::fill(begin(sol), end(sol), 0);
        vector_view view_y{ sol.data(), {Nx, 2*Ny}};
        for (auto i = 0; i < x.dofs(); ++ i) {
            for (auto j = 0; j < y.dofs(); ++ j) {
                view_y(i % Nx, j % Ny) += u(i, j);
                view_y(i % Nx, Ny + j % Ny) += u2(i, j);
            }
        }

        // ads_solve(u, buffer, x.data(), dim_data{Ky, y.ctx});
        // ads_solve(u2, buffer, x.data(), dim_data{Ky, y.ctx});
        // Expanded ADS call:
        // NEW: periodic mass matrix instead of a standard one
        // lin::solve_with_factorized(x.data().M, view_y, x.data().ctx);
        lin::solve_with_factorized(Mx, view_y, Mx_ctx);
        auto F2 = lin::cyclic_transpose(view_y, buf.data());
        lin::solve_with_factorized(Ay, F2, Ay_ctx);
        lin::cyclic_transpose(F2, view_y);

        // ... and copy the solutions back to the separate vectors
        for (auto i = 0; i < x.dofs(); ++ i) {
            for (auto j = 0; j < y.dofs(); ++ j) {
                u(i, j) = view_y(i % Nx, j % Ny);
                u2(i, j) = view_y(i % Nx, Ny + j % Ny);
            }
        }

    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % 1 == 0) {
            std::cout << "Iter " << iter << std::endl;
            output.to_file(u, "out1_%d.data", iter + 1);
            output.to_file(u2, "out2_%d.data", iter + 1);
        }
    }

    void compute_rhs_1() {
        auto& rhs = u;
        auto& rhs2 = u2;

        zero(rhs);
        zero(rhs2);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto U2 = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);
                    value_type u2 = eval_fun(u2_prev, e, q);

                    double gradient_prod = u.dy * v.dy;
                    double val = u.val * v.val - 0.5 * steps.dt * gradient_prod + 0.5 * steps.dt * s * u.dy * v.val;
                    U(aa[0], aa[1]) += val * w * J;

                    gradient_prod = u2.dy * v.dy;
                    val = u2.val * v.val - 0.5 * steps.dt * gradient_prod + 0.5 * steps.dt * s * u2.dy * v.val;
                    U2(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
                update_global_rhs(rhs2, U2, e);
            });
        });
    }

    void compute_rhs_2() {
        auto& rhs = u;
        auto& rhs2 = u2;

        zero(rhs);
        zero(rhs2);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto U2 = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);
                    value_type u2 = eval_fun(u2_prev, e, q);

                    double gradient_prod = u.dx * v.dx;
                    double val = u.val * v.val - 0.5 * steps.dt * gradient_prod + 0.5 * steps.dt * s * u.dx * v.val;
                    U(aa[0], aa[1]) += val * w * J;

                    gradient_prod = u2.dx * v.dx;
                    val = u2.val * v.val - 0.5 * steps.dt * gradient_prod + 0.5 * steps.dt * s * u2.dx * v.val;
                    U2(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
                update_global_rhs(rhs2, U2, e);
            });
        });
    }

};

}

#endif /* PROBLEMS_IMPLICIT_HPP_*/
