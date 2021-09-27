// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ERIKKSON_ERIKKSON_INVERSE_HPP
#define ERIKKSON_ERIKKSON_INVERSE_HPP

#include <galois/Timer.h>

#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/solver/mumps.hpp"
#include "erikkson_base.hpp"

namespace ads {

class erikkson_inverse : public erikkson_base {
private:
    using Base = erikkson_base;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;

    lin::band_matrix MVx, MVy;
    lin::band_matrix KVx, KVy;

    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;

    vector_type u, rhs;
    vector_type u_buffer;
    std::vector<double> full_rhs;

    int save_every = 1;

    double h;

    double peclet = 1e2;
    double epsilon = 1 / peclet;

    point_type c_diff{{epsilon, epsilon}};

    point_type beta;

    galois::StatTimer solver_timer{"solver"};
    mumps::solver solver;

    output_manager<2> output;

public:
    erikkson_inverse(double bx, double by,                                //
                     dimension const& trial_x, dimension const& trial_y,  //
                     dimension const& test_x, dimension const& test_y)
    : Base{test_x, test_y, timesteps_config{1, 0.0}}
    , Ux{trial_x}
    , Uy{trial_y}
    , Vx{x}
    , Vy{y}
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUVx{Vx.dofs(), Ux.dofs()}
    , MUVy{Vy.dofs(), Uy.dofs()}
    , KUVx{Vx.dofs(), Ux.dofs()}
    , KUVy{Vy.dofs(), Uy.dofs()}
    , AUVx{Vx.dofs(), Ux.dofs()}
    , AUVy{Vy.dofs(), Uy.dofs()}
    , u{{Ux.dofs(), Uy.dofs()}}
    , rhs{{Vx.dofs(), Vy.dofs()}}
    , u_buffer{{Ux.dofs(), Uy.dofs()}}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , h{element_diam(Ux, Uy)}
    , beta{bx, by}
    , output{Ux.B, Uy.B, 500} { }

private:
    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    void assemble_problem(mumps::problem& problem, double /*dt*/) {
        auto N = Vx.dofs() * Vy.dofs();
        auto hh = h * h;

        // Gram matrix - upper left
        for (auto i : internal_dofs(Vx, Vy)) {
            for (auto j : overlapping_internal_dofs(i, Vx, Vy)) {
                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double val =
                    kron(MVx, MVy, i, j) + hh * (kron(KVx, MVy, i, j) + kron(MVx, KVy, i, j));
                problem.add(ii, jj, val);
            }
        }

        // B, B^T
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : dofs(Ux, Uy)) {
                double val = 0;
                val += c_diff[0] * kron(KUVx, MUVy, i, j) + beta[0] * kron(AUVx, MUVy, i, j);
                val += c_diff[1] * kron(MUVx, KUVy, i, j) + beta[1] * kron(MUVx, AUVy, i, j);

                if (val != 0) {
                    int ii = linear_index(i, Vx, Vy) + 1;
                    int jj = linear_index(j, Ux, Uy) + 1;

                    if (!is_boundary(i, Vx, Vy)) {
                        problem.add(ii, N + jj, -val);
                    }
                    if (!is_boundary(i, Vx, Vy) && !is_boundary(j, Ux, Uy)) {
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

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);
        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);

        gram_matrix_1d(MUVx, Ux.basis, Vx.basis);
        gram_matrix_1d(MUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KUVx, Ux.basis, Vx.basis);
        stiffness_matrix_1d(KUVy, Uy.basis, Vy.basis);

        advection_matrix_1d(AUVx, Ux.basis, Vx.basis);
        advection_matrix_1d(AUVy, Uy.basis, Vy.basis);
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min((dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();

        stationary_bc(u, Ux, Uy);
    }

    void step(int /*iter*/, double t) override {
        mumps::problem problem{full_rhs};

        std::cout << "Assembling matrix" << std::endl;
        assemble_problem(problem, steps.dt);

        std::cout << "Computing RHS" << std::endl;
        zero(rhs);

        std::fill(begin(full_rhs), end(full_rhs), 0);
        vector_view view_in{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view view_out{full_rhs.data() + view_in.size(), {Ux.dofs(), Uy.dofs()}};

        for (auto i : dofs(Vx, Vy)) {
            view_in(i[0], i[1]) = -rhs(i[0], i[1]);
        }
        stationary_bc(view_out, Ux, Uy);

        std::cout << "Solving" << std::endl;
        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) = view_out(i[0], i[1]);
        }

        std::cout << "  solver time:       " << static_cast<double>(solver_timer.get()) << " ms"
                  << std::endl;
        std::cout << "  assembly    FLOPS: " << solver.flops_assembly() << std::endl;
        std::cout << "  elimination FLOPS: " << solver.flops_elimination() << std::endl;

        std::cout << "Error: L2 = " << errorL2(t) << "%, H1 =  " << errorH1(t) << "%" << std::endl;
    }

    void after() override {
        output.to_file(u, "values.data");
        print_solution("coeffs.data", u, Ux, Uy);
    }

    double errorL2(double /*t*/) const {
        auto const exact = [&](auto x) { return erikkson_exact(x[0], x[1], epsilon); };
        return error_relative_L2(u, Ux, Uy, exact) * 100;
    }

    double errorH1(double /*t*/) const {
        auto const exact = [&](auto x) { return erikkson_exact(x[0], x[1], epsilon); };
        return error_relative_H1(u, Ux, Uy, exact) * 100;
    }
};

}  // namespace ads

#endif  // ERIKKSON_ERIKKSON_INVERSE_HPP
