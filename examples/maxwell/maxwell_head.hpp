// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_HEAD_HPP
#define MAXWELL_MAXWELL_HEAD_HPP

#include <iostream>
#include <utility>

#include "ads/form_matrix.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"
#include "maxwell_base.hpp"
#include "maxwell_head_problem.hpp"
#include "spaces.hpp"
#include "state.hpp"

class maxwell_head : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = maxwell_head_problem;

    space V;
    space_set U;

    ads::mumps::problem E1_1, E2_1, E3_1;
    ads::mumps::problem E1_2, E2_2, E3_2;

    state prev, half, now;

    Problem problem;
    ads::mumps::solver solver;

    bool avg_material_data;

    ads::output_manager<3> output;

public:
    maxwell_head(const ads::config_3d& config, std::string_view data_file, bool avg_material_data)
    : Base{config}
    , V{x, y, z}
    , U{V, V, V, V, V, V}
    , E1_1{nullptr, V.dofs()}
    , E2_1{nullptr, V.dofs()}
    , E3_1{nullptr, V.dofs()}
    , E1_2{nullptr, V.dofs()}
    , E2_2{nullptr, V.dofs()}
    , E3_2{nullptr, V.dofs()}
    , prev{vector_shape(V)}
    , half{vector_shape(V)}
    , now{vector_shape(V)}
    , problem{data_file}
    , avg_material_data{avg_material_data}
    , output{V.x.B, V.y.B, V.z.B, 50} { }

private:
    template <typename BC>
    void matrix(ads::mumps::problem& A, double tau, double cx, double cy, double cz, BC&& bc) {
        auto shape_loc = local_matrix_shape(V);

        for (auto e : elements(V.x, V.y, V.z)) {
            auto J = jacobian(e, V.x, V.y, V.z);
            using buffer_type = ads::lin::tensor<double, 6>;
            auto loc = buffer_type{shape_loc};

            for (auto q : quad_points(V.x, V.y, V.z)) {
                auto W = weight(q, V.x, V.y, V.z);
                auto x = point(e, q, V.x, V.y, V.z);

                for (auto i : dofs_on_element(e, V.x, V.y, V.z)) {
                    auto il = dof_global_to_local(e, i, V.x, V.y, V.z);
                    auto v = eval_basis(e, q, i, V.x, V.y, V.z);

                    double a = stiffness_coeff(tau, x, i);

                    for (auto j : dofs_on_element(e, V.x, V.y, V.z)) {
                        auto jl = dof_global_to_local(e, j, V.x, V.y, V.z);
                        auto u = eval_basis(e, q, j, V.x, V.y, V.z);

                        auto M = u.val * v.val;
                        auto S = cx * u.dx * v.dx + cy * u.dy * v.dy + cz * u.dz * v.dz;
                        auto form = M + a * S;
                        loc(il[0], il[1], il[2], jl[0], jl[1], jl[2]) += form * W * J;
                    }
                }
            }
            // Update global matrix
            for (auto i : dofs_on_element(e, V.x, V.y, V.z)) {
                if (bc(i))
                    continue;
                auto il = dof_global_to_local(e, i, V.x, V.y, V.z);
                int ii = linear_index(i, V.x, V.y, V.z) + 1;
                for (auto j : dofs_on_element(e, V.x, V.y, V.z)) {
                    if (bc(j))
                        continue;
                    auto jl = dof_global_to_local(e, j, V.x, V.y, V.z);
                    int jj = linear_index(j, V.x, V.y, V.z) + 1;

                    auto val = loc(il[0], il[1], il[2], jl[0], jl[1], jl[2]);
                    A.add(ii, jj, val);
                }
            }
        }
        apply_dirichlet_bc(A, bc);
    }

    auto stiffness_coeff(double tau, point_type x, index_type test) const -> double {
        point_type xx;
        if (avg_material_data) {
            auto const [rx, ry, rz] = dof_support(test, V);
            xx = point_type{(rx.a + rx.b) / 2, (ry.a + ry.b) / 2, (rz.a + rz.b) / 2};
        } else {
            xx = x;
        }
        double eps = problem.eps(xx);
        double mu = problem.mu(xx);
        return tau * tau / (4 * eps * mu);
    }

    template <typename BC>
    auto apply_dirichlet_bc(ads::mumps::problem& A, BC&& bc) const -> void {
        for (auto i : dofs(V.x, V.y, V.z)) {
            if (bc(i)) {
                int ii = linear_index(i, V.x, V.y, V.z) + 1;
                A.add(ii, ii, 1.0);
            }
        }
    }

    void prepare_matrices() {
        set_boundary_conditions(U);

        factorize_matrices(U);
        factorize_matrices(V);

        auto tau = steps.dt;

        // E x n = 0
        matrix(E1_1, tau, 0, 1, 0,
               [this](auto i) { return is_boundary(i[1], V.y) || is_boundary(i[2], V.z); });
        matrix(E2_1, tau, 0, 0, 1,
               [this](auto i) { return is_boundary(i[0], V.x) || is_boundary(i[2], V.z); });
        matrix(E3_1, tau, 1, 0, 0,
               [this](auto i) { return is_boundary(i[0], V.x) || is_boundary(i[1], V.y); });

        matrix(E1_2, tau, 0, 0, 1,
               [this](auto i) { return is_boundary(i[1], V.y) || is_boundary(i[2], V.z); });
        matrix(E2_2, tau, 1, 0, 0,
               [this](auto i) { return is_boundary(i[0], V.x) || is_boundary(i[2], V.z); });
        matrix(E3_2, tau, 0, 1, 0,
               [this](auto i) { return is_boundary(i[0], V.x) || is_boundary(i[1], V.y); });
    }

    void before() override {
        prepare_matrices();
        set_init_state(now, U, problem);
        after_step(-1, -steps.dt);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

    auto substep1_solve_E(state& rhs) -> void {
        E1_1.rhs(rhs.E1.data());
        solver.solve(E1_1);

        E2_1.rhs(rhs.E2.data());
        solver.solve(E2_1);

        E3_1.rhs(rhs.E3.data());
        solver.solve(E3_1);
    }

    auto substep2_solve_E(state& rhs) -> void {
        E1_2.rhs(rhs.E1.data());
        solver.solve(E1_2);

        E2_2.rhs(rhs.E2.data());
        solver.solve(E2_2);

        E3_2.rhs(rhs.E3.data());
        solver.solve(E3_2);
    }

    auto solve_H(state& rhs, vector_type& buffer) -> void {
        ads_solve(rhs.H1, buffer, U.H1.x.data(), U.H1.y.data(), U.H1.z.data());
        ads_solve(rhs.H2, buffer, U.H2.x.data(), U.H2.y.data(), U.H2.z.data());
        ads_solve(rhs.H3, buffer, U.H3.x.data(), U.H3.y.data(), U.H3.z.data());
    }

    void step(int /*iter*/, double /*t*/) override {
        const auto tau = steps.dt;
        const auto a = [this, tau](auto x) { return tau / (2 * problem.eps(x)); };
        const auto b = [this, tau](auto x) { return tau * tau / (4 * problem.eps(x)); };
        const auto c = [this, tau](auto x) { return tau / (2 * problem.mu(x)); };

        const auto shape = vector_shape(V);

        // Buffer large enough for all the RHS
        auto buffer = vector_type{shape};

        auto mid = state{shape};

        // First substep
        substep1_fill_E(mid, prev, U, a, b);
        substep1_solve_E(mid);

        substep1_fill_H(mid, prev, mid, U, c);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        substep2_solve_E(now);

        substep2_fill_H(now, mid, now, U, c);
        solve_H(now, buffer);
    }

    void after_step(int iter, double t) override {
        const auto i = iter + 1;
        const auto tt = t + steps.dt;

        if (i % 10 == 0)
            output_solution(output, i, now);

        const auto res = compute_norms(now, U, problem, tt);
        std::cout << "After step " << i << ", t = " << tt << '\n';
        print_result_info(res);
    }
};

#endif  // MAXWELL_MAXWELL_HEAD_HPP
