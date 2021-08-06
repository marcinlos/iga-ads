// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_GALERKIN_HPP
#define MAXWELL_MAXWELL_GALERKIN_HPP

#include <iostream>
#include <utility>

#include "ads/form_matrix.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "maxwell_base.hpp"
#include "problems.hpp"
#include "spaces.hpp"
#include "state.hpp"

class maxwell_galerkin : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = maxwell_manufactured1;

    space V;
    space_set U;

    state prev, half, now;

    ads::lin::band_matrix Bx, By, Bz;
    ads::lin::solver_ctx Bx_ctx, By_ctx, Bz_ctx;

    Problem problem{1, 1};

    ads::output_manager<3> output;

public:
    explicit maxwell_galerkin(const ads::config_3d& config)
    : Base{config}
    , V{x, y, z}
    , U{V, V, V, V, V, V}
    , prev{vector_shape(V)}
    , half{vector_shape(V)}
    , now{vector_shape(V)}
    , Bx{V.x.p, V.x.p, V.x.dofs()}
    , By{V.y.p, V.y.p, V.y.dofs()}
    , Bz{V.z.p, V.z.p, V.z.dofs()}
    , Bx_ctx{Bx}
    , By_ctx{By}
    , Bz_ctx{Bz}
    , output{V.x.B, V.y.B, V.z.B, 50} { }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

private:
    void prepare_matrices() {
        set_boundary_conditions(U);

        factorize_matrices(U);
        factorize_matrices(V);

        Bx.zero();
        By.zero();
        Bz.zero();

        auto tau = steps.dt;
        auto h = tau * tau / (4 * problem.eps({0.5, 0.5}) * problem.mu({0.5, 0.5}));
        auto form = [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; };

        form_matrix(Bx, V.x.basis, form);
        form_matrix(By, V.y.basis, form);
        form_matrix(Bz, V.z.basis, form);

        fix_dof(0, V.x, Bx);
        fix_dof(V.x.dofs() - 1, V.x, Bx);
        fix_dof(0, V.y, By);
        fix_dof(V.y.dofs() - 1, V.y, By);
        fix_dof(0, V.z, Bz);
        fix_dof(V.z.dofs() - 1, V.z, Bz);

        ads::lin::factorize(Bx, Bx_ctx);
        ads::lin::factorize(By, By_ctx);
        ads::lin::factorize(Bz, Bz_ctx);
    }

    void before() override {
        prepare_matrices();
        set_init_state(now, U, problem);
        after_step(-1, -steps.dt);
    }

    auto substep1_solve_E(state& rhs, vector_type& buffer) -> void {
        ads_solve(rhs.E1, buffer, U.E1.x.data(), ads::dim_data{By, By_ctx}, U.E1.z.data());
        ads_solve(rhs.E2, buffer, U.E2.x.data(), U.E2.y.data(), ads::dim_data{Bz, Bz_ctx});
        ads_solve(rhs.E3, buffer, ads::dim_data{Bx, Bx_ctx}, U.E3.y.data(), U.E3.z.data());
    }

    auto substep2_solve_E(state& rhs, vector_type& buffer) -> void {
        ads_solve(rhs.E1, buffer, U.E1.x.data(), U.E1.y.data(), ads::dim_data{Bz, Bz_ctx});
        ads_solve(rhs.E2, buffer, ads::dim_data{Bx, Bx_ctx}, U.E2.y.data(), U.E2.z.data());
        ads_solve(rhs.E3, buffer, U.E3.x.data(), ads::dim_data{By, By_ctx}, U.E3.z.data());
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
        substep1_solve_E(mid, buffer);

        substep1_fill_H(mid, prev, mid, U, c);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        substep2_solve_E(now, buffer);

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

#endif  // MAXWELL_MAXWELL_GALERKIN_HPP
