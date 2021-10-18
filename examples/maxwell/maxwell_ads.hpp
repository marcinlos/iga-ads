// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_ADS_HPP
#define MAXWELL_MAXWELL_ADS_HPP

#include <cassert>
#include <string_view>
#include <utility>

#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "maxwell_base.hpp"
#include "maxwell_head_problem.hpp"
#include "spaces.hpp"
#include "state.hpp"

class solver_phase {
private:
    using matrix = ads::lin::band_matrix;
    using context = ads::lin::solver_ctx;

    struct entry {
        matrix mat;
        context ctx;

        entry()
        : ctx{mat} { }

        explicit entry(matrix mat)
        : mat{std::move(mat)}
        , ctx{this->mat} { }
    };

    using entry_table = ads::lin::tensor<entry, 2>;

    entry_table entries_;

public:
    solver_phase(int n, int m)
    : entries_{{n, m}} { }

    auto set_matrix(int i, int j, matrix mat) -> void {
        assert(i < entries_.size(0));
        assert(j < entries_.size(1));

        auto e = entry{std::move(mat)};
        ads::lin::factorize(e.mat, e.ctx);

        entries_(i, j) = std::move(e);
    }

    template <typename Rhs>
    auto operator()(Rhs& rhs) -> void {
        auto* data = rhs.data();
        auto const rhs_size = rhs.size(0);

        auto offset = 0;
        for (int j = 0; j < entries_.size(1); ++j) {
            for (int i = 0; i < entries_.size(0); ++i) {
                auto& entry = entries_(i, j);
                solve_with_factorized(entry.mat, data + offset, entry.ctx, 1);
                offset += rhs_size;
            }
        }
    }
};

class maxwell_ads : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = maxwell_head_problem;

    space V;
    space_set U;

    state prev, half, now;

    ads::lin::tensor<double, 3> stiffness_coeffs;
    solver_phase Bx, By, Bz;

    Problem problem;

    ads::output_manager<3> output;

public:
    maxwell_ads(ads::config_3d const& config, std::string_view data_file)
    : Base{config}
    , V{x, y, z}
    , U{V, V, V, V, V, V}
    , prev{vector_shape(V)}
    , half{vector_shape(V)}
    , now{vector_shape(V)}
    , stiffness_coeffs{vector_shape(V)}
    , Bx{V.y.dofs(), V.z.dofs()}
    , By{V.z.dofs(), V.y.dofs()}
    , Bz{V.x.dofs(), V.y.dofs()}
    , problem{data_file}
    , output{V.x.B, V.y.B, V.z.B, 50} { }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

private:
    auto compute_stiffness_coeffs() -> void {
        auto const tau = steps.dt;
        for (auto const i : dofs(V.x, V.y, V.z)) {
            auto const [rx, ry, rz] = dof_support(i, V);
            auto const x = point_type{(rx.a + rx.b) / 2, (ry.a + ry.b) / 2, (rz.a + rz.b) / 2};
            auto const eps = problem.eps(x);
            auto const mu = problem.mu(x);
            stiffness_coeffs(i[0], i[1], i[2]) = tau * tau / (4 * eps * mu);
        }
    }

    template <typename Form>
    void form_matrix(ads::lin::band_matrix& M, ads::basis_data const& d, Form&& form) {
        for (auto e = 0; e < d.elements; ++e) {
            for (auto q = 0; q < d.quad_order; ++q) {
                auto const first = d.first_dof(e);
                auto const last = d.last_dof(e);
                for (auto a = 0; a + first <= last; ++a) {
                    for (auto b = 0; b + first <= last; ++b) {
                        int const ia = a + first;
                        int const ib = b + first;
                        auto const va = d.b[e][q][0][a];
                        auto const vb = d.b[e][q][0][b];
                        auto const da = d.b[e][q][1][a];
                        auto const db = d.b[e][q][1][b];
                        auto const fa = ads::function_value_1d{va, da};
                        auto const fb = ads::function_value_1d{vb, db};
                        M(ia, ib) += form(fb, fa, ia) * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    template <typename CoeffFun>
    auto fill_solver_phase(solver_phase& phase, ads::dimension const& special_dim,
                           ads::dimension const& dim2, ads::dimension const& dim3, CoeffFun&& coeff)
        -> void {
        for (int i = 0; i < dim2.dofs(); ++i) {
            for (int j = 0; j < dim3.dofs(); ++j) {
                auto M = ads::lin::band_matrix{special_dim.p, special_dim.p, special_dim.dofs()};

                auto const form = [i, j, &coeff](auto u, auto v, auto dof) {
                    auto const h = coeff(dof, i, j);
                    return u.val * v.val + h * u.dx * v.dx;
                };
                form_matrix(M, special_dim.basis, form);

                fix_dof(0, special_dim, M);
                fix_dof(special_dim.dofs() - 1, special_dim, M);

                phase.set_matrix(i, j, M);
            }
        }
    }

    auto fill_solver_phases() -> void {
        fill_solver_phase(Bx, V.x, V.y, V.z,
                          [this](int ix, int iy, int iz) { return stiffness_coeffs(ix, iy, iz); });
        fill_solver_phase(By, V.y, V.z, V.x,
                          [this](int iy, int iz, int ix) { return stiffness_coeffs(ix, iy, iz); });
        fill_solver_phase(Bz, V.z, V.x, V.y,
                          [this](int iz, int ix, int iy) { return stiffness_coeffs(ix, iy, iz); });
    }

    void prepare_matrices() {
        set_boundary_conditions(U);

        factorize_matrices(U);
        factorize_matrices(V);

        compute_stiffness_coeffs();
        fill_solver_phases();
    }

    void before() override {
        prepare_matrices();
        set_init_state(now, U, problem);
        after_step(-1, -steps.dt);
    }

    auto substep1_solve_E(state& rhs, vector_type& buffer) -> void {
        ads_solve(rhs.E1, buffer, U.E1.x.data(), By, U.E1.z.data());
        ads_solve(rhs.E2, buffer, U.E2.x.data(), U.E2.y.data(), Bz);
        ads_solve(rhs.E3, buffer, Bx, U.E3.y.data(), U.E3.z.data());
    }

    auto substep2_solve_E(state& rhs, vector_type& buffer) -> void {
        ads_solve(rhs.E1, buffer, U.E1.x.data(), U.E1.y.data(), Bz);
        ads_solve(rhs.E2, buffer, Bx, U.E2.y.data(), U.E2.z.data());
        ads_solve(rhs.E3, buffer, U.E3.x.data(), By, U.E3.z.data());
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
        apply_bc_E(mid, U);
        substep1_solve_E(mid, buffer);

        substep1_fill_H(mid, prev, mid, U, c);
        apply_bc_H(mid, U);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        apply_bc_E(now, U);
        substep2_solve_E(now, buffer);

        substep2_fill_H(now, mid, now, U, c);
        apply_bc_H(now, U);
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

#endif  // MAXWELL_MAXWELL_ADS_HPP
