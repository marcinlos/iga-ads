// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_CAUCHY_HEAD_HPP
#define MAXWELL_MAXWELL_CAUCHY_HEAD_HPP

#include <iostream>
#include <utility>

#include "ads/experimental/all.hpp"
#include "ads/form_matrix.hpp"
#include "ads/simulation.hpp"
#include "antenna.hpp"
#include "head_antenna_problem.hpp"
#include "maxwell_base.hpp"
#include "solver_phase.hpp"
#include "spaces.hpp"
#include "state.hpp"
#include "utils.hpp"

class maxwell_cauchy_head : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = head_antenna_problem;

    space V;
    space_set U;

    ads::regular_mesh3 mesh_;
    ads::quadrature3 quad_;
    ads::space3 space_;

    state pprev, prev, half, now;

    ads::lin::tensor<double, 3> stiffness_coeffs;
    solver_phase Bx, By, Bz;

    Problem problem;
    antenna source;

public:
    explicit maxwell_cauchy_head(ads::config_3d const& config, ads::regular_mesh3 const& mesh,
                                 ads::quadrature3 const& quad, ads::space3& space,
                                 std::string_view data_file)
    : Base{config}
    , V{x, y, z}
    , U{V, V, V, V, V, V}
    , mesh_{mesh}
    , quad_{quad}
    , space_{space}
    , pprev{vector_shape(V)}
    , prev{vector_shape(V)}
    , half{vector_shape(V)}
    , now{vector_shape(V)}
    , stiffness_coeffs{vector_shape(V)}
    , Bx{V.y.dofs(), V.z.dofs()}
    , By{V.z.dofs(), V.y.dofs()}
    , Bz{V.x.dofs(), V.y.dofs()}
    , problem{data_file}
    , source{problem.omega, problem.tau} { }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(prev, pprev);
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

    void step(int /*iter*/, double t) override {
        const auto tau = steps.dt;
        const auto a = [this, tau](auto x) { return tau / (2 * problem.eps(x)); };
        const auto b = [this, tau](auto x) {
            return tau * tau / (4 * problem.eps(x) * problem.mu(x));
        };
        const auto c = [this, tau](auto x) { return tau / (2 * problem.mu(x)); };

        const auto shape = vector_shape(V);

        // Buffer large enough for all the RHS
        auto buffer = vector_type{shape};

        auto mid = state{shape};

        // First substep
        substep1_fill_E(mid, prev, U, a, b);
        source.apply_forcing(t, mid, V);
        substep1_boundary_E(t, mid, prev, pprev);
        substep1_solve_E(mid, buffer);

        substep1_fill_H(mid, prev, mid, U, c);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        source.apply_forcing(t, now, V);
        substep2_boundary_E(t, now, mid, prev);
        substep2_solve_E(now, buffer);

        substep2_fill_H(now, mid, now, U, c);
        solve_H(now, buffer);

    }

    auto substep1_boundary_E(double t, state& rhs, state& prev, state& pprev) -> void {
        const auto tau = steps.dt;
        const auto a = [this, tau](auto x) {
            return -tau * tau / (4 * problem.eps(x) * problem.mu(x));
        };
        const auto mu = [this](auto x) { return problem.mu(x); };
        const auto b = [this](auto x) { return problem.mu(x) / problem.eta; };

        auto out = [](auto& buf) { return [&buf](int J, double val) { buf.data()[J] += val; }; };

        auto E1_prev = ads::bspline_function3(&space_, prev.E1.data());
        auto E2_prev = ads::bspline_function3(&space_, prev.E2.data());
        auto E3_prev = ads::bspline_function3(&space_, prev.E3.data());

        auto E1_pprev = ads::bspline_function3(&space_, pprev.E1.data());
        auto E2_pprev = ads::bspline_function3(&space_, pprev.E2.data());
        auto E3_pprev = ads::bspline_function3(&space_, pprev.E3.data());

        auto dE1_dt = [&](auto x) { return (E1_prev(x) - E1_pprev(x)) / tau; };
        auto dE2_dt = [&](auto x) { return (E2_prev(x) - E2_pprev(x)) / tau; };
        auto dE3_dt = [&](auto x) { return (E3_prev(x) - E3_pprev(x)) / tau; };

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E1),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<1>(face.normal));
                         auto const f = mu(x) * problem.U1(x, t) + b(x) * dE1_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<2>(face.normal));
                         auto const f = mu(x) * problem.U2(x, t) + b(x) * dE2_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<0>(face.normal));
                         auto const f = mu(x) * problem.U3(x, t) + b(x) * dE3_dt(xx);
                         return ok * a(x) * f * v.val;
                     });
    }

    auto substep2_boundary_E(double t, state& rhs, state& prev, state& pprev) -> void {
        const auto tau = steps.dt;
        const auto a = [this, tau](auto x) {
            return -tau * tau / (4 * problem.eps(x) * problem.mu(x));
        };
        const auto mu = [this](auto x) { return problem.mu(x); };
        const auto b = [this](auto x) { return problem.mu(x) / problem.eta; };

        auto out = [](auto& buf) { return [&buf](int J, double val) { buf.data()[J] += val; }; };

        auto E1_prev = ads::bspline_function3(&space_, prev.E1.data());
        auto E2_prev = ads::bspline_function3(&space_, prev.E2.data());
        auto E3_prev = ads::bspline_function3(&space_, prev.E3.data());

        auto E1_pprev = ads::bspline_function3(&space_, pprev.E1.data());
        auto E2_pprev = ads::bspline_function3(&space_, pprev.E2.data());
        auto E3_pprev = ads::bspline_function3(&space_, pprev.E3.data());

        auto dE1_dt = [&](auto x) { return (E1_prev(x) - E1_pprev(x)) / (tau / 2); };
        auto dE2_dt = [&](auto x) { return (E2_prev(x) - E2_pprev(x)) / (tau / 2); };
        auto dE3_dt = [&](auto x) { return (E3_prev(x) - E3_pprev(x)) / (tau / 2); };

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E1),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<2>(face.normal));
                         auto const f = mu(x) * problem.U1(x, t) + b(x) * dE1_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<0>(face.normal));
                         auto const f = mu(x) * problem.U2(x, t) + b(x) * dE2_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const ok = std::abs(std::get<1>(face.normal));
                         auto const f = mu(x) * problem.U3(x, t) + b(x) * dE3_dt(xx);
                         return ok * a(x) * f * v.val;
                     });
    }

    auto as_array(ads::point3_t x) const -> point_type {
        return {std::get<0>(x), std::get<1>(x), std::get<2>(x)};
    }

    void after_step(int iter, double t) override {
        const auto i = iter + 1;
        const auto tt = t + steps.dt;

        if (i % 1 == 0) {
            save(i);
        }

        auto const res = compute_norms(now, U, problem, tt);
        auto const q_sar = compute_q_sar(now, U);
        std::cout << "After step " << i << ", t = " << tt << '\n';
        std::cout << "q_sar = " << q_sar << '\n';
        print_result_info(res);
    }

    auto save(int iter) -> void {
        const auto name = fmt::format("out_{}.vti", iter);

        auto E1 = ads::bspline_function3(&space_, now.E1.data());
        auto E2 = ads::bspline_function3(&space_, now.E2.data());
        auto E3 = ads::bspline_function3(&space_, now.E3.data());

        auto H1 = ads::bspline_function3(&space_, now.H1.data());
        auto H2 = ads::bspline_function3(&space_, now.H2.data());
        auto H3 = ads::bspline_function3(&space_, now.H3.data());

        maxwell_to_file(name, E1, E2, E3, H1, H2, H3);
    }

    auto compute_q_sar(state const& s, space_set const& U) const -> double {
        double val = 0;

        for (auto const e : elements(U.E1.x, U.E1.y, U.E1.z)) {
            double const J = jacobian(e, U.E1.x, U.E1.y, U.E1.z);
            for (auto const q : quad_points(U.E1.x, U.E1.y, U.E1.z)) {
                double const w = weight(q, U.E1.x, U.E1.y, U.E1.z);
                auto const x = point(e, q, U.E1.x, U.E1.y, U.E1.z);
                auto const sigma = problem.empty(x) ? 0 : problem.sigma(x);

                auto const E1 = eval(s.E1, e, q, U.E1.x, U.E1.y, U.E1.z);
                auto const E2 = eval(s.E2, e, q, U.E2.x, U.E2.y, U.E2.z);
                auto const E3 = eval(s.E3, e, q, U.E3.x, U.E3.y, U.E1.z);

                auto const EE = E1.val * E1.val + E2.val * E2.val + E3.val * E3.val;
                val += 0.5 * EE * sigma * w * J;
            }
        }
        return val;
    }
};

#endif  // MAXWELL_MAXWELL_CAUCHY_HEAD_HPP
