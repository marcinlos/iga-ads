// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_CAUCHY_HPP
#define MAXWELL_MAXWELL_CAUCHY_HPP

#include <iostream>
#include <utility>

#include "ads/experimental/all.hpp"
#include "ads/form_matrix.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "maxwell_base.hpp"
#include "scattering_problem.hpp"
#include "spaces.hpp"
#include "state.hpp"

class maxwell_cauchy : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = scattering_problem;

    space V;
    space_set U;

    ads::regular_mesh3 mesh_;
    ads::quadrature3 quad_;
    ads::space3 space_;

    state pprev, prev, half, now;

    ads::lin::band_matrix Bx, By, Bz, Bz_E1;
    ads::lin::solver_ctx Bx_ctx, By_ctx, Bz_ctx, Bz_E1_ctx;

    Problem problem;

    ads::output_manager<3> output;

public:
    explicit maxwell_cauchy(ads::config_3d const& config, ads::regular_mesh3 const& mesh,
                            ads::quadrature3 const& quad, ads::space3& space)
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
    , Bx{V.x.p, V.x.p, V.x.dofs()}
    , By{V.y.p, V.y.p, V.y.dofs()}
    , Bz{V.z.p, V.z.p, V.z.dofs()}
    , Bz_E1{V.z.p, V.z.p, V.z.dofs()}
    , Bx_ctx{Bx}
    , By_ctx{By}
    , Bz_ctx{Bz}
    , Bz_E1_ctx{Bz_E1}
    , output{V.x.B, V.y.B, V.z.B, 2, 2, 400} { }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(prev, pprev);
        swap(now, prev);
    }

private:
    void prepare_matrices() {
        auto const zero = [](auto& dim) {
            dim.fix_left();
            dim.fix_right();
        };
        zero(U.E2.x);
        zero(U.E3.x);

        U.E1.z.fix_right();

        factorize_matrices(U);
        factorize_matrices(V);

        auto tau = steps.dt;
        auto h = tau * tau / (4 * problem.eps({0.5, 0.5}) * problem.mu({0.5, 0.5}));
        auto form = [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; };

        form_matrix(Bx, V.x.basis, form);
        form_matrix(By, V.y.basis, form);
        form_matrix(Bz, V.z.basis, form);
        form_matrix(Bz_E1, V.z.basis, form);

        fix_dof(0, V.x, Bx);
        fix_dof(V.x.dofs() - 1, V.x, Bx);

        fix_dof(V.z.dofs() - 1, V.z, Bz_E1);

        ads::lin::factorize(Bx, Bx_ctx);
        ads::lin::factorize(By, By_ctx);
        ads::lin::factorize(Bz, Bz_ctx);
        ads::lin::factorize(Bz_E1, Bz_E1_ctx);
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
        ads_solve(rhs.E1, buffer, U.E1.x.data(), U.E1.y.data(), ads::dim_data{Bz_E1, Bz_E1_ctx});
        ads_solve(rhs.E2, buffer, ads::dim_data{Bx, Bx_ctx}, U.E2.y.data(), U.E2.z.data());
        ads_solve(rhs.E3, buffer, U.E3.x.data(), ads::dim_data{By, By_ctx}, U.E3.z.data());
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
        substep1_boundary_E(t, mid, prev, pprev);
        zero_sides("x", mid.E2, U.E2);
        zero_sides("x", mid.E3, U.E3);
        zero_z_up(mid.E1, U.E1);
        substep1_solve_E(mid, buffer);

        substep1_fill_H(mid, prev, mid, U, c);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        substep2_boundary_E(t, now, mid, prev);
        zero_sides("x", now.E2, U.E2);
        zero_sides("x", now.E3, U.E3);
        zero_z_up(now.E1, U.E1);
        substep2_solve_E(now, buffer);

        substep2_fill_H(now, mid, now, U, c);
        solve_H(now, buffer);

        // zero(now.E1);
        // zero(now.E2);
        // zero(now.E3);
        // integrate_U(t, now);

        // zero(now.E1);
        // zero(now.E2);
        // zero(now.E3);
        // zero(now.H1);
        // zero(now.H2);
        // zero(now.H3);

        // project(now.E1, U.E1, problem.E1_val_at(t));
        // project(now.E2, U.E2, problem.E2_val_at(t));
        // project(now.E3, U.E3, problem.E3_val_at(t));

        // project(now.H1, U.H1, problem.H1_val_at(t));
        // project(now.H2, U.H2, problem.H2_val_at(t));
        // project(now.H3, U.H3, problem.H3_val_at(t));
    }

    auto zero_z_up(vector_type& rhs, space const& U) const -> void {
        auto const nz = U.z.dofs();
        for (int i = 0; i < U.x.dofs(); ++i) {
            for (int j = 0; j < U.x.dofs(); ++j) {
                rhs(i, j, nz - 1) = 0;
            }
        }
    }

    auto integrate_U(double t, state& rhs) -> void {
        auto out = [](auto& buf) { return [&buf](int J, double val) { buf.data()[J] += val; }; };

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E1),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const f = problem.U1(x, t) * problem.eta;
                         return f * v.val;
                     });
        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const f = problem.U2(x, t) * problem.eta;
                         return f * v.val;
                     });
        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         auto const f = problem.U3(x, t) * problem.eta;
                         return f * v.val;
                     });
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

        // assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E1),
        //              [&](auto v, auto xx, auto const& face) {
        //                  auto const x = as_array(xx);
        //                  auto const ok = std::abs(std::get<1>(face.normal));
        //                  auto const f = mu(x) * problem.U1(x, t) + b(x) * dE1_dt(xx);
        //                  return ok * a(x) * f * v.val;
        //              });

        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
                     [&](auto v, auto xx, auto const& face) {
                         auto const x = as_array(xx);
                         // auto const ok = std::abs(std::get<2>(face.normal));
                         auto const ok = static_cast<double>(std::get<2>(face.normal) == -1);
                         auto const f = mu(x) * problem.U2(x, t) + b(x) * dE2_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        // assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
        //              [&](auto v, auto xx, auto const& face) {
        //                  auto const x = as_array(xx);
        //                  auto const ok = std::abs(std::get<0>(face.normal));
        //                  auto const f = mu(x) * problem.U3(x, t) + b(x) * dE3_dt(xx);
        //                  return ok * a(x) * f * v.val;
        //              });
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
                         // auto const ok = std::abs(std::get<2>(face.normal));
                         auto const ok = static_cast<double>(std::get<2>(face.normal) == -1);
                         auto const f = mu(x) * problem.U1(x, t) + b(x) * dE1_dt(xx);
                         return ok * a(x) * f * v.val;
                     });

        // assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
        //              [&](auto v, auto xx, auto const& face) {
        //                  auto const x = as_array(xx);
        //                  auto const ok = std::abs(std::get<0>(face.normal));
        //                  auto const f = mu(x) * problem.U2(x, t) + b(x) * dE2_dt(xx);
        //                  return ok * a(x) * f * v.val;
        //              });

        // assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
        //              [&](auto v, auto xx, auto const& face) {
        //                  auto const x = as_array(xx);
        //                  auto const ok = std::abs(std::get<1>(face.normal));
        //                  auto const f = mu(x) * problem.U3(x, t) + b(x) * dE3_dt(xx);
        //                  return ok * a(x) * f * v.val;
        //              });
    }

    auto as_array(ads::point3_t x) const -> point_type {
        return {std::get<0>(x), std::get<1>(x), std::get<2>(x)};
    }

    void after_step(int iter, double t) override {
        const auto i = iter + 1;
        const auto tt = t + steps.dt;

        if (i % 1 == 0)
            output_solution(output, i, now);

        const auto res = compute_norms(now, U, problem, tt);
        std::cout << "After step " << i << ", t = " << tt << '\n';
        print_result_info(res);

        constexpr int n = 4;
        auto print_face = [this, n, t](auto make_point) {
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= n; ++j) {
                    auto const a = static_cast<double>(i) / (n + 1);
                    auto const b = static_cast<double>(j) / (n + 1);
                    auto const s = problem.eta / problem.omega;
                    auto const p = make_point(a, b);
                    auto const U1 = problem.U1(p, t);
                    auto const U2 = problem.U2(p, t);
                    auto const U3 = problem.U3(p, t);
                    fmt::print("({:.2f}, {:.2f}, {:.2f}) -> {} x [{}, {}, {}]\n", p[0], p[1], p[2],
                               1 / s, U1 * s, U2 * s, U3 * s);
                }
            }
        };
        fmt::print(">>>>> n = -X\n");
        print_face([](double a, double b) { return point_type{0, a, b}; });
        fmt::print(">>>>> n = +X\n");
        print_face([](double a, double b) { return point_type{1, a, b}; });
        fmt::print(">>>>> n = -Y\n");
        print_face([](double a, double b) { return point_type{a, 0, b}; });
        fmt::print(">>>>> n = +Y\n");
        print_face([](double a, double b) { return point_type{a, 1, b}; });
        fmt::print(">>>>> n = -Z\n");
        print_face([](double a, double b) { return point_type{a, b, 0}; });
        fmt::print(">>>>> n = +Z\n");
        print_face([](double a, double b) { return point_type{a, b, 1}; });
    }
};

#endif  // MAXWELL_MAXWELL_CAUCHY_HPP
