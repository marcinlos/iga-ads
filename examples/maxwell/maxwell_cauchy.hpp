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
#include "plane_wave_problem.hpp"
#include "scattering_problem.hpp"
#include "spaces.hpp"
#include "state.hpp"

auto radius(double x, double y, double z) -> double {
    return std::hypot(x, y, z);
}

auto theta(double x, double y, double z) -> double {
    auto const r = std::hypot(x, y, z);
    return std::acos(z / r);
}

auto phi(double x, double y, double z) -> double {
    return std::atan2(y, x);
}

auto spherical(ads::point3_t p, ads::point3_t v) -> ads::point3_t {
    auto const [x, y, z] = p;
    auto const [vx, vy, vz] = v;

    auto const r = std::hypot(x, y, z);
    auto const theta = std::acos(z / r);
    auto const phi = std::atan2(y, x);

    auto const sphi = std::sin(phi);
    auto const cphi = std::cos(phi);
    auto const sth = std::sin(theta);
    auto const cth = std::cos(theta);

    double const A[3][3] = {
        {sth * cphi, sth * sphi, cth},
        {cth * cphi, cth * sphi, -sth},
        {-sphi, cphi, 0},
    };
    auto const sr = A[0][0] * vx + A[0][1] * vy + A[0][2] * vz;
    auto const st = A[1][0] * vx + A[1][1] * vy + A[1][2] * vz;
    auto const sp = A[2][0] * vx + A[2][1] * vy + A[2][2] * vz;

    return {sr, st, sp};
}

template <typename Vx, typename Vy, typename Vz>
auto maxwell_to_file(const std::string& path,    //
                     Vx&& Ex, Vy&& Ey, Vz&& Ez,  //
                     Vx&& Hx, Vy&& Hy, Vz&& Hz   //
                     ) -> void {
    constexpr auto res = 50;
    auto extent = fmt::format("0 {0} 0 {0} 0 {0}", res);
    auto const rx = ads::interval{0, 2};
    auto const ry = ads::interval{0, 2};
    auto const rz = ads::interval{0, 2};

    auto const for_all_points = [&](auto&& fun) {
        for (auto z : ads::evenly_spaced(rx, res)) {
            for (auto y : ads::evenly_spaced(ry, res)) {
                for (auto x : ads::evenly_spaced(rz, res)) {
                    const auto X = ads::point3_t{x, y, z};
                    fun(X);
                }
            }
        }
    };

    auto out = fmt::output_file(path);
    out.print("<?xml version=\"1.0\"?>\n");
    out.print("<VTKFile type=\"ImageData\" version=\"0.1\">\n");
    out.print("  <ImageData WholeExtent=\"{}\" origin=\"0 0 0\" spacing=\"1 1 1\">\n", extent);
    out.print("    <Piece Extent=\"{}\">\n", extent);
    out.print("      <PointData Scalars=\"r\" Vectors=\"E\">\n", extent);

    out.print("        <DataArray Name=\"E\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"3\">\n");
    for_all_points([&](auto X) { out.print("{:.7} {:.7} {:.7}\n", Ex(X), Ey(X), Ez(X)); });
    out.print("        </DataArray>\n");

    out.print("        <DataArray Name=\"H\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"3\">\n");
    for_all_points([&](auto X) { out.print("{:.7} {:.7} {:.7}\n", Hx(X), Hy(X), Hz(X)); });
    out.print("        </DataArray>\n");

    // out.print("        <DataArray Name=\"radius\" type=\"Float32\" format=\"ascii\" "
    //           "NumberOfComponents=\"1\">\n");
    // for_all_points([&](auto X) { out.print("{:.7}\n", radius(Ex(X), Ey(X), Ez(X))); });
    // out.print("        </DataArray>\n");

    // out.print("        <DataArray Name=\"theta\" type=\"Float32\" format=\"ascii\" "
    //           "NumberOfComponents=\"1\">\n");
    // for_all_points([&](auto X) { out.print("{:.7}\n", theta(Ex(X), Ey(X), Ez(X))); });
    // out.print("        </DataArray>\n");

    // out.print("        <DataArray Name=\"phi\" type=\"Float32\" format=\"ascii\" "
    //           "NumberOfComponents=\"1\">\n");
    // for_all_points([&](auto X) { out.print("{:.7}\n", phi(Ex(X), Ey(X), Ez(X))); });
    // out.print("        </DataArray>\n");

    out.print("        <DataArray Name=\"spherical\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"3\">\n");
    for_all_points([&](auto X) {
        auto const [r, theta, phi] = spherical(X, {Ex(X), Ey(X), Ez(X)});
        out.print("{:.7} {:.7} {:.7}\n", r, theta, phi);
    });
    out.print("        </DataArray>\n");

    out.print("      </PointData>\n");
    out.print("    </Piece>\n");
    out.print("  </ImageData>\n");
    out.print("</VTKFile>\n");
}

class maxwell_cauchy : public maxwell_base {
private:
    using Base = maxwell_base;
    using Problem = plane_wave_problem;

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
    , output{V.x.B, V.y.B, V.z.B, 40, 40, 40} { }

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
        // z Ex dir
        // zero(U.E2.x);
        // zero(U.E3.x);

        // x Ey dir
        // zero(U.E1.y);
        // zero(U.E3.y);

        // y Ez dir
        zero(U.E1.z);
        zero(U.E2.z);

        // U.E1.z.fix_right();

        factorize_matrices(U);
        factorize_matrices(V);

        auto tau = steps.dt;
        auto h = tau * tau / (4 * problem.eps({0.5, 0.5}) * problem.mu({0.5, 0.5}));
        auto form = [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; };

        form_matrix(Bx, V.x.basis, form);
        form_matrix(By, V.y.basis, form);
        form_matrix(Bz, V.z.basis, form);
        // form_matrix(Bz_E1, V.z.basis, form);

        // fix_dof(0, V.x, Bx);
        // fix_dof(V.x.dofs() - 1, V.x, Bx);

        // fix_dof(V.z.dofs() - 1, V.z, Bz_E1);

        ads::lin::factorize(Bx, Bx_ctx);
        ads::lin::factorize(By, By_ctx);
        ads::lin::factorize(Bz, Bz_ctx);
        // ads::lin::factorize(Bz_E1, Bz_E1_ctx);
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
        // ads_solve(rhs.E1, buffer, U.E1.x.data(), U.E1.y.data(), ads::dim_data{Bz_E1, Bz_E1_ctx});
        ads_solve(rhs.E1, buffer, U.E1.x.data(), U.E1.y.data(), ads::dim_data{Bz, Bz_ctx});
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
        apply_forcing(t, mid, V);
        substep1_boundary_E(t, mid, prev, pprev);
        // z Ex dir
        // zero_sides("x", mid.E2, U.E2);
        // zero_sides("x", mid.E3, U.E3);
        // x Ey dir
        // zero_sides("y", mid.E1, U.E1);
        // zero_sides("y", mid.E3, U.E3);
        // y Ez dir
        // zero_sides("z", mid.E1, U.E1);
        // zero_sides("z", mid.E2, U.E2);
        // zero_z_up(mid.E1, U.E1);
        substep1_solve_E(mid, buffer);

        substep1_fill_H(mid, prev, mid, U, c);
        solve_H(mid, buffer);

        // Second substep
        substep2_fill_E(now, mid, U, a, b);
        apply_forcing(t, now, V);
        substep2_boundary_E(t, now, mid, prev);
        // z Ex dir
        // zero_sides("x", now.E2, U.E2);
        // zero_sides("x", now.E3, U.E3);
        // x Ey dir
        // zero_sides("y", now.E1, U.E1);
        // zero_sides("y", now.E3, U.E3);
        // y Ez dir
        // zero_sides("z", now.E1, U.E1);
        // zero_sides("z", now.E2, U.E2);
        // zero_z_up(now.E1, U.E1);
        substep2_solve_E(now, buffer);

        substep2_fill_H(now, mid, now, U, c);
        solve_H(now, buffer);

        // zero(now.E1);
        // zero(now.E2);
        // zero(now.E3);
        // apply_forcing(t, now, V);
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
                     [&](auto v, auto xx, auto const&) {
                         auto const x = as_array(xx);
                         auto const f = problem.U1(x, t) * problem.eta;
                         return f * v.val;
                     });
        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E2),
                     [&](auto v, auto xx, auto const&) {
                         auto const x = as_array(xx);
                         auto const f = problem.U2(x, t) * problem.eta;
                         return f * v.val;
                     });
        assemble_rhs(mesh_.boundary_facets(), space_, quad_, out(rhs.E3),
                     [&](auto v, auto xx, auto const&) {
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
        // auto dE1_dt = [&](auto x) { return problem.dE1(x, t); };
        // auto dE2_dt = [&](auto x) { return problem.dE2(x, t); };
        // auto dE3_dt = [&](auto x) { return problem.dE3(x, t); };

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
        // auto dE1_dt = [&](auto x) { return problem.dE1(x, t); };
        // auto dE2_dt = [&](auto x) { return problem.dE2(x, t); };
        // auto dE3_dt = [&](auto x) { return problem.dE3(x, t); };

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

    auto apply_forcing(double t, state& rhs, space const& V) -> void {
        auto const l = 0.1 / 5;
        auto const x0 = point_type{1.0, 1.0, 1.0 - l / 2};
        auto const x1 = point_type{1.0, 1.0, 1.0 + l / 2};
        auto const len = dist(x0, x1);

        auto const q = 20;
        auto const points = quad_points(x0, x1, q);
        auto const* weigths = ads::quad::gauss::Ws[q];
        auto const scale = len / 2;

        for (int i = 0; i < q; ++i) {
            auto const& x = points[i];
            auto const W = weigths[i] * scale;
            auto const F = forcing(t, x);

            for (auto const e : elements(V.x, V.y, V.z)) {
                if (!is_inside(x, e, V))
                    continue;

                for (auto const a : dofs_on_element(e, V.x, V.y, V.z)) {
                    auto const v = eval_basis_at(x, a, V);
                    rhs.E1(a[0], a[1], a[2]) += F[0] * v.val * W;
                    rhs.E2(a[0], a[1], a[2]) += F[1] * v.val * W;
                    rhs.E3(a[0], a[1], a[2]) += F[2] * v.val * W;
                }
            }
        }
    }

    auto quad_points(point_type const x0, point_type const x1, int q) -> std::vector<point_type> {
        auto points = std::vector<point_type>(q);
        for (int i = 0; i < q; ++i) {
            const auto t = ads::quad::gauss::Xs[q][i];
            const auto s = (t + 1) / 2;
            const auto x = ads::lerp(s, x0[0], x1[0]);
            const auto y = ads::lerp(s, x0[1], x1[1]);
            const auto z = ads::lerp(s, x0[2], x1[2]);
            points[i] = point_type{x, y, z};
        }
        return points;
    }

    auto dist(point_type a, point_type b) const -> double {
        auto const dx = a[0] - b[0];
        auto const dy = a[1] - b[1];
        auto const dz = a[2] - b[2];
        return std::hypot(dx, dy, dz);
    }

    auto is_inside(point_type x, index_type e, space const& V) const -> bool {
        return V.x.B.points[e[0]] <= x[0] && x[0] <= V.x.B.points[e[0] + 1]
            && V.y.B.points[e[1]] <= x[1] && x[1] <= V.y.B.points[e[1] + 1]
            && V.z.B.points[e[2]] <= x[2] && x[2] <= V.z.B.points[e[2] + 1];
    }

    auto eval_basis_at(point_type p, index_type dof, space const& V) const -> value_type {
        const auto spanx = ads::bspline::find_span(p[0], V.x.B);
        const auto spany = ads::bspline::find_span(p[1], V.y.B);
        const auto spanz = ads::bspline::find_span(p[2], V.z.B);

        ads::bspline::eval_ders_ctx cx{x.p, 1};
        ads::bspline::eval_ders_ctx cy{y.p, 1};
        ads::bspline::eval_ders_ctx cz{z.p, 1};

        double** bvx = cx.basis_vals();
        double** bvy = cy.basis_vals();
        double** bvz = cz.basis_vals();

        eval_basis_with_derivatives(spanx, p[0], V.x.B, bvx, 1, cx);
        eval_basis_with_derivatives(spany, p[1], V.y.B, bvy, 1, cy);
        eval_basis_with_derivatives(spanz, p[2], V.z.B, bvz, 1, cz);

        int offsetx = spanx - V.x.p;
        int offsety = spany - V.y.p;
        int offsetz = spanz - V.z.p;

        int ix = dof[0] - offsetx;
        int iy = dof[1] - offsety;
        int iz = dof[2] - offsetz;

        auto value = bvx[0][ix] * bvy[0][iy] * bvz[0][iz];
        auto dx = bvx[1][ix] * bvy[0][iy] * bvz[0][iz];
        auto dy = bvx[0][ix] * bvy[1][iy] * bvz[0][iz];
        auto dz = bvx[0][ix] * bvy[0][iy] * bvz[1][iz];

        return {value, dx, dy, dz};
    }

    auto forcing(double t, point_type /*x*/) -> point_type {
        auto const s = problem.excitation(t) * std::sin(problem.omega * t);
        return {0, 0, s};
    }

    auto as_array(ads::point3_t x) const -> point_type {
        return {std::get<0>(x), std::get<1>(x), std::get<2>(x)};
    }

    void after_step(int iter, double t) override {
        const auto i = iter + 1;
        const auto tt = t + steps.dt;

        if (i % 1 == 0) {
            // output_solution(output, i, now);
            save(i);
        }

        const auto res = compute_norms(now, U, problem, tt);
        std::cout << "After step " << i << ", t = " << tt << '\n';
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
};

#endif  // MAXWELL_MAXWELL_CAUCHY_HPP
