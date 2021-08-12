// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MAXWELL_BASE_HPP
#define MAXWELL_MAXWELL_BASE_HPP

#include <array>
#include <cmath>
#include <string_view>

#include "ads/executor/galois.hpp"
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/output_manager.hpp"
#include "ads/projection.hpp"
#include "ads/simulation.hpp"
#include "ads/simulation/dimension.hpp"
#include "results.hpp"
#include "spaces.hpp"
#include "state.hpp"

auto fix_dof(int k, ads::dimension const& dim, ads::lin::band_matrix& K) -> void;

auto zero_sides(std::string_view dims, ads::lin::tensor<double, 3>& rhs, space const& U) -> void;

auto set_boundary_conditions(space_set& s) -> void;

struct interval {
    double a, b;
};

auto dof_support(int dof, ads::dimension const& V) -> interval;

template <typename Fun>
auto project(ads::lin::tensor<double, 3>& rhs, space& V, Fun&& fun) -> void {
    auto const f = [&](double x, double y, double z) { return fun({x, y, z}); };
    auto const shape = vector_shape(V);
    auto buffer = ads::lin::tensor<double, 3>{shape};
    compute_projection(rhs, V.x.basis, V.y.basis, V.z.basis, f);
    ads_solve(rhs, buffer, V.x.data(), V.y.data(), V.z.data());
}

template <typename Problem>
auto set_init_state(state& s, space_set& U, Problem const& problem) -> void {
    project(s.E1, U.E1, problem.init_E1());
    project(s.E2, U.E2, problem.init_E2());
    project(s.E3, U.E3, problem.init_E3());

    project(s.H1, U.H1, problem.init_H1());
    project(s.H2, U.H2, problem.init_H2());
    project(s.H3, U.H3, problem.init_H3());
}

class maxwell_base : public ads::simulation_3d {
private:
    using Base = ads::simulation_3d;

    ads::galois_executor executor{4};

protected:
    explicit maxwell_base(ads::config_3d const& config)
    : Base{config} { }

    auto dof_support(index_type dof, space const& V) const -> std::array<interval, 3> {
        auto const [ix, iy, iz] = dof;
        using ::dof_support;
        return {
            dof_support(ix, V.x),
            dof_support(iy, V.y),
            dof_support(iz, V.z),
        };
    }

    template <typename A, typename B>
    auto substep1_fill_E(state& rhs, state& prev, space_set const& U, A&& a, B&& b) -> void {
        constexpr int X = 0;
        constexpr int Y = 1;
        constexpr int Z = 2;

        compute_rhs(rhs.E1, prev, prev, U, U.E1, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[X].val + a(x) * (H[Z].dy - H[Y].dz)) * v.val + b(x) * E[Y].dx * v.dy;
        });

        compute_rhs(rhs.E2, prev, prev, U, U.E2, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Y].val + a(x) * (H[X].dz - H[Z].dx)) * v.val + b(x) * E[Z].dy * v.dz;
        });

        compute_rhs(rhs.E3, prev, prev, U, U.E3, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Z].val + a(x) * (H[Y].dx - H[X].dy)) * v.val + b(x) * E[X].dz * v.dx;
        });

        zero_sides("yz", rhs.E1, U.E1);
        zero_sides("xz", rhs.E2, U.E2);
        zero_sides("xy", rhs.E3, U.E3);
    }

    template <typename C>
    auto substep1_fill_H(state& rhs, state& prev, state& mid, space_set const& U, C&& c) -> void {
        constexpr int X = 0;
        constexpr int Y = 1;
        constexpr int Z = 2;

        compute_rhs(rhs.H1, prev, mid, U, U.H1, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[X].val - c(x) * (E[Z].dy - En[Y].dz)) * v.val;
        });

        compute_rhs(rhs.H2, prev, mid, U, U.H2, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Y].val - c(x) * (E[X].dz - En[Z].dx)) * v.val;
        });

        compute_rhs(rhs.H3, prev, mid, U, U.H3, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Z].val - c(x) * (E[Y].dx - En[X].dy)) * v.val;
        });

        zero_sides("x", rhs.H1, U.H1);
        zero_sides("y", rhs.H2, U.H2);
        zero_sides("z", rhs.H3, U.H3);
    }

    template <typename A, typename B>
    auto substep2_fill_E(state& rhs, state& prev, space_set const& U, A&& a, B&& b) -> void {
        constexpr int X = 0;
        constexpr int Y = 1;
        constexpr int Z = 2;

        compute_rhs(rhs.E1, prev, prev, U, U.E1, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[X].val + a(x) * (H[Z].dy - H[Y].dz)) * v.val + b(x) * E[Z].dx * v.dz;
        });

        compute_rhs(rhs.E2, prev, prev, U, U.E2, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Y].val + a(x) * (H[X].dz - H[Z].dx)) * v.val + b(x) * E[X].dy * v.dx;
        });

        compute_rhs(rhs.E3, prev, prev, U, U.E3, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Z].val + a(x) * (H[Y].dx - H[X].dy)) * v.val + b(x) * E[Y].dz * v.dy;
        });

        zero_sides("yz", rhs.E1, U.E1);
        zero_sides("xz", rhs.E2, U.E2);
        zero_sides("xy", rhs.E3, U.E3);
    }

    template <typename C>
    auto substep2_fill_H(state& rhs, state& prev, state& mid, space_set const& U, C&& c) -> void {
        constexpr int X = 0;
        constexpr int Y = 1;
        constexpr int Z = 2;

        compute_rhs(rhs.H1, prev, mid, U, U.H1, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[X].val - c(x) * (En[Z].dy - E[Y].dz)) * v.val;
        });

        compute_rhs(rhs.H2, prev, mid, U, U.H2, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Y].val - c(x) * (En[X].dz - E[Z].dx)) * v.val;
        });

        compute_rhs(rhs.H3, prev, mid, U, U.H3, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Z].val - c(x) * (En[Y].dx - E[X].dy)) * v.val;
        });

        zero_sides("x", rhs.H1, U.H1);
        zero_sides("y", rhs.H2, U.H2);
        zero_sides("z", rhs.H3, U.H3);
    }

    template <typename RHS, typename Form>
    void compute_rhs(RHS& rhs, state const& prev, state const& mid, space_set const& U,
                     space const& V, Form&& form) {
        zero(rhs);
        auto const shape = ::local_shape(V);

        executor.for_each(elements(V.x, V.y, V.z), [&](auto const e) {
            auto loc = vector_type{shape};

            auto const J = jacobian(e, V.x, V.y, V.z);
            for (auto const q : quad_points(V.x, V.y, V.z)) {
                auto const W = weight(q, V.x, V.y, V.z);
                auto const x = point(e, q, V.x, V.y, V.z);

                using vec = std::array<value_type, 3>;

                auto const E1 = eval(prev.E1, e, q, U.E1.x, U.E1.y, U.E1.z);
                auto const E2 = eval(prev.E2, e, q, U.E2.x, U.E2.y, U.E2.z);
                auto const E3 = eval(prev.E3, e, q, U.E3.x, U.E3.y, U.E3.z);
                auto const E = vec{E1, E2, E3};

                auto const E1_n = eval(mid.E1, e, q, U.E1.x, U.E1.y, U.E1.z);
                auto const E2_n = eval(mid.E2, e, q, U.E2.x, U.E2.y, U.E2.z);
                auto const E3_n = eval(mid.E3, e, q, U.E3.x, U.E3.y, U.E3.z);
                auto const E_n = vec{E1_n, E2_n, E3_n};

                auto const H1 = eval(prev.H1, e, q, U.H1.x, U.H1.y, U.H1.z);
                auto const H2 = eval(prev.H2, e, q, U.H2.x, U.H2.y, U.H2.z);
                auto const H3 = eval(prev.H3, e, q, U.H3.x, U.H3.y, U.H3.z);
                auto const H = vec{H1, H2, H3};

                for (auto const a : dofs_on_element(e, V.x, V.y, V.z)) {
                    auto const aa = dof_global_to_local(e, a, V.x, V.y, V.z);
                    auto const v = eval_basis(e, q, a, V.x, V.y, V.z);

                    auto const val = form(E, E_n, H, v, x);
                    loc(aa[0], aa[1], aa[2]) += val * W * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, loc, e, V.x, V.y, V.z); });
        });
    }

    template <typename Problem>
    auto compute_norms(state const& s, space_set const& U, Problem const& problem, double t) const
        -> maxwell_result_info {
        return {
            compute_norms_E(s, U, problem, t),
            compute_norms_H(s, U, problem, t),
        };
    }

    template <typename Problem>
    auto compute_norms_E(state const& s, space_set const& U, Problem const& problem, double t) const
        -> vector_result_info {
        auto const norm_x = scalar_norm{
            normL2(s.E1, U.E1.x, U.E1.y, U.E1.z),
            normH1(s.E1, U.E1.x, U.E1.y, U.E1.z),
        };
        auto const norm_y = scalar_norm{
            normL2(s.E2, U.E2.x, U.E2.y, U.E2.z),
            normH1(s.E2, U.E2.x, U.E2.y, U.E2.z),
        };
        auto const norm_z = scalar_norm{
            normL2(s.E3, U.E3.x, U.E3.y, U.E3.z),
            normH1(s.E3, U.E3.x, U.E3.y, U.E3.z),
        };
        auto const rot = norm_rot(s.E1, s.E2, s.E3,        //
                                  U.E1.x, U.E1.y, U.E1.z,  //
                                  U.E2.x, U.E2.y, U.E2.z,  //
                                  U.E3.x, U.E3.y, U.E3.z);

        auto const div = norm_div(s.E1, s.E2, s.E3,        //
                                  U.E1.x, U.E1.y, U.E1.z,  //
                                  U.E2.x, U.E2.y, U.E2.z,  //
                                  U.E3.x, U.E3.y, U.E3.z);

        auto const norm = vector_norm{norm_x, norm_y, norm_z, rot, div};

        auto const err_x = scalar_norm{
            errorL2(s.E1, U.E1.x, U.E1.y, U.E1.z, problem.E1_at(t)),
            errorH1(s.E1, U.E1.x, U.E1.y, U.E1.z, problem.E1_at(t)),
        };
        auto const err_y = scalar_norm{
            errorL2(s.E2, U.E2.x, U.E2.y, U.E2.z, problem.E2_at(t)),
            errorH1(s.E2, U.E2.x, U.E2.y, U.E2.z, problem.E2_at(t)),
        };
        auto const err_z = scalar_norm{
            errorL2(s.E3, U.E3.x, U.E3.y, U.E3.z, problem.E3_at(t)),
            errorH1(s.E3, U.E3.x, U.E3.y, U.E3.z, problem.E3_at(t)),
        };
        auto const err_rot = error_rot(s.E1, s.E2, s.E3,        //
                                       U.E1.x, U.E1.y, U.E1.z,  //
                                       U.E2.x, U.E2.y, U.E2.z,  //
                                       U.E3.x, U.E3.y, U.E3.z,  //
                                       problem.E1_at(t), problem.E2_at(t), problem.E3_at(t));
        auto const err_div = div - 0.0;  // exact fields are divergence free

        auto const err = vector_norm{err_x, err_y, err_z, err_rot, err_div};

        return {norm, err};
    }

    template <typename Problem>
    auto compute_norms_H(state const& s, space_set const& U, Problem const& problem, double t) const
        -> vector_result_info {
        auto const norm_x = scalar_norm{
            normL2(s.H1, U.H1.x, U.H1.y, U.H1.z),
            normH1(s.H1, U.H1.x, U.H1.y, U.H1.z),
        };
        auto const norm_y = scalar_norm{
            normL2(s.H2, U.H2.x, U.H2.y, U.H2.z),
            normH1(s.H2, U.H2.x, U.H2.y, U.H2.z),
        };
        auto const norm_z = scalar_norm{
            normL2(s.H3, U.H3.x, U.H3.y, U.H3.z),
            normH1(s.H3, U.H3.x, U.H3.y, U.H3.z),
        };
        auto const rot = norm_rot(s.H1, s.H2, s.H3,        //
                                  U.H1.x, U.H1.y, U.H1.z,  //
                                  U.H2.x, U.H2.y, U.H2.z,  //
                                  U.H3.x, U.H3.y, U.H3.z);

        auto const div = norm_div(s.H1, s.H2, s.H3,        //
                                  U.H1.x, U.H1.y, U.H1.z,  //
                                  U.H2.x, U.H2.y, U.H2.z,  //
                                  U.H3.x, U.H3.y, U.H3.z);

        auto const norm = vector_norm{norm_x, norm_y, norm_z, rot, div};

        auto const err_x = scalar_norm{
            errorL2(s.H1, U.H1.x, U.H1.y, U.H1.z, problem.H1_at(t)),
            errorH1(s.H1, U.H1.x, U.H1.y, U.H1.z, problem.H1_at(t)),
        };
        auto const err_y = scalar_norm{
            errorL2(s.H2, U.H2.x, U.H2.y, U.H2.z, problem.H2_at(t)),
            errorH1(s.H2, U.H2.x, U.H2.y, U.H2.z, problem.H2_at(t)),
        };
        auto const err_z = scalar_norm{
            errorL2(s.H3, U.H3.x, U.H3.y, U.H3.z, problem.H3_at(t)),
            errorH1(s.H3, U.H3.x, U.H3.y, U.H3.z, problem.H3_at(t)),
        };
        auto const err_rot = error_rot(s.H1, s.H2, s.H3,        //
                                       U.H1.x, U.H1.y, U.H1.z,  //
                                       U.H2.x, U.H2.y, U.H2.z,  //
                                       U.H3.x, U.H3.y, U.H3.z,  //
                                       problem.H1_at(t), problem.H2_at(t), problem.H3_at(t));

        auto const err_div = div - 0.0;  // exact fields are divergence free

        auto const err = vector_norm{err_x, err_y, err_z, err_rot, err_div};

        return {norm, err};
    }

    auto output_solution(ads::output_manager<3>& out, int iter, state const& s) -> void {
        out.to_file("out_%d.vti", iter,  //
                    out.evaluate(s.E1),  //
                    out.evaluate(s.E2),  //
                    out.evaluate(s.E3),  //
                    out.evaluate(s.H1),  //
                    out.evaluate(s.H2),  //
                    out.evaluate(s.H3));
    }
};

#endif  // MAXWELL_MAXWELL_BASE_HPP
