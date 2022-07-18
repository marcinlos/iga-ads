// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_ANTENNA_PROBLEM_HPP
#define MAXWELL_ANTENNA_PROBLEM_HPP

#include <cmath>

#include "excitation.hpp"
#include "problems.hpp"
#include "utils.hpp"

class antenna_problem : public maxwell_problem<antenna_problem> {
public:
    static constexpr double pi = M_PI;

    static constexpr double eps0 = 8.854e-12;
    static constexpr double mu0 = 12.556e-7;
    static constexpr double c = 1 / std::sqrt(eps0 * mu0);
    static constexpr double eta = c * mu0;

    static constexpr double f = 2 * c;
    static constexpr double omega = 2 * pi * f;
    static constexpr double k = omega / c;
    static constexpr double T0 = 1 / f;
    static constexpr double tau = 2 * T0;
    static constexpr double t0 = 4 * tau;

    using vec3 = ads::math::vec<3>;

    static constexpr ads::point3_t center{1, 1, 1};

    static constexpr tapered_excitation excitation{tau};

    auto eps(point_type /*x*/) const -> double { return eps0; }
    auto mu(point_type /*x*/) const -> double { return mu0; }

    auto E(point_type p, double t) const -> vec3 {
        auto const [r, theta, phi] = to_spherical(p);
        auto const pi = M_PI;
        auto const r2 = r * r;
        auto const r3 = r * r * r;

        auto const g = g_val(t - r / c);
        auto const g0 = g_val(-r / c);
        auto const dg = dg_val(t - r / c);

        auto const s = 1 / (4 * pi * eps0);

        auto const r_coeff = (2 / r3 * (g - g0) + 2 / (c * r2) * g) * std::cos(theta);
        auto const theta_coeff =
            (1 / r3 * (g - g0) + 1 / (c * r2) * g + 1 / (c * c * r) * dg) * std::sin(theta);

        auto const v = ads::point3_t{s * r_coeff, s * theta_coeff, 0};
        return as_vec3(tangent_spherical_to_cartesian({r, theta, phi}, v));
    }

    auto H(point_type p, double t) const -> vec3 {
        auto const [r, theta, phi] = to_spherical(p);
        auto const pi = M_PI;

        auto const g = excitation.value(t - r / c);
        auto const dg = excitation.deriv(t - r / c);
        auto const s = 1 / (4 * pi);

        auto const phi_coeff = (1 / (c * r) * dg + 1 / (r * r) * g) * std::sin(theta);

        auto const v = ads::point3_t{0, 0, s * phi_coeff};
        return as_vec3(tangent_spherical_to_cartesian({r, theta, phi}, v));
    }

    auto g_val(double t) const -> double { return excitation.value(t) * std::sin(omega * t); }

    auto dg_val(double t) const -> double {
        return excitation.deriv(t) * std::sin(omega * t)
             + omega * excitation.value(t) * std::cos(omega * t);
    }

    auto to_spherical(point_type p) const -> ads::point3_t {
        return cartesian_to_spherical({p[0], p[1], p[2]}, center);
    }

    auto as_vec3(ads::point3_t v) const -> vec3 {
        auto const [x, y, z] = v;
        return {x, y, z};
    }

    auto E1(point_type x, double t) const -> value_type { return {E(x, t).x, 0, 0, 0}; }
    auto E2(point_type x, double t) const -> value_type { return {E(x, t).y, 0, 0, 0}; }
    auto E3(point_type x, double t) const -> value_type { return {E(x, t).z, 0, 0, 0}; }

    // auto dE1(point_type x, double t) const -> double { return dE(x, t).x; }
    // auto dE2(point_type x, double t) const -> double { return dE(x, t).y; }
    // auto dE3(point_type x, double t) const -> double { return dE(x, t).z; }

    auto H1(point_type x, double t) const -> value_type { return {H(x, t).x, 0, 0, 0}; }
    auto H2(point_type x, double t) const -> value_type { return {H(x, t).y, 0, 0, 0}; }
    auto H3(point_type x, double t) const -> value_type { return {H(x, t).z, 0, 0, 0}; }

    auto U1(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto U2(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto U3(point_type /*x*/, double /*t*/) const -> double { return 0; };
};

#endif  // MAXWELL_ANTENNA_PROBLEM_HPP
