// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_SCATTERING_PROBLEM_HPP
#define MAXWELL_SCATTERING_PROBLEM_HPP

#include <cmath>

#include "ads/util/math/vec.hpp"
#include "problems.hpp"

auto near(double a, double b) -> bool {
    return std::abs(a - b) < 1e-8;
}

class scattering_problem : public maxwell_problem<scattering_problem> {
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

    auto eps(point_type /*x*/) const -> double { return eps0; }
    auto mu(point_type /*x*/) const -> double { return mu0; }

    auto E1(point_type x, double t) const -> value_type {
        auto const z = x[2];
        auto const val = std::cos(omega * t - k * z) * excitation(t - z / c);
        return {val, 0, 0, 0};
    }

    auto E2(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto E3(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto H1(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto H2(point_type x, double t) const -> value_type {
        auto const z = x[2];
        auto const val = 1 / eta * std::cos(omega * t - k * z) * excitation(t - z / c);
        return {val, 0, 0, 0};
    }

    auto H3(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto excitation(double t) const -> double {
        return t > 0 ? (1 - std::exp(-t / tau)) : 0;
        // return 1.0;
    }

    auto dexcitation(double t) const -> double {
        return t > 0 ? (1 / tau * std::exp(-t / tau)) : 0;
        // return 0;
    }

    auto A(point_type x, double t) const -> double {
        auto const z = x[2];
        auto const sin = std::sin(omega * t - k * z);
        auto const cos = std::cos(omega * t - k * z);
        return omega * sin * excitation(t - z / c) - cos * dexcitation(t - z / c);
    }

    auto normal(point_type x) const -> ads::math::vec<3> {
        if (near(x[0], 0))
            return {-1, 0, 0};
        else if (near(x[0], 1))
            return {1, 0, 0};
        else if (near(x[1], 0))
            return {0, -1, 0};
        else if (near(x[1], 1))
            return {0, 1, 0};
        else if (near(x[2], 0))
            return {0, 0, -1};
        else if (near(x[2], 1))
            return {0, 0, 1};
        else
            return {0, 0, 0};
    }

    auto U(point_type x, double t) const -> ads::math::vec<3> {
        auto const xx = ads::math::vec<3>{1, 0, 0};
        auto const yy = ads::math::vec<3>{0, 1, 0};
        auto const n = normal(x);

        return (cross(n, yy) - cross(n, cross(n, xx))) * A(x, t) / eta;
    }

    auto U1(point_type x, double t) const -> double { return U(x, t).x; };
    auto U2(point_type x, double t) const -> double { return U(x, t).y; };
    auto U3(point_type x, double t) const -> double { return U(x, t).z; };
};

#endif  // MAXWELL_SCATTERING_PROBLEM_HPP
