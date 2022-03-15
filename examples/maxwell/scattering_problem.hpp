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

    static constexpr double c = 3 * 1e8;
    static constexpr double eps0 = 8.854e-12;
    static constexpr double mu0 = 12.556e-7;
    static constexpr double Z = c * mu0;  // 120 * pi;
    static constexpr double eta = Z;      // std::sqrt(mu0 * eps0);

    static constexpr double f = 3 * 1e8;
    static constexpr double omega = 2 * pi * f;
    static constexpr double k = omega / c;

    auto eps(point_type /*x*/) const -> double { return eps0; }
    auto mu(point_type /*x*/) const -> double { return mu0; }

    auto E1(point_type x, double t) const -> value_type {
        auto const val = std::cos(omega * t - k * x[2]) * excitation(t);
        return {val, 0, 0, 0};
    }

    auto E2(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto E3(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto H1(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto H2(point_type x, double t) const -> value_type {
        auto const val = 1 / Z * std::cos(omega * t - k * x[2]) * excitation(t);
        return {val, 0, 0, 0};
    }

    auto H3(point_type /*x*/, double /*t*/) const -> value_type { return {}; }

    auto Uval(point_type x, double t) const -> double {
        return omega / eta * std::sin(omega * t - k * x[2]);
    }

    auto excitation(double t) const -> double {
        constexpr double tau = 1 / omega;
        return (1 - std::exp(-t / tau)) * std::sin(omega * t);
        // return 1.0;
    }

    auto dexcitation(double t) const -> double {
        constexpr double tau = 1 / omega;
        return 1 / tau * std::exp(-t / tau) * std::sin(omega * t)
             + (1 - std::exp(-t / tau)) * omega * std::cos(omega * t);
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
        auto const e = excitation(t);
        auto const de = dexcitation(t);
        auto const sin = std::sin(omega * t - k * x[2]);
        auto const cos = std::cos(omega * t - k * x[2]);

        auto const curlE = k * sin;
        auto const dE_dt = -omega / eta * sin * e + cos * de;

        return cross(n, yy) * curlE / mu0 + cross(n, cross(n, xx)) * dE_dt;
    }

    auto U1(point_type x, double t) const -> double { return U(x, t).x; };
    auto U2(point_type x, double t) const -> double { return U(x, t).y; };
    auto U3(point_type x, double t) const -> double { return U(x, t).z; };
};

#endif  // MAXWELL_SCATTERING_PROBLEM_HPP
