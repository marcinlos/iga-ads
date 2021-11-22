// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_PLANE_WAVE_PROBLEM_HPP
#define MAXWELL_PLANE_WAVE_PROBLEM_HPP

#include <cmath>

#include "problems.hpp"

class plane_wave_problem : public maxwell_problem<plane_wave_problem> {
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

    static constexpr vec3 E0 = {1, 0, 0};
    static constexpr vec3 k_hat = {0, 0, 1};

    auto eps(point_type /*x*/) const -> double { return eps0; }
    auto mu(point_type /*x*/) const -> double { return mu0; }

    auto amplitude(double s, double t) const -> double {
        return std::cos(omega * t - k * s) * excitation(t - s / c);
    }

    auto ds_amplitude(double s, double t) const -> double {
        return k * std::sin(omega * t - k * s) * excitation(t - s / c)
             - std::cos(omega * t - k * s) * dexcitation(t - s / c) / c;
    }

    auto dt_amplitude(double s, double t) const -> double {
        return -omega * std::sin(omega * t - k * s) * excitation(t - s / c)
             + std::cos(omega * t - k * s) * dexcitation(t - s / c);
    }

    auto E(point_type x, double t) const -> vec3 {
        auto const s = dot(k_hat, x);
        auto const val = amplitude(s, t);
        return E0 * val;
    }

    auto H(point_type x, double t) const -> vec3 {
        auto const s = dot(k_hat, x);
        auto const val = amplitude(s, t);
        return cross(k_hat, E0) / eta * val;
    }

    auto dE(point_type x, double t) const -> vec3 {
        auto const s = dot(k_hat, x);
        auto const val = dt_amplitude(s, t);
        return E0 * val;
    }

    auto curlE(point_type x, double t) const -> vec3 {
        auto const s = dot(k_hat, x);
        auto const val = ds_amplitude(s, t);
        return cross(k_hat, E0) * val;
    }

    auto E1(point_type x, double t) const -> value_type { return {E(x, t).x, 0, 0, 0}; }
    auto E2(point_type x, double t) const -> value_type { return {E(x, t).y, 0, 0, 0}; }
    auto E3(point_type x, double t) const -> value_type { return {E(x, t).z, 0, 0, 0}; }

    auto dE1(point_type x, double t) const -> double { return dE(x, t).x; }
    auto dE2(point_type x, double t) const -> double { return dE(x, t).y; }
    auto dE3(point_type x, double t) const -> double { return dE(x, t).z; }

    auto H1(point_type x, double t) const -> value_type { return {H(x, t).x, 0, 0, 0}; }
    auto H2(point_type x, double t) const -> value_type { return {H(x, t).y, 0, 0, 0}; }
    auto H3(point_type x, double t) const -> value_type { return {H(x, t).z, 0, 0, 0}; }

    auto excitation(double t) const -> double {
        // return t > 0 ? (1 - std::exp(-t / tau)) : 0;
        return t > 0 ? std::exp(-0.5 * std::pow((t - t0) / tau, 2)) : 0;
        // return 1.0;
    }

    auto dexcitation(double t) const -> double {
        // return t > 0 ? (1 / tau * std::exp(-t / tau)) : 0;
        return t > 0 ? -(t - t0) / (tau * tau) * std::exp(-0.5 * std::pow((t - t0) / tau, 2)) : 0;
        // return 0;
    }

    auto U(point_type x, double t) const -> vec3 {
        auto const n = normal(x);
        auto const cE = curlE(x, t);
        auto const dtE = dE(x, t);
        return cross(n, cE) / mu0 + cross(n, cross(n, dtE)) / eta;
    }

    auto U1(point_type x, double t) const -> double { return U(x, t).x; };
    auto U2(point_type x, double t) const -> double { return U(x, t).y; };
    auto U3(point_type x, double t) const -> double { return U(x, t).z; };
};

#endif  // MAXWELL_PLANE_WAVE_PROBLEM_HPP
