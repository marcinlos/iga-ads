// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_EXCITATION_HPP
#define MAXWELL_EXCITATION_HPP

#include <cmath>

struct tapered_excitation {
    double tau;

    auto value(double t) const -> double { return t > 0 ? (1 - std::exp(-t / tau)) : 0; }

    auto deriv(double t) const -> double { return t > 0 ? (1 / tau * std::exp(-t / tau)) : 0; };

    auto integral(double t) const -> double { return t > 0 ? (t + tau * std::exp(-t / tau)) : 0; };
};

struct pulse_excitation {
    double t0;
    double tau;

    auto value(double t) const -> double {
        return t > 0 ? std::exp(-0.5 * std::pow((t - t0) / tau, 2)) : 0;
    }

    auto deriv(double t) const -> double{
        return t > 0 ? -(t - t0) / (tau * tau) * std::exp(-0.5 * std::pow((t - t0) / tau, 2)) : 0;
    };

    auto integral(double t) const -> double {
        auto const s = -tau * std::sqrt(M_PI / 2);
        auto const r2 = std::sqrt(2);
        return t > 0 ? (s * (std::erf((t0 - t) / (r2 * tau)) - std::erf(t0 / (r2 * tau)))) : 0;
    };
};

#endif  // MAXWELL_EXCITATION_HPP
