// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cmath>
#include <tuple>

#include "ads/util/function_value.hpp"
#include "params.hpp"

constexpr double scale = 100.0;

inline auto f(double x, int n) -> double {
    auto const s = 2 * M_PI * n / scale;
    return 1 - std::cos(s * x);
}

inline auto df(double x, int n) -> double {
    auto const s = 2 * M_PI * n / scale;
    return s * std::sin(s * x);
}

inline auto ddf(double x, int n) -> double {
    auto const s = 2 * M_PI * n / scale;
    return s * s * std::cos(s * x);
}

struct solution {
    fire_params params;
    double bx;
    double by;

    double base = 300;

    // nx, ny, scale, lambda
    using entry = std::tuple<int, int, double, double>;
    std::array<entry, 3> coeffs = {{
        {1, 1, 80.0, 3.0},
        {2, 1, 30, 5.0},
        {2, 2, 110, 1.5},
    }};

    auto exact(double x, double y, double t) const -> ads::function_value_2d {
        auto a = ads::function_value_2d{base, 0.0, 0.0};

        for (auto const [nx, ny, c, lambda] : coeffs) {
            auto const fx = f(x, nx);
            auto const fy = f(y, ny);
            auto const dfx = df(x, nx);
            auto const dfy = df(y, ny);
            auto const s = std::exp(-lambda * t);

            a.val += c * fx * fy * s;
            a.dx = c * dfx * fy * s;
            a.dy = c * fx * dfy * s;
        }

        return a;
    }

    auto forcing(double x, double y, double t) const -> double {
        auto T = ads::function_value_2d{base, 0.0, 0.0};
        auto dt = 0.0;
        auto dxx = 0.0;
        auto dyy = 0.0;

        for (auto const [nx, ny, c, lambda] : coeffs) {
            auto const fx = f(x, nx);
            auto const fy = f(y, ny);
            auto const s = std::exp(-lambda * t);

            auto const dfx = df(x, nx);
            auto const dfy = df(y, ny);
            auto const ddfx = ddf(x, nx);
            auto const ddfy = ddf(y, ny);
            auto const ds = -lambda * std::exp(-lambda * t);

            T.val += c * fx * fy * s;
            T.dx = c * dfx * fy * s;
            T.dy = c * fx * dfy * s;
            dxx += c * ddfx * fy * s;
            dyy += c * fx * ddfy * s;
            dt += c * fx * fy * ds;
        }

        auto const lap = dxx + dyy;
        auto const b_grad = bx * T.dx + by * T.dy;

        auto const& p = params;

        double delta = T.val > p.Tig ? 1.0 : 0.0;
        double r = delta * p.Ar * T.val * std::exp(-p.Ta / T.val);
        double Rc = -1e3 * p.rho * p.ch * p.hc * p.M / p.M1 * r;
        double div_qr = -4 * p.sigma * p.eps * p.delta_x * std::pow(T.val, 3) * lap;
        double Qconv = p.xi * (p.T0 - T.val);
        double Qrz = p.sigma * p.eps / p.delta_z * (std::pow(p.T0, 4) - std::pow(T.val, 4));
        double rhs = Rc + Qconv + Qrz - div_qr;

        return
            // rho * cp * dT/dt
            params.rho * params.cp * dt
            //- kappa /\T
            - params.kappa * lap
            //+ b * cw * rho * \/T
            + params.rho * params.cw * b_grad
            // rhs
            - rhs;
    }
};
