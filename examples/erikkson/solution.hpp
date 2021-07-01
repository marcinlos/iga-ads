// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_ERIKKSON_SOLUTION_HPP_
#define PROBLEMS_ERIKKSON_SOLUTION_HPP_

#include <cmath>

#include "ads/util/function_value.hpp"


namespace ads {

    inline function_value_2d erikkson_exact(double x, double y, double eps) {
        using std::exp;
        constexpr double pi = M_PI;
        auto lambda = pi * eps;
        auto del = std::sqrt(1 + 4 * lambda * lambda);
        auto r1 = (1 + del) / (2 * eps);
        auto r2 = (1 - del) / (2 * eps);

        auto norm = exp(-r1) - exp(-r2);
        auto alpha = (exp(r1 * (x - 1)) - exp(r2 * (x - 1))) / norm;
        double val = alpha * std::sin(pi * y);

        return {
            val,
            (r1 * exp(r1 * (x - 1)) - r2 * exp(r2 * (x - 1))) / norm * std::sin(pi * y),
            alpha * pi * std::cos(pi * y)
        };
    }

    inline double erikkson_forcing(double x, double y, double eps, double t) {
        constexpr double pi = M_PI;
        auto s = [](double a) { return std::sin(pi * a); };
        auto c = [](double a) { return std::cos(pi * a); };
        return pi*s(x)*s(y)*c(t) + 2*pi*pi*eps*s(x)*s(y)*s(t) + pi*c(x)*s(y)*s(t);
    }

    inline function_value_2d erikkson_nonstationary_exact(double x, double y, double t) {
        constexpr double pi = M_PI;
        auto s = [](double a) { return std::sin(pi * a); };
        auto c = [](double a) { return std::cos(pi * a); };

        return function_value_2d{ s(t)*s(x)*s(y), pi*s(t)*c(x)*s(y), pi*s(t)*s(x)*c(y) };
    }

    inline function_value_2d erikkson2_exact(double x, double y, double eps) {
        using std::exp;
        const double Pe = 1 / eps;

        auto g  = [Pe](double t) { return t + (exp(Pe * t) - 1) / (1 - exp(Pe)); };
        auto dg = [Pe](double t) { return 1 + Pe * exp(Pe * t) / (1 - exp(Pe)); };

        return { g(x) * g(y), dg(x) * g(y), g(x) * dg(y) };
    }

    inline double erikkson2_forcing(double x, double y, double eps) {
        using std::exp;
        const double Pe = 1 / eps;
        const double c = 1 / (1 - exp(Pe));

        auto g   = [Pe,c](double t) { return t + c * (exp(Pe * t) - 1); };
        auto dg  = [Pe,c](double t) { return 1 + c *Pe * exp(Pe * t); };
        auto ddg = [Pe,c](double t) { return c * Pe * Pe * exp(Pe * t); };

        return -eps * (ddg(x) * g(y) + g(x) * ddg(y)) + dg(x) * g(y) + g(x) * dg(y);
    }



}

#endif // PROBLEMS_ERIKKSON_SOLUTION_HPP_
