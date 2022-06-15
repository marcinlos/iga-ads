// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MT_PROBLEM_HPP
#define MAXWELL_MT_PROBLEM_HPP

#include <cmath>

#include "ads/experimental/all.hpp"
#include "problems.hpp"

class mt_problem : public maxwell_problem<mt_problem> {
public:
    static constexpr double pi = M_PI;

    static constexpr double eps0 = 8.854e-12;
    static constexpr double mu0 = 12.556e-7;
    static constexpr double c = 1 / std::sqrt(eps0 * mu0);
    static constexpr double eta = c * mu0;
    static constexpr double sigma1 = 1e-1;
    static constexpr double sigma_air = 1e-10;

    auto eps(point_type /*x*/) const -> double { return eps0; }

    auto mu(point_type /*x*/) const -> double { return mu0; }

    auto sigma(point_type x) const -> double {
        auto const z = x[2];
        return z < 60e3 ? sigma1 : sigma_air;
    }

    auto empty(point_type /*x*/) const -> bool { return true; }

    constexpr static value_type unknown{NAN, NAN, NAN, NAN};

    auto E1(point_type, double) const -> value_type { return unknown; }
    auto E2(point_type, double) const -> value_type { return unknown; }
    auto E3(point_type, double) const -> value_type { return unknown; }

    auto H1(point_type, double) const -> value_type { return unknown; }
    auto H2(point_type, double) const -> value_type { return unknown; }
    auto H3(point_type, double) const -> value_type { return unknown; }

    auto zero() const {
        return [](point_type /*x*/) { return 0.0; };
    }

    auto init_E1() const { return zero(); }
    auto init_E2() const { return zero(); }
    auto init_E3() const { return zero(); }

    auto init_H1() const { return zero(); }
    auto init_H2() const { return zero(); }
    auto init_H3() const { return zero(); }

    auto U1(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto U2(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto U3(point_type /*x*/, double /*t*/) const -> double { return 0; };

    auto in_source(point_type x) const -> bool {
        auto const [xx, yy, zz] = x;
        return 150e3 < xx && xx < 300e3  //
            && 150e3 < yy && yy < 300e3  //
            && 100e3 < zz && zz < 110e3;
    }

    auto J1(point_type x, double t) const -> double{
        return (t == 0 && in_source(x)) ? 1.0 : 0.0;
    };

    auto J2(point_type x, double t) const -> double {
        return (t == 0 && in_source(x)) ? 1.0 : 0.0;
    };

    auto J3(point_type /*x*/, double /*t*/) const -> double { return 0; };
};

#endif  // MAXWELL_MT_PROBLEM_HPP
