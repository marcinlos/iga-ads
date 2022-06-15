// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_HEAD_ANTENNA_PROBLEM_HPP
#define MAXWELL_HEAD_ANTENNA_PROBLEM_HPP

#include <cmath>
#include <string_view>

#include "ads/experimental/all.hpp"
#include "head_data.hpp"
#include "problems.hpp"

class head_antenna_problem : public maxwell_problem<head_antenna_problem> {
private:
    head_data data_;

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

    using point_type = head_data::point_type;

    explicit head_antenna_problem(std::string_view data_path)
    : data_{read_density_data(data_path)} { }

    auto eps(point_type x) const -> double { return eps0 * data_.eps(map(x)); }

    auto mu(point_type x) const -> double { return mu0 * data_.mu(map(x)); }

    auto sigma(point_type x) const -> double { return eps0 * omega * data_.sigma(map(x)); }

    auto empty(point_type x) const -> bool { return data_.empty(x); }

    auto map(point_type x) const -> point_type {
        auto const rx = ads::interval{0.4, 0.8};
        auto const ry = ads::interval{0.8, 1.2};
        auto const rz = ads::interval{0.8, 1.2};

        if (inside(x[0], rx) && inside(x[1], ry) && inside(x[2], rz)) {
            auto const tx = (x[0] - rx.left) / length(rx);
            auto const ty = (x[1] - ry.left) / length(ry);
            auto const tz = (x[2] - rz.left) / length(rz);
            return {tx, ty, tz};
        } else {
            return {0, 0, 0};
        }
    }

    auto inside(double x, ads::interval s) const -> bool { return x >= s.left && x <= s.right; }

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

    auto J1(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto J2(point_type /*x*/, double /*t*/) const -> double { return 0; };
    auto J3(point_type /*x*/, double /*t*/) const -> double { return 0; };
};

#endif  // MAXWELL_HEAD_ANTENNA_PROBLEM_HPP
