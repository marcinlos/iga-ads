// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ERIKKSON_INVERSE_HPP
#define ERIKKSON_INVERSE_HPP

#include <array>
#include <cmath>
#include <vector>

#include "ads/bspline/bspline.hpp"

inline auto shishkin_const(int n, double eps) -> double {
    return std::log(n) / std::log(2) * eps;
}

auto choose_nodes(bool left, bool right, double width) -> std::array<double, 4>;

auto make_points(int elems, bool left, bool right, double width) -> std::vector<double>;

auto basis_from_points(std::vector<double> const& points, int p, int c) -> ads::bspline::basis;

#endif  // ERIKKSON_INVERSE_HPP
