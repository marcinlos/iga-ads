// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CG_SHISHKIN_HPP
#define CG_SHISHKIN_HPP

#include <cmath>

#include "ads/bspline/bspline.hpp"
#include "ads/util.hpp"

inline double shishkin_const(int n, double eps) {
    return std::log(n) / std::log(2) * eps;
}

ads::bspline::basis create_basis(double a, double b, int p, int elements, int repeated_nodes,
                                 bool adapt, double d);

#endif  // CG_SHISHKIN_HPP
