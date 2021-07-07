// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MATH_VEC_OPERATORS_HPP
#define ADS_UTIL_MATH_VEC_OPERATORS_HPP

#include "ads/util/math/vec/vec_fwd.hpp"

namespace ads::math {

template <std::size_t D>
vec<D> operator+(vec<D> x, const vec<D>& v) {
    x += v;
    return x;
}

template <std::size_t D>
vec<D> operator-(vec<D> x, const vec<D>& v) {
    x -= v;
    return x;
}

template <std::size_t D>
vec<D> operator*(double a, vec<D> u) {
    u *= a;
    return u;
}

template <std::size_t D>
vec<D> operator*(vec<D> u, double a) {
    u *= a;
    return u;
}

template <std::size_t D>
vec<D> operator/(vec<D> u, double a) {
    u /= a;
    return u;
}

}  // namespace ads::math

#endif  // ADS_UTIL_MATH_VEC_OPERATORS_HPP
