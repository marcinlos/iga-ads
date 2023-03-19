// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MATH_VEC_VEC_2D_HPP
#define ADS_UTIL_MATH_VEC_VEC_2D_HPP

#include "ads/util/math/vec/vec_fwd.hpp"

namespace ads::math {

template <>
struct vec<2> {
    using vec_type = vec<2>;

    double x, y;

    constexpr vec_type& operator+=(const vec_type& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    constexpr vec_type& operator-=(const vec_type& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    constexpr vec_type operator-() const { return {-x, -y}; }

    constexpr vec_type& operator*=(double a) {
        x *= a;
        y *= a;
        return *this;
    }

    constexpr vec_type& operator/=(double a) {
        double inv = 1 / a;
        return (*this) *= inv;
    }

    constexpr double dot(const vec_type& v) const { return x * v.x + y * v.y; }

    constexpr double norm_sq() const { return x * x + y * y; }
};

}  // namespace ads::math

#endif  // ADS_UTIL_MATH_VEC_VEC_2D_HPP
