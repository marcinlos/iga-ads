// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MATH_VEC_FUNCTIONS_HPP
#define ADS_UTIL_MATH_VEC_FUNCTIONS_HPP

#include "ads/util/math/vec/vec_2d.hpp"
#include "ads/util/math/vec/vec_3d.hpp"
#include "ads/util/math/vec/vec_fwd.hpp"

namespace ads::math {

template <std::size_t D>
constexpr double dot(const vec<D>& u, const vec<D>& v) {
    return u.dot(v);
}

template <std::size_t D>
constexpr double norm_sq(const vec<D>& u) {
    return u.norm_sq();
}

template <std::size_t D>
constexpr double norm(const vec<D>& u) {
    return std::sqrt(norm_sq(u));
}

constexpr inline double cross(const vec<2>& u, const vec<2>& v) {
    return u.x * v.y - u.y * v.x;
}

constexpr inline vec<3> cross(const vec<3>& u, const vec<3>& v) {
    return {
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x,
    };
}

template <std::size_t D>
constexpr vec<D> normalized(vec<D> u) {
    return u / norm(u);
}

}  // namespace ads::math

#endif  // ADS_UTIL_MATH_VEC_FUNCTIONS_HPP
