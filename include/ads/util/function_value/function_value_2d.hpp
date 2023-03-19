// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_2D_HPP
#define ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_2D_HPP

namespace ads {

struct function_value_2d {
    double val;
    double dx, dy;

    constexpr function_value_2d(double val, double dx, double dy) noexcept
    : val{val}
    , dx{dx}
    , dy{dy} { }

    constexpr function_value_2d() noexcept
    : function_value_2d{0, 0, 0} { }

    function_value_2d& operator+=(const function_value_2d& v) {
        val += v.val;
        dx += v.dx;
        dy += v.dy;
        return *this;
    }

    function_value_2d& operator-=(const function_value_2d& v) {
        val -= v.val;
        dx -= v.dx;
        dy -= v.dy;
        return *this;
    }

    function_value_2d operator-() const { return {-val, -dx, -dy}; }

    function_value_2d& operator*=(double a) {
        val *= a;
        dx *= a;
        dy *= a;
        return *this;
    }

    function_value_2d& operator/=(double a) {
        val /= a;
        dx /= a;
        dy /= a;
        return *this;
    }
};

inline function_value_2d operator+(function_value_2d x, const function_value_2d& v) {
    x += v;
    return x;
}

inline function_value_2d operator-(function_value_2d x, const function_value_2d& v) {
    x -= v;
    return x;
}

inline function_value_2d operator*(double a, function_value_2d u) {
    u *= a;
    return u;
}

inline function_value_2d operator*(function_value_2d u, double a) {
    u *= a;
    return u;
}

inline function_value_2d operator/(function_value_2d u, double a) {
    u /= a;
    return u;
}

}  // namespace ads

#endif  // ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_2D_HPP
