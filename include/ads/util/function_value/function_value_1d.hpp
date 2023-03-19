// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_1D_HPP
#define ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_1D_HPP

namespace ads {

struct function_value_1d {
    double val;
    double dx;

    constexpr function_value_1d(double val, double dx) noexcept
    : val{val}
    , dx{dx} { }

    constexpr function_value_1d() noexcept
    : function_value_1d{0, 0} { }

    function_value_1d& operator+=(const function_value_1d& v) {
        val += v.val;
        dx += v.dx;
        return *this;
    }

    function_value_1d& operator-=(const function_value_1d& v) {
        val -= v.val;
        dx -= v.dx;
        return *this;
    }

    function_value_1d operator-() const { return {-val, -dx}; }

    function_value_1d& operator*=(double a) {
        val *= a;
        dx *= a;
        return *this;
    }

    function_value_1d& operator/=(double a) {
        val /= a;
        dx /= a;
        return *this;
    }
};

inline function_value_1d operator+(function_value_1d x, const function_value_1d& v) {
    x += v;
    return x;
}

inline function_value_1d operator-(function_value_1d x, const function_value_1d& v) {
    x -= v;
    return x;
}

inline function_value_1d operator*(double a, function_value_1d u) {
    u *= a;
    return u;
}

inline function_value_1d operator*(function_value_1d u, double a) {
    u *= a;
    return u;
}

inline function_value_1d operator/(function_value_1d u, double a) {
    u /= a;
    return u;
}

}  // namespace ads

#endif  // ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_1D_HPP
