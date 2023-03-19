// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_HPP
#define ADS_UTIL_HPP

#include <vector>

namespace ads {

template <typename T, typename TMin, typename TMax>
auto clamp(T v, TMin min, TMax max) {
    return v < min ? min : v > max ? max : v;
}

template <typename Num, typename Val>
Val lerp(Num t, Val a, Val b) {
    return (1 - t) * a + t * b;
}

template <typename Int1, typename Int2, typename Num>
Num lerp(Int1 i, Int2 n, Num a, Num b) {
    Num t = static_cast<Num>(i) / static_cast<Num>(n);
    return lerp(t, a, b);
}

template <typename Num>
inline std::vector<Num> linspace(Num a, Num b, std::size_t n) {
    std::vector<Num> xs(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        xs[i] = lerp(i, n, a, b);
    }
    return xs;
}

template <typename T>
auto as_signed(T a) {
    return std::make_signed_t<T>(a);
}

template <typename T>
auto as_unsigned(T a) {
    return std::make_unsigned_t<T>(a);
}

template <typename T, typename U>
constexpr T narrow_cast(U&& u) noexcept {
    return static_cast<T>(std::forward<U>(u));
}

}  // namespace ads

#endif  // ADS_UTIL_HPP
