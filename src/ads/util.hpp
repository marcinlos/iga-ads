#ifndef ADS_UTIL_HPP_
#define ADS_UTIL_HPP_

#include <vector>


namespace ads {

template <typename T>
T clamp(T x, T min, T max) {
    return x < min ? min : x > max ? max : x;
}

template <typename Num>
Num lerp(Num t, Num a, Num b) {
    return (1 - t) * a + t * b;
}

template <typename Int1, typename Int2, typename Num>
Num lerp(Int1 i, Int2 n, Num a, Num b) {
    Num t = static_cast<Num>(i) / n;
    return lerp(t, a, b);
}

template <typename Num>
inline std::vector<Num> linspace(Num a, Num b, std::size_t n) {
    std::vector<Num> xs(n + 1);
    for (std::size_t i = 0; i <= n; ++ i) {
        xs[i] = lerp(i, n, a, b);
    }
    return xs;
}

}

#endif /* ADS_UTIL_HPP_ */
