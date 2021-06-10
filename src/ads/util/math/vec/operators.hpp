#ifndef ADS_UTIL_MATH_VEC_OPERATORS_HPP_
#define ADS_UTIL_MATH_VEC_OPERATORS_HPP_

#include "ads/util/math/vec/vec_fwd.hpp"


namespace ads::math {


template <std::size_t D>
vec<D> operator + (vec<D> x, const vec<D>& v) {
    x += v;
    return x;
}

template <std::size_t D>
vec<D> operator - (vec<D> x, const vec<D>& v) {
    x -= v;
    return x;
}

template <std::size_t D>
vec<D> operator * (double a, vec<D> u) {
    u *= a;
    return u;
}

template <std::size_t D>
vec<D> operator * (vec<D> u, double a) {
    u *= a;
    return u;
}

template <std::size_t D>
vec<D> operator / (vec<D> u, double a) {
    u /= a;
    return u;
}

}

#endif /* ADS_UTIL_MATH_VEC_OPERATORS_HPP_ */
