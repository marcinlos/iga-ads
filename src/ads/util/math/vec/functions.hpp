#ifndef ADS_UTIL_MATH_VEC_FUNCTIONS_HPP_
#define ADS_UTIL_MATH_VEC_FUNCTIONS_HPP_


#include "ads/util/math/vec/vec_fwd.hpp"
#include "ads/util/math/vec/vec_2d.hpp"
#include "ads/util/math/vec/vec_3d.hpp"


namespace ads {
namespace math {


template <std::size_t D>
double dot(const vec<D>& u, const vec<D>& v) {
    return u.dot(v);
}

template <std::size_t D>
double norm_sq(const vec<D>& u) {
    return u.norm_sq();
}


template <std::size_t D>
double norm(const vec<D>& u) {
    return std::sqrt(norm_sq(u));
}


inline double cross(const vec<2>& u, const vec<2>& v) {
    return u.x * v.y - u.y * v.x;
}


inline vec<3> cross(const vec<3>& u, const vec<3>& v) {
    return {
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x
    };
}

template <std::size_t D>
vec<D> normalized(vec<D> u) {
    return u /= norm(u);
}



}
}


#endif /* ADS_UTIL_MATH_VEC_FUNCTIONS_HPP_ */
