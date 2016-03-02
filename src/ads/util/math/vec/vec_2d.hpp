#ifndef ADS_UTIL_MATH_VEC_VEC_2D_HPP_
#define ADS_UTIL_MATH_VEC_VEC_2D_HPP_

#include "ads/util/math/vec/vec_fwd.hpp"


namespace ads {
namespace math {

template <>
struct vec<2> {

    using vec_type = vec<2>;

    double x, y;


    vec_type& operator += (const vec_type& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    vec_type& operator -= (const vec_type& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    vec_type operator - () const {
        return { -x, -y };
    }

    vec_type& operator *= (double a) {
        x *= a;
        y *= a;
        return *this;
    }

    vec_type& operator /= (double a) {
        double inv = 1 / a;
        return (*this) *= inv;
    }

    double dot(const vec_type& v) const {
        return x * v.x + y * v.y;
    }

    double norm_sq() const {
        return x * x + y * y;
    }
};


}
}


#endif /* ADS_UTIL_MATH_VEC_VEC_2D_HPP_ */
