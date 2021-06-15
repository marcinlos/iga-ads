#ifndef ADS_UTIL_MATH_VEC_VEC_3D_HPP_
#define ADS_UTIL_MATH_VEC_VEC_3D_HPP_


namespace ads::math {

template <>
struct vec<3> {

    using vec_type = vec<3>;

    double x, y, z;


    vec_type& operator += (const vec_type& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vec_type& operator -= (const vec_type& v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    vec_type operator - () const {
        return { -x, -y, -z };
    }

    vec_type& operator *= (double a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }

    vec_type& operator /= (double a) {
        double inv = 1 / a;
        return (*this) *= inv;
    }

    double dot(const vec_type& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    double norm_sq() const {
        return x * x + y * y + z * z;
    }
};

}

#endif /* ADS_UTIL_MATH_VEC_VEC_3D_HPP_ */
