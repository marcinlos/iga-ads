#ifndef ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_3D_HPP_
#define ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_3D_HPP_


namespace ads {

struct function_value_3d {
    double val;
    double dx, dy, dz;

    constexpr function_value_3d(double val, double dx, double dy, double dz) noexcept
    : val{ val }
    , dx{ dx }, dy{ dy }, dz{ dz }
    { }

    constexpr function_value_3d() noexcept
    : function_value_3d{ 0, 0, 0, 0 }
    { }

    function_value_3d& operator += (const function_value_3d& v) {
        val += v.val;
        dx  += v.dx;
        dy  += v.dy;
        dz  += v.dz;
        return *this;
    }

    function_value_3d& operator -= (const function_value_3d& v) {
        val -= v.val;
        dx  -= v.dx;
        dy  -= v.dy;
        dz  -= v.dz;
        return *this;
    }

    function_value_3d operator - () const {
        return { -val, -dx, -dy, -dz };
    }

    function_value_3d& operator *= (double a) {
        val *= a;
        dx  *= a;
        dy  *= a;
        dz  *= a;
        return *this;
    }

    function_value_3d& operator /= (double a) {
        val /= a;
        dx  /= a;
        dy  /= a;
        dz  /= a;
        return *this;
    }
};

inline function_value_3d operator + (function_value_3d x, const function_value_3d& v) {
    x += v;
    return x;
}

inline function_value_3d operator - (function_value_3d x, const function_value_3d& v) {
    x -= v;
    return x;
}

inline function_value_3d operator * (double a, function_value_3d u) {
    u *= a;
    return u;
}

inline function_value_3d operator * (function_value_3d u, double a) {
    u *= a;
    return u;
}

inline function_value_3d operator / (function_value_3d u, double a) {
    u /= a;
    return u;
}

}

#endif /* ADS_UTIL_FUNCTION_VALUE_FUNCTION_VALUE_3D_HPP_ */
