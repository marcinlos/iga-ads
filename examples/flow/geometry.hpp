// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef FLOW_GEOMETRY_HPP
#define FLOW_GEOMETRY_HPP

#include <cmath>
#include <ostream>

#include "ads/util.hpp"


namespace ads {

struct vec3d {
    double x;
    double y;
    double z;
};

inline std::ostream& operator << (std::ostream& os, const vec3d& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

inline vec3d operator - (const vec3d& a, const vec3d& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline vec3d operator + (const vec3d& b, const vec3d& a) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline vec3d operator * (double c, const vec3d& a) {
    return {c * a.x, c * a.y, c * a.z};
}

inline double dot(const vec3d& a, const vec3d& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double len_sq(const vec3d& a) {
    return dot(a, a);
}

inline double len(const vec3d& a) {
    return std::sqrt(len_sq(a));
}


inline double falloff(double r, double R, double t) {
    if (t < r) return 1.0;
    if (t > R) return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

// r < R in [0, 1]
inline double bump(double r, double R, double x, double y, double z) {
    double dx = x - 0.5;
    double dy = y - 0.5;
    double dz = z - 0.5;
    double t = std::sqrt(dx * dx + dy * dy + dz * dz);
    return falloff(r / 2, R / 2, t);
}


inline double dist_from_segment(const vec3d& p, const vec3d& a, const vec3d& b) {
    auto s = b - a;
    double proj = dot(p - a, s) / len_sq(s);
    auto closest = lerp(clamp(proj, 0, 1), a, b);
    return len(p - closest);
}



}


#endif // FLOW_GEOMETRY_HPP
