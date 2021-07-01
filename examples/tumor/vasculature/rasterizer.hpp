// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_TUMOR_VASCULATURE_RASTERIZER_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_RASTERIZER_HPP_

#include "defs.hpp"


namespace tumor::vasc {

inline int octant(vector a) {
    if (a.x > 0) {
        if (a.y > 0) {
            return a.x > a.y ? 0 : 1;
        } else {
            return a.x > -a.y ? 7 : 6;
        }
    } else {
        if (a.y > 0) {
            return -a.x > a.y ? 3 : 2;
        } else {
            return -a.x > -a.y ? 4 : 5;
        }
    }
}

inline vector to_base_octant(int octant, vector v) {
    switch (octant) {
    case 0: return {  v.x,  v.y };
    case 1: return {  v.y,  v.x };
    case 2: return {  v.y, -v.x };
    case 3: return { -v.x,  v.y };
    case 4: return { -v.x, -v.y };
    case 5: return { -v.y, -v.x };
    case 6: return { -v.y,  v.x };
    case 7: return {  v.x, -v.y };
    }
    return v;
}

inline vector from_base_octant(int octant, vector v) {
    switch (octant) {
    case 0: return {  v.x,  v.y };
    case 1: return {  v.y,  v.x };
    case 2: return { -v.y,  v.x };
    case 3: return { -v.x,  v.y };
    case 4: return { -v.x, -v.y };
    case 5: return { -v.y, -v.x };
    case 6: return {  v.y, -v.x };
    case 7: return {  v.x, -v.y };
    }
    return v;
}

template <typename Array, typename Value>
void draw_segment(vector a, vector b, Array& v, Value val) {
    int oct = octant(b - a);
    vector aa = to_base_octant(oct, a);
    vector bb = to_base_octant(oct, b);

    auto s = v.sizes();
    double sx = 1.0 / (s[0] - 1);
    double sy = 1.0 / (s[1] - 1);

    int x1 = static_cast<int>(aa.x / sx);
    int y1 = static_cast<int>(aa.y / sy);
    int x2 = static_cast<int>(bb.x / sx);
    int y2 = static_cast<int>(bb.y / sy);

    if (x1 == x2) {
        auto p = from_base_octant(oct, {x1, y1});
        int px = static_cast<int>(p.x);
        int py = static_cast<int>(p.y);
        v(px, py) = val;
        return;
    }

    double dx = x2 - x1;
    double dy = y2 - y1;

    double e = 0;
    double de = std::abs(dy / dx);

    int y = y1;
    for (int x = x1; x <= x2; ++ x) {
        auto p = from_base_octant(oct, {x, y});
        int px = static_cast<int>(p.x);
        int py = static_cast<int>(p.y);
        v(px, py) = val;
        e += de;
        while (e >= 0.5) {
            y += 1;
            e -= 1;
        }
    }
}

}

#endif /* PROBLEMS_TUMOR_VASCULATURE_RASTERIZER_HPP_ */
