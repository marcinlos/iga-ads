// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_3D_RASTERIZER_HPP
#define TUMOR_3D_RASTERIZER_HPP

#include <iostream>

#include "ads/lin/tensor.hpp"
#include "ads/util/math/vec.hpp"

namespace tumor {

class rasterizer {
public:
    using point_type = ads::math::vec<3>;
    using canvas = ads::lin::tensor<double, 3>;

    void draw(point_type a, point_type b, double w, canvas& c) const {
        using std::abs;
        using std::max;

        auto sizes = c.sizes();
        int sx = sizes[0], sy = sizes[1], sz = sizes[2];

        auto d = b - a;

        int nx = static_cast<int>(abs(d.x) * sx);
        int ny = static_cast<int>(abs(d.y) * sy);
        int nz = static_cast<int>(abs(d.z) * sz);
        int nsteps = max(nx, max(ny, nz));

        point_type step = d / nsteps;

        point_type p = a;
        for (int i = 0; i < nsteps; ++i) {
            p += step;
            int ix = coord(p.x, sx);
            int iy = coord(p.y, sy);
            int iz = coord(p.z, sz);
            c(ix, iy, iz) = w;
        }
    }

private:
    int coord(double t, int s) const { return std::min(static_cast<int>(t * s), s - 1); }
};

}  // namespace tumor

#endif  // TUMOR_3D_RASTERIZER_HPP
