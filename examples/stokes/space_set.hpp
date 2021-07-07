// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef STOKES_SPACE_SET_HPP
#define STOKES_SPACE_SET_HPP

#include "ads/simulation.hpp"

namespace ads {

struct space_set {
    dimension U1x, U1y;
    dimension U2x, U2y;
    dimension Px, Py;
};

inline int total_dimension(const space_set& s) {
    auto dimU1 = s.U1x.dofs() * s.U1y.dofs();
    auto dimU2 = s.U2x.dofs() * s.U2y.dofs();
    auto dimP = s.Px.dofs() * s.Py.dofs();
    return dimU1 + dimU2 + dimP;
}

}  // namespace ads

#endif  // STOKES_SPACE_SET_HPP
