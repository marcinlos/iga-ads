// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_SPACES_HPP
#define MAXWELL_SPACES_HPP

#include <array>

#include "ads/simulation/dimension.hpp"

struct space {
    ads::dimension x, y, z;

    auto dofs() const noexcept -> int { return x.dofs() * y.dofs() * z.dofs(); }
};

inline auto vector_shape(const space& s) noexcept -> std::array<int, 3> {
    return {s.x.dofs(), s.y.dofs(), s.z.dofs()};
}

inline auto local_shape(const space& s) noexcept -> std::array<int, 3> {
    return {
        s.x.basis.dofs_per_element(),
        s.y.basis.dofs_per_element(),
        s.z.basis.dofs_per_element(),
    };
}

inline auto local_matrix_shape(const space& s) noexcept -> std::array<int, 6> {
    auto shape = local_shape(s);
    return {
        shape[0], shape[1], shape[2],  //
        shape[0], shape[1], shape[2],  //
    };
}

inline auto factorize_matrices(space& s) -> void {
    s.x.factorize_matrix();
    s.y.factorize_matrix();
    s.z.factorize_matrix();
}

struct space_set {
    space E1, E2, E3;
    space H1, H2, H3;
};

inline auto factorize_matrices(space_set& s) -> void {
    factorize_matrices(s.E1);
    factorize_matrices(s.E2);
    factorize_matrices(s.E3);

    factorize_matrices(s.H1);
    factorize_matrices(s.H2);
    factorize_matrices(s.H3);
}

#endif  // MAXWELL_SPACES_HPP
