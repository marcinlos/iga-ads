// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "maxwell_base.hpp"

#include <algorithm>

void validate_args(lyra::cli const& cli, lyra::parse_result const& result, bool show_help) {
    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }
    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }
}

auto fix_dof(int k, ads::dimension const& dim, ads::lin::band_matrix& K) -> void {
    auto const last = dim.dofs() - 1;

    auto const low = std::clamp(k - dim.p, 0, last);
    auto const high = std::clamp(k + dim.p, 0, last);

    for (int i = low; i <= high; ++i) {
        K(k, i) = 0;
    }
    K(k, k) = 1;
}

auto zero_sides(std::string_view dims, ads::lin::tensor<double, 3>& rhs, space const& U) -> void {
    auto const nx = U.x.dofs();
    auto const ny = U.y.dofs();
    auto const nz = U.z.dofs();

    for (const char dim : dims) {
        if (dim == 'x') {
            for (int i = 0; i < ny; ++i) {
                for (int j = 0; j < nz; ++j) {
                    rhs(0, i, j) = 0;
                    rhs(nx - 1, i, j) = 0;
                }
            }
        } else if (dim == 'y') {
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < nz; ++j) {
                    rhs(i, 0, j) = 0;
                    rhs(i, ny - 1, j) = 0;
                }
            }
        } else if (dim == 'z') {
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    rhs(i, j, 0) = 0;
                    rhs(i, j, nz - 1) = 0;
                }
            }
        }
    }
}

auto set_boundary_conditions(space_set& s) -> void {
    auto const zero = [](auto& dim) {
        dim.fix_left();
        dim.fix_right();
    };

    // E x n = 0
    zero(s.E1.y);
    zero(s.E1.z);
    zero(s.E2.x);
    zero(s.E2.z);
    zero(s.E3.x);
    zero(s.E3.y);

    // H * n = 0
    zero(s.H1.x);
    zero(s.H2.y);
    zero(s.H3.z);
}

auto dof_support(int dof, ads::dimension const& V) -> interval {
    auto const [first, last] = V.basis.element_ranges[dof];
    auto const a = V.B.points[first];
    auto const b = V.B.points[last + 1];
    return {a, b};
}
