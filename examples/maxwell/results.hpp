// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_RESULTS_HPP
#define MAXWELL_RESULTS_HPP

#include <cmath>

struct scalar_norm {
    double L2;
    double H1;
};

struct vector_norm {
    scalar_norm x, y, z;
    double rot;
    double div;

    auto total() const -> scalar_norm {
        return {
            std::sqrt(x.L2 * x.L2 + y.L2 * y.L2 + z.L2 * z.L2),
            std::sqrt(x.H1 * x.H1 + y.H1 * y.H1 + z.H1 * z.H1),
        };
    }
};

struct vector_result_info {
    vector_norm norm;
    vector_norm error;

    auto relative_error() const -> scalar_norm {
        auto const total_norm = norm.total();
        auto const total_error = error.total();
        return {
            total_error.L2 / total_norm.L2,
            total_error.H1 / total_norm.H1,
        };
    }
};

struct maxwell_result_info {
    vector_result_info E;
    vector_result_info H;
};

auto print_result_info(maxwell_result_info const& res) -> void;

#endif  // MAXWELL_RESULTS_HPP
