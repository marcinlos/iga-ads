// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_SOLVER_PHASE_HPP
#define MAXWELL_SOLVER_PHASE_HPP

#include <cassert>
#include <utility>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"
#include "ads/lin/solver_ctx.hpp"
#include "ads/lin/tensor.hpp"

class solver_phase {
private:
    using matrix = ads::lin::band_matrix;
    using context = ads::lin::solver_ctx;

    struct entry {
        matrix mat;
        context ctx;

        entry()
        : ctx{mat} { }

        explicit entry(matrix mat)
        : mat{std::move(mat)}
        , ctx{this->mat} { }
    };

    using entry_table = ads::lin::tensor<entry, 2>;

    entry_table entries_;

public:
    solver_phase(int n, int m)
    : entries_{{n, m}} { }

    auto set_matrix(int i, int j, matrix mat) -> void {
        assert(i < entries_.size(0));
        assert(j < entries_.size(1));

        auto e = entry{std::move(mat)};
        ads::lin::factorize(e.mat, e.ctx);

        entries_(i, j) = std::move(e);
    }

    template <typename Rhs>
    auto operator()(Rhs& rhs) -> void {
        auto* data = rhs.data();
        auto const rhs_size = rhs.size(0);

        auto offset = 0;
        for (int j = 0; j < entries_.size(1); ++j) {
            for (int i = 0; i < entries_.size(0); ++i) {
                auto& entry = entries_(i, j);
                solve_with_factorized(entry.mat, data + offset, entry.ctx, 1);
                offset += rhs_size;
            }
        }
    }
};

#endif  // MAXWELL_SOLVER_PHASE_HPP
