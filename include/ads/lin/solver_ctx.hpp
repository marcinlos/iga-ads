// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_SOLVER_CTX_HPP
#define ADS_LIN_SOLVER_CTX_HPP

#include <vector>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/dense_matrix.hpp"

namespace ads::lin {

struct solver_ctx {
    std::vector<int> pivot_vector;
    int info = 0;
    int lda;

    solver_ctx(int n, int lda)
    : pivot_vector(n)
    , lda(lda) { }

    explicit solver_ctx(const band_matrix& a)
    : solver_ctx(a.rows, 2 * a.kl + a.ku + 1) { }

    explicit solver_ctx(const dense_matrix& a)
    : solver_ctx(a.rows(), a.rows()) { }

    int* pivot() { return pivot_vector.data(); }
};

}  // namespace ads::lin

#endif  // ADS_LIN_SOLVER_CTX_HPP
