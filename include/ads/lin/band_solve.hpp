// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_BAND_SOLVE_HPP
#define ADS_LIN_BAND_SOLVE_HPP

#include <iostream>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/lapack.hpp"
#include "ads/lin/solver_ctx.hpp"
#include "ads/lin/tensor.hpp"

namespace ads::lin {

inline void factorize(band_matrix& a, solver_ctx& ctx) {
    dgbtrf_(&a.rows, &a.cols, &a.kl, &a.ku, a.full_buffer(), &ctx.lda, ctx.pivot(), &ctx.info);
}

template <typename Rhs>
inline void solve_with_factorized(const band_matrix& a, Rhs& b, solver_ctx& ctx) {
    int nrhs = b.size() / b.size(0);
    const char* trans = "No transpose";

    dgbtrs_(trans, &a.cols, &a.kl, &a.ku, &nrhs, a.full_buffer(), &ctx.lda, ctx.pivot(), b.data(),
            &a.cols, &ctx.info);
}

template <typename Rhs>
inline void solve(band_matrix& a, Rhs& b, solver_ctx& ctx) {
    factorize(a, ctx);
    solve_with_factorized(a, b, ctx);
}

}  // namespace ads::lin

#endif  // ADS_LIN_BAND_SOLVE_HPP
