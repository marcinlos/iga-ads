#ifndef ADS_LIN_DENSE_SOLVE_HPP_
#define ADS_LIN_DENSE_SOLVE_HPP_

#include <iostream>

#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/lapack.hpp"
#include "ads/lin/solver_ctx.hpp"


namespace ads::lin {

inline void factorize(dense_matrix& a, solver_ctx& ctx) {
    int rows = a.rows();
    int cols = a.cols();
    dgetrf_(&rows, &cols, a.data(), &ctx.lda, ctx.pivot(), &ctx.info);
}

template <typename Rhs>
inline void solve_with_factorized(const dense_matrix& a, Rhs& b, solver_ctx& ctx) {
    int nrhs = b.size() / b.size(0);
    int n = a.rows();
    dgetrs_("N", &n, &nrhs, a.data(), &ctx.lda, ctx.pivot(), b.data(), &n, &ctx.info);
}

}

#endif /* ADS_LIN_DENSE_SOLVE_HPP_ */
