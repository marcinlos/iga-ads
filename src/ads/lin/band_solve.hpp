#ifndef ADS_LIN_BAND_SOLVE_HPP_
#define ADS_LIN_BAND_SOLVE_HPP_

#include <iostream>
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/lin/lapack.hpp"

namespace ads {
namespace lin {

struct solver_ctx {
    std::vector<int> pivot_vector;
    int info = 0;
    int ldab;

    solver_ctx(int n, int kl, int ku)
    : pivot_vector(n)
    , ldab(2 * kl + ku + 1)
    { }

    solver_ctx(const band_matrix& a)
    : solver_ctx(a.rows, a.kl, a.ku)
    { }

    int* pivot() {
        return pivot_vector.data();
    }
};

inline void factorize(band_matrix& a, solver_ctx& ctx) {
    dgbtrf_(&a.rows, &a.cols, &a.kl, &a.ku, a.full_buffer(), &ctx.ldab, ctx.pivot(), &ctx.info);
}

template <typename Rhs>
inline void solve_with_factorized(const band_matrix& a, Rhs& b, solver_ctx& ctx) {
    int nrhs = b.size() / b.size(0);
    const char* trans = "No transpose";

    dgbtrs_(trans, &a.cols, &a.kl, &a.ku, &nrhs, a.full_buffer(), &ctx.ldab, ctx.pivot(), b.data(), &a.cols, &ctx.info);
}

template <typename Rhs>
inline void solve(band_matrix& a, Rhs& b, solver_ctx& ctx) {
    factorize(a, ctx);
    solve_with_factorized(a, b, ctx);
}


}
}

#endif /* ADS_LIN_BAND_SOLVE_HPP_ */
