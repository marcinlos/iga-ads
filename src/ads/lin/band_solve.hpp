#ifndef ADS_LIN_BAND_SOLVE_HPP_
#define ADS_LIN_BAND_SOLVE_HPP_

#include <iostream>
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/tensor.hpp"

// LAPACK routines
extern "C" {

using in_int = const int*;
using in_int_array = const int*;
using out_int = int*;
using out_int_array = int*;

int dgbtrf_(
        in_int m,
        in_int n,
        in_int kl,
        in_int ku,
        double* ab,
        in_int ldab,
        out_int_array ipiv,
        out_int info);

int dgbtrs_(
        const char* trans,
        in_int n,
        in_int kl,
        in_int ku,
        in_int nrhs,
        const double* ab,
        in_int ldab,
        in_int_array ipiv,
        double* b,
        in_int ldb,
        out_int info);
}

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
    : solver_ctx(a.n, a.kl, a.ku)
    { }

    int* pivot() {
        return pivot_vector.data();
    }
};

inline void factorize(band_matrix& a, solver_ctx& ctx) {
    dgbtrf_(&a.n, &a.n, &a.kl, &a.ku, a.data(), &ctx.ldab, ctx.pivot(), &ctx.info);
}

template <typename Rhs>
inline void solve_with_factorized(const band_matrix& a, Rhs& b, solver_ctx& ctx) {
    int nrhs = b.size() / b.size(0);
    const char* trans = "No transpose";

    dgbtrs_(trans, &a.n, &a.kl, &a.ku, &nrhs, a.data(), &ctx.ldab, ctx.pivot(), b.data(), &a.n, &ctx.info);
}

template <typename Rhs>
inline void solve(band_matrix& a, Rhs& b, solver_ctx& ctx) {
    factorize(a, ctx);
    solve_with_factorized(a, b, ctx);
}


}
}

#endif /* ADS_LIN_BAND_SOLVE_HPP_ */
