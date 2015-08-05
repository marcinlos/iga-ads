#ifndef ADS_LIN_BAND_SOLVE_HPP_
#define ADS_LIN_BAND_SOLVE_HPP_

#include <iostream>
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/tensor.hpp"

// LAPACK routines
extern "C" {

int dgbtrf_(int* m, int* n, int* kl, int* ku, double* ab, int* ldab, int* ipiv, int* info);

int dgbtrs_(const char* trans, int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab,
        int* ipiv, double* b, int* ldb, int* info);
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
inline void solve_with_factorized(band_matrix& a, Rhs& b, solver_ctx& ctx) {
    int nrhs = b.size(1);
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
