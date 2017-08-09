#ifndef ADS_LIN_SOLVER_CTX_HPP_
#define ADS_LIN_SOLVER_CTX_HPP_

#include <vector>
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/dense_matrix.hpp"

namespace ads {
namespace lin {

struct solver_ctx {
    std::vector<int> pivot_vector;
    int info = 0;
    int lda;

    solver_ctx(int n, int lda)
    : pivot_vector(n)
    , lda(lda)
    { }

    solver_ctx(const band_matrix& a)
    : solver_ctx(a.rows, 2 * a.kl +  a.ku + 1)
    { }

    solver_ctx(const dense_matrix& a)
    : solver_ctx(a.rows(), a.rows())
    { }

    int* pivot() {
        return pivot_vector.data();
    }
};

}
}

#endif /* ADS_LIN_SOLVER_CTX_HPP_ */
