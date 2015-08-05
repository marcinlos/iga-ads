#ifndef ADS_SOLVER_HPP_
#define ADS_SOLVER_HPP_

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"


namespace ads {

struct dim_data {
    const lin::band_matrix& M;
    lin::solver_ctx& ctx;
};

template <typename Rhs>
void ads_solve(Rhs& rhs, const dim_data& dim) {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
}

template <typename Rhs>
void ads_solve(Rhs& rhs, Rhs& buf, const dim_data& dim1, const dim_data& dim2) {
    double* buffer2 = buf.data();

    lin::solve_with_factorized(dim1.M, rhs, dim1.ctx);
    auto F2 = lin::cyclic_transpose(rhs, buffer2);

    lin::solve_with_factorized(dim2.M, F2, dim2.ctx);
    lin::cyclic_transpose(F2, rhs);

    using std::swap;
    swap(buf, rhs);
}


template <typename Rhs>
void ads_solve(Rhs& rhs, Rhs& buf, const dim_data& dim1, const dim_data& dim2, const dim_data& dim3) {
    double* buffer1 = rhs.data();
    double* buffer2 = buf.data();

    lin::solve_with_factorized(dim1.M, rhs, dim1.ctx);
    auto F2 = lin::cyclic_transpose(rhs, buffer2);

    lin::solve_with_factorized(dim2.M, F2, dim2.ctx);
    auto F3 = lin::cyclic_transpose(F2, buffer1);

    lin::solve_with_factorized(dim3.M, F3, dim3.ctx);
    lin::cyclic_transpose(F3, buf);

    using std::swap;
    swap(buf, rhs);
}


}


#endif /* ADS_SOLVER_HPP_ */
