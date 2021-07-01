// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SOLVER_HPP_
#define ADS_SOLVER_HPP_

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"


namespace ads {

struct dim_data {
    const lin::band_matrix& M;
    lin::solver_ctx& ctx;
};

template <typename Rhs, typename Buf, typename Dim, typename... Dims>
void solve_worker(Rhs& rhs, Buf& buf, const Dim& dim, const Dims&... dims) {

    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
    auto F = lin::cyclic_transpose(rhs, buf.data());

    solve_worker(F, rhs, dims...);
}

template <typename Rhs, typename Buf>
void solve_worker(Rhs&, Buf&) {
    // base case
}


template <typename Rhs>
void ads_solve(Rhs& rhs, dim_data dim) {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
}


template <typename Rhs>
void ads_solve(Rhs& rhs, Rhs& /*buf*/, dim_data dim) {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
}

template <typename Rhs, typename... Dims>
void ads_solve(Rhs& rhs, Rhs& buf, Dims... dims) {
    solve_worker(rhs, buf, dims...);

    constexpr std::size_t N = sizeof...(Dims);

    // for N > 1 we perform transposition N times
    // transposition swaps rhs and buffer, so for odd N we need one more swap
    if (N % 2 == 1) {
        using std::swap;
        swap(buf, rhs);
    }
}

template <typename Rhs, std::size_t D, std::size_t... Idx>
void ads_solve_aux(Rhs& rhs, Rhs& buf, const std::array<dim_data, D>& dims,
                   std::index_sequence<Idx...>) {
    ads_solve(rhs, buf, dims[Idx]...);
}


template <typename Rhs, std::size_t D>
void ads_solve(Rhs& rhs, Rhs& buf, const std::array<dim_data, D>& dims) {
    ads_solve_aux(rhs, buf, dims, std::make_index_sequence<D>{});
}

}

#endif /* ADS_SOLVER_HPP_ */
