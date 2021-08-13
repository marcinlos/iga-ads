// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SOLVER_HPP
#define ADS_SOLVER_HPP

#include <type_traits>
#include <utility>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/util.hpp"

namespace ads {

struct dim_data {
    lin::band_matrix const& M;
    lin::solver_ctx& ctx;
};

namespace detail {

// Recursive implementation of the standard (Kronecker product) ADS.
//
// For each dimension:
// 1. Solve all the right-hand sides with a corresponding 1D matrix
// 2. Cyclically transpose dimensions to bring the next one "to the front"

template <typename Rhs, typename Buf>
auto solve_with_dim_data(Rhs&, Buf&) -> void {
    // base case
}

template <typename Rhs, typename Buf, typename Dim, typename... Dims>
auto solve_with_dim_data(Rhs& rhs, Buf& buf, Dim&& dim, Dims&&... dims) -> void {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
    auto F = lin::cyclic_transpose(rhs, buf.data());

    solve_with_dim_data(F, rhs, std::forward<Dims>(dims)...);
}

// Auxiliary functions that compute the index of the dimension requiring custom handling in the
// generalized version of the ADS

// SFINAE needed since Fun&& can be a better match for dim_data arguments
template <typename Fun, typename... Dims,  //
          std::enable_if_t<!std::is_convertible_v<Fun, dim_data>, int> = 0>
auto find_special_dim(Fun&&, Dims&&...) -> int {
    return 0;
}

template <typename... Dims>
auto find_special_dim(dim_data const&, Dims&&... dims) -> int {
    return 1 + find_special_dim(std::forward<Dims>(dims)...);
}

// Recursive implementation of the generalized ADS, where one of the dimensions receives special
// treatment - each right-hand side in this dimension can be solved with a different matrix.  In
// this implementation the user specifies what to do with these right hand sides by supplying a
// callable object that receives the entire right-hand side in this step.
//
// The generalized ADS requires that the special dimension be solved first. Transposing the RHS to
// bring it to the front is done by suitable ads_solve_impl overload, so that solve_with_special_dim
// is called with correctly reordered RHS (tensor_view), but dims are passed as given to
// ads_solve_impl, with no order adjustment. To account for this, solve_with_special_dim moves all
// the leading dim_data arguments to the end. Only when it finally encounters a suitable callable
// object as the first dimension describing object does it do its job.
//
// Right-hand side is passed in rhs, and after calling this function, rhs holds the solution. This
// solution resides in memory originally pointed to by either rhs or buf, depending on the number of
// dimensions (after the call, memory of rhs and buf may be swapped).

template <typename Rhs, typename Fun, typename... Dims,
          std::enable_if_t<std::is_invocable_v<Fun, Rhs&>, int> = 0>
auto solve_with_special_dim(Rhs& rhs, Rhs& buf, Fun&& fun, Dims&&... dims) -> void {
    fun(rhs);
    auto transposed = lin::cyclic_transpose(rhs, buf.data());

    constexpr auto N = sizeof...(Dims);

    // rhs is in the memory of buf, then solve_with_dim_data
    // swaps it N times between buf and rhs
    solve_with_dim_data(transposed, rhs, std::forward<Dims>(dims)...);

    // ensure rhs point to the memory with solution
    if (N % 2 == 0) {
        using std::swap;
        swap(rhs, buf);
    }
}

template <typename Rhs, typename... Dims>
auto solve_with_special_dim(Rhs& rhs, Rhs& buf, dim_data const& dim, Dims&&... dims) -> void {
    solve_with_special_dim(rhs, buf, std::forward<Dims>(dims)..., dim);
}

template <typename T, std::size_t Rank>
auto do_n_cyclic_transpose(lin::tensor_view<T, Rank>& rhs, T* buf, int n)
    -> lin::tensor_view<T, Rank> {
    if (n > 0) {
        auto F = lin::cyclic_transpose(rhs, buf);
        return do_n_cyclic_transpose(F, rhs.data(), n - 1);
    } else {
        return rhs;
    }
}

// Cyclically transposes rhs n times, returns a view of the transposed tensor.
//
// Transposed tensor physically resides in the memory originally pointed to by either rhs of buf.
// After the call, rhs points to the memory containing the transposed tensor and buf to the other
// memory region.
template <typename T, std::size_t Rank>
auto n_cyclic_transpose(lin::tensor_view<T, Rank>& rhs, lin::tensor_view<T, Rank>& buf, int n)
    -> lin::tensor_view<T, Rank> {
    if (n > 0) {
        auto res = do_n_cyclic_transpose(rhs, buf.data(), n);

        // do_n_cyclic_transpose copies data between memory pointed to by rhs and buf n times
        // Transposed data is physically either in memory of rhs or in memory of buf
        // Ultimately, we want rhs to point to it
        if (n % 2 == 1) {
            using std::swap;
            swap(rhs, buf);
        }
        return res;
    } else {
        return rhs;
    }
}

// Tag structs used to dispatch ads_solve to its correct implementation
struct only_dim_data { };
struct with_special_dim { };

// Chooses appropriate implementation tag based on whether all the dimension describing objects are
// of type dim_data.
template <typename... Dims>
using choose_impl = std::conditional_t<                          //
    std::conjunction_v<std::is_convertible<Dims, dim_data>...>,  //
    only_dim_data,                                               //
    with_special_dim                                             //
    >;

// Standard ADS implementation. All dims arguments are of type dim_data.
// This is only called for the number of dimensions > 1.
template <typename Rhs, typename... Dims>
auto ads_solve_impl(Rhs& rhs, Rhs& buf, only_dim_data, Dims&&... dims) -> void {
    solve_with_dim_data(rhs, buf, std::forward<Dims>(dims)...);

    constexpr auto N = sizeof...(Dims);

    // We perform transposition N times
    // Transposition swaps rhs and buffer, so for odd N we need one more swap
    if (N % 2 == 1) {
        using std::swap;
        swap(buf, rhs);
    }
}
//
// Specialized implementation that avoids needless transpositions.
template <typename Rhs>
void ads_solve_impl(Rhs& rhs, Rhs& /*buf*/, only_dim_data, dim_data const& dim) {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
}

// Generalized ADS implementation. One of the dims arguments should be a callable object that
// accepts Rhs.
template <typename Rhs, typename... Dims>
auto ads_solve_impl(Rhs& rhs, Rhs& buf, with_special_dim, Dims&&... dims) -> void {
    auto const n = detail::find_special_dim(std::forward<Dims>(dims)...);

    // solve_with_special_dim operates on tensor views
    auto rhs_view = lin::as_tensor(rhs.data(), rhs.sizes());
    auto buf_view = lin::as_tensor(buf.data(), buf.sizes());

    auto rhs_trans = n_cyclic_transpose(rhs_view, buf_view, n);
    auto buf_trans = lin::as_tensor(buf_view.data(), rhs_trans.sizes());

    solve_with_special_dim(rhs_trans, buf_trans, std::forward<Dims>(dims)...);

    // if n = 0, special dimension was first and no additional transpositions
    // were necessary before calling solve_with_special_dim
    if (n != 0) {
        constexpr auto N = sizeof...(Dims);
        n_cyclic_transpose(rhs_trans, buf_trans, narrow_cast<int>(N - n));
    }

    // transpositions move data between rhs and buf, here we ensure it ends up in rhs
    if (rhs_trans.data() != rhs.data()) {
        using std::swap;
        swap(rhs, buf);
    }
}

}  // namespace detail

// For one dimension no auxiliary buffer is necessary
template <typename Rhs>
auto ads_solve(Rhs& rhs, dim_data const& dim) -> void {
    lin::solve_with_factorized(dim.M, rhs, dim.ctx);
}

/**
 * @brief Solve the system of linear equations using ADS.
 *
 * There are two modes the ADS can operate in:
 *
 * - If each of @c dims is a @c dim_data object, the matrix of the system being solved is simply the
 *   Kronecker product of 1D matrices in @c dims (standard ADS)
 *
 * - Otherwise, exactly one of @c dims should be a callable object that accepts @c Rhs
 *
 * @tparam Rhs type of the right hand side and buffer
 * @tparam Dims types of arguments describing dimensions
 *
 * @param rhs right-hand side of the system
 * @param buf auxiliary buffer large enough to store @p rhs
 * @param dims objects describing dimensions the full matrix is decomposed into
 */
template <typename Rhs, typename... Dims>
auto ads_solve(Rhs& rhs, Rhs& buf, Dims&&... dims) -> void {
    using impl_tag = detail::choose_impl<Dims...>;
    detail::ads_solve_impl(rhs, buf, impl_tag{}, std::forward<Dims>(dims)...);
}

}  // namespace ads

#endif  // ADS_SOLVER_HPP
