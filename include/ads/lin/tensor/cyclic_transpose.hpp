// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP
#define ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP

#include <cstddef>

#include "ads/lin/tensor/base.hpp"
#include "ads/lin/tensor/view.hpp"

namespace ads::lin {

namespace detail {

template <typename T, typename S, std::size_t Rank, typename Impl1, typename Impl2, std::size_t I,
          typename... Indices>
struct cyclic_transpose_helper {
    using Input = const tensor_base<T, Rank, Impl1>&;
    using Output = tensor_base<S, Rank, Impl2>&;

    using Next = cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, I + 1, int, Indices...>;
    static constexpr std::size_t Index = Rank - I - 1;

    static void do_transpose(Input a, Output b, Indices... indices) {
        for (int i = 0; i < a.size(Index); ++i) {
            Next::do_transpose(a, b, i, indices...);
        }
    }
};

template <typename T, typename S, std::size_t Rank, typename Impl1, typename Impl2,
          typename... Indices>
struct cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, Rank, int, Indices...> {
    using Input = const tensor_base<T, Rank, Impl1>&;
    using Output = tensor_base<S, Rank, Impl2>&;

    static void do_transpose(Input a, Output b, int i, Indices... indices) {
        b(indices..., i) = a(i, indices...);
    }
};

template <std::size_t N>
auto cyclic_transpose_sizes(const std::array<int, N>& sizes) {
    // TODO: use std::rotate
    std::array<int, N> new_sizes;
    for (std::size_t i = 0; i < N - 1; ++i) {
        new_sizes[i] = sizes[i + 1];
    }
    new_sizes[N - 1] = sizes[0];
    return new_sizes;
}

}  // namespace detail

template <typename T, typename S, std::size_t Rank, typename Impl1, typename Impl2>
void cyclic_transpose(const tensor_base<T, Rank, Impl1>& a, tensor_base<S, Rank, Impl2>& out) {
    detail::cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, 0>::do_transpose(a, out);
}

template <typename T, std::size_t Rank, typename Impl>
tensor_view<T, Rank> cyclic_transpose(const tensor_base<T, Rank, Impl>& a, T* out) {
    tensor_view<T, Rank> view{out, detail::cyclic_transpose_sizes(a.sizes())};
    detail::cyclic_transpose_helper<T, T, Rank, Impl, decltype(view), 0>::do_transpose(a, view);
    return view;
}

}  // namespace ads::lin

#endif  // ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP
