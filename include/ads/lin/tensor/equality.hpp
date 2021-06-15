#ifndef ADS_LIN_TENSOR_EQUALITY_HPP_
#define ADS_LIN_TENSOR_EQUALITY_HPP_

#include <cmath>
#include <cstddef>

#include "ads/lin/tensor/base.hpp"


namespace ads::lin {

namespace impl {

template <
    typename T,
    typename S,
    typename Eps,
    std::size_t Rank,
    typename Impl1,
    typename Impl2,
    std::size_t I,
    typename... Indices
>
struct equal_helper {

    using Tensor1 = const tensor_base<T, Rank, Impl1>&;
    using Tensor2 = const tensor_base<S, Rank, Impl2>&;

    using Next = equal_helper<T, S, Eps, Rank, Impl1, Impl2, I + 1, Indices..., std::size_t>;

    static bool approx_equal(Tensor1 a, Tensor2 b, Eps eps, Indices... indices) {
        for (std::size_t i = 0; i < a.size(I); ++ i) {
            if (! Next::approx_equal(a, b, eps, indices..., i)) {
                return false;
            }
        }
        return true;
    }

};

template <
    typename T,
    typename S,
    typename Eps,
    std::size_t Rank,
    typename Impl1,
    typename Impl2,
    typename... Indices
>
struct equal_helper<T, S, Eps, Rank, Impl1, Impl2, Rank, Indices...> {

    using Tensor1 = const tensor_base<T, Rank, Impl1>&;
    using Tensor2 = const tensor_base<S, Rank, Impl2>&;

    static bool approx_equal(Tensor1 a, Tensor2 b, Eps eps, Indices... indices) {
        return std::abs(a(indices...) - b(indices...)) <= eps;
    }
};

}


template <typename T, typename S, typename Eps, std::size_t Rank, typename Impl1, typename Impl2>
inline bool approx_equal(const tensor_base<T, Rank, Impl1>& a, const tensor_base<S, Rank, Impl2>& b, Eps eps) {
    for (std::size_t i = 0; i < Rank; ++ i) {
        if (a.size(i) != b.size(i)) {
            return false;
        }
    }
    return impl::equal_helper<T, S, Eps, Rank, Impl1, Impl2, 0>::approx_equal(a, b, eps);
}

template <typename T, typename S, std::size_t Rank, typename Impl1, typename Impl2>
inline bool operator ==(const tensor_base<T, Rank, Impl1>& a, const tensor_base<S, Rank, Impl2>& b) {
    return approx_equal(a, b, T { 0 });
}

}

#endif /* ADS_LIN_TENSOR_EQUALITY_HPP_ */
