#ifndef ADS_LIN_TENSOR_FOR_EACH_HPP_
#define ADS_LIN_TENSOR_FOR_EACH_HPP_

#include <functional>

#include "ads/lin/tensor/base.hpp"


namespace ads::lin {

namespace impl {

template <
    typename T,
    std::size_t Rank,
    std::size_t I,
    typename Tensor,
    typename... Indices
>
struct tensor_helper {

    using Next = tensor_helper<T, Rank, I + 1, Tensor, std::size_t, Indices...>;

    template <typename F>
    static void for_each_multiindex(F&& fun, const Tensor& t, Indices... indices) {
        auto size = t.size(I);
        for (std::size_t i = 0; i < size; ++ i) {
            Next::for_each_multiindex(std::forward<F>(fun), t, i, indices...);
        }
    }
};

template <
    typename T,
    std::size_t Rank,
    typename Tensor,
    typename... Indices
>
struct tensor_helper<T, Rank, Rank, Tensor, Indices...> {

    template <typename F>
    static void for_each_multiindex(F&& fun, const Tensor&, Indices... indices) {
        fun(indices...);
    }
};

}


template <typename T, std::size_t Rank, typename Impl, typename F>
void for_each_multiindex(F&& fun, const lin::tensor_base<T, Rank, Impl>& t) {
    using tensor_type = lin::tensor_base<T, Rank, Impl>;
    using helper_type = impl::tensor_helper<T, Rank, 0, tensor_type>;

    helper_type::for_each_multiindex(std::forward<F>(fun), t);
}

}

#endif /* ADS_LIN_TENSOR_FOR_EACH_HPP_ */
