#ifndef ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP_
#define ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP_

#include <cstddef>
#include "ads/lin/tensor/base.hpp"
#include "ads/lin/tensor/view.hpp"


namespace ads {
namespace lin {


namespace impl {

template <
    typename T,
    typename S,
    std::size_t Rank,
    typename Impl1,
    typename Impl2,
    std::size_t I,
    typename... Indices
>
struct cyclic_transpose_helper {

    using Input = const tensor_base<T, Rank, Impl1>&;
    using Output = tensor_base<S, Rank, Impl2>&;

    using Next = cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, I + 1, std::size_t, Indices...>;
    static constexpr std::size_t Index = Rank - I - 1;

    static void do_transpose(Input a, Output b, Indices... indices) {
        for (std::size_t i = 0; i < a.size(Index); ++ i) {
            Next::do_transpose(a, b, i, indices...);
        }
    }

};

template <
    typename T,
    typename S,
    std::size_t Rank,
    typename Impl1,
    typename Impl2,
    typename... Indices
>
struct cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, Rank, std::size_t, Indices...> {

    using Input = const tensor_base<T, Rank, Impl1>&;
    using Output = tensor_base<S, Rank, Impl2>&;

    static void do_transpose(Input a, Output b, size_t i, Indices... indices) {
        b(indices..., i) = a(i, indices...);
    }

};

}


template <typename T, typename S, std::size_t Rank, typename Impl1, typename Impl2>
void cyclic_transpose(const tensor_base<T, Rank, Impl1>& a, tensor_base<S, Rank, Impl2>& out) {
    impl::cyclic_transpose_helper<T, S, Rank, Impl1, Impl2, 0>::do_transpose(a, out);
}


template <typename T, std::size_t Rank, typename Impl>
tensor_view<T, Rank> cyclic_transpose(const tensor_base<T, Rank, Impl>& a, double* out) {
    tensor_view<T, Rank> view { out };
    impl::cyclic_transpose_helper<T, T, Rank, Impl, decltype(view), 0>::do_transpose(a, out);
    return view;
}



}
}


#endif /* ADS_LIN_TENSOR_CYCLIC_TRANSPOSE_HPP_ */
