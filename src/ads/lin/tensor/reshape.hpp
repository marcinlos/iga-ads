#ifndef ADS_LIN_TENSOR_RESHAPE_HPP_
#define ADS_LIN_TENSOR_RESHAPE_HPP_

#include "ads/lin/tensor/base.hpp"
#include "ads/lin/tensor/view.hpp"


namespace ads {
namespace lin {

template <typename T, std::size_t Rank, typename Impl, typename... Sizes>
tensor_view<T, sizeof...(Sizes)> reshape(tensor_base<T, Rank, Impl>& tensor, Sizes... sizes) {
    return { tensor.data(), sizes... };
}

template <typename T, std::size_t Rank, typename Impl, typename... Sizes>
tensor_view<T, sizeof...(Sizes)> reshape(const tensor_base<T, Rank, Impl>& tensor, Sizes... sizes) {
    return { tensor.data(), sizes... };
}


}
}


#endif /* ADS_LIN_TENSOR_RESHAPE_HPP_ */
