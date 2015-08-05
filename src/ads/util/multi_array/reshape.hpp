#ifndef ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP_
#define ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP_

#include "ads/util/multi_array/wrapper.hpp"


namespace ads {

template <
    std::size_t DestDim,
    template <std::size_t> class DestOrder = standard_ordering,
    typename T,
    std::size_t D,
    typename Buffer,
    template <std::size_t> class Order,
    typename... Sizes
>
multi_array_wrapper<T, DestDim, Buffer, DestOrder> reshape(multi_array_wrapper<T, D, Buffer, Order> a, Sizes... sizes) {
    return { a.data, { sizes... } };
}


}


#endif /* ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP_ */
