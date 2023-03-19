// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP
#define ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP

#include "ads/util/multi_array/wrapper.hpp"

namespace ads {

template <                                                       //
    std::size_t DestDim,                                         //
    template <std::size_t> class DestOrder = standard_ordering,  //
    typename T,                                                  //
    std::size_t D,                                               //
    typename Buffer,                                             //
    template <std::size_t> class Order,                          //
    typename... Sizes                                            //
    >
multi_array_wrapper<T, DestDim, Buffer, DestOrder>
reshape(multi_array_wrapper<T, D, Buffer, Order> a, Sizes... sizes) {
    return {a.data, {sizes...}};
}

}  // namespace ads

#endif  // ADS_UTIL_MULTI_ARRAY_RESHAPE_HPP
