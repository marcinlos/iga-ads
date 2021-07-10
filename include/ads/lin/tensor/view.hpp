// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_VIEW_HPP
#define ADS_LIN_TENSOR_VIEW_HPP

#include "ads/lin/tensor/base.hpp"

namespace ads::lin {

template <typename T, std::size_t Rank>
struct tensor_view : tensor_base<T, Rank, tensor_view<T, Rank>> {
private:
    using Self = tensor_view<T, Rank>;
    using Base = tensor_base<T, Rank, Self>;
    using size_array = typename Base::size_array;

    T* data_;

public:
    tensor_view(T* data, const size_array& sizes)
    : Base{sizes}
    , data_{data} { }

    T* data() { return data_; }

    const T* data() const { return data_; }
};

template <typename T, typename... Sizes>
tensor_view<T, sizeof...(Sizes)> as_tensor(T* data, Sizes... sizes) {
    constexpr auto N = sizeof...(Sizes);
    std::array<int, N> sz{sizes...};
    return {data, sz};
}

}  // namespace ads::lin

#endif  // ADS_LIN_TENSOR_VIEW_HPP
