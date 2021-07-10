// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_TENSOR_HPP
#define ADS_LIN_TENSOR_TENSOR_HPP

#include <vector>

#include "ads/lin/tensor/base.hpp"

namespace ads::lin {

namespace detail {

template <std::size_t N>
int product(const std::array<int, N>& ns) {
    int p = 1;
    for (auto n : ns) {
        p *= n;
    }
    return p;
}

}  // namespace detail

template <typename T, std::size_t Rank>
struct tensor : tensor_base<T, Rank, tensor<T, Rank>> {
private:
    using Self = tensor<T, Rank>;
    using Base = tensor_base<T, Rank, Self>;
    using size_array = typename Base::size_array;

    std::vector<T> buffer_;

public:
    explicit tensor(const size_array& sizes)
    : Base{sizes}
    , buffer_(detail::product(sizes)) { }

    T* data() { return buffer_.data(); }

    const T* data() const { return buffer_.data(); }

    void fill_with_zeros() { buffer_.assign(buffer_.size(), T{}); }
};

template <typename T, std::size_t Rank>
void zero(tensor<T, Rank>& tensor) {
    tensor.fill_with_zeros();
}

}  // namespace ads::lin

#endif  // ADS_LIN_TENSOR_TENSOR_HPP
