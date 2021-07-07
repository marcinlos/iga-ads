// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_BASE_HPP
#define ADS_LIN_TENSOR_BASE_HPP

#include <cstddef>

#include "ads/util/multi_array.hpp"

namespace ads::lin {

template <typename T, std::size_t Rank, typename Impl>
struct tensor_base : multi_array_base<T, Rank, tensor_base<T, Rank, Impl>, reverse_ordering> {
private:
    using Self = tensor_base<T, Rank, Impl>;
    using Base = multi_array_base<T, Rank, Self, reverse_ordering>;

protected:
    using size_array = std::array<std::size_t, Rank>;

public:
    explicit tensor_base(const size_array& sizes)
    : Base{sizes} { }

    T* data() { return self_()->data(); }

    const T* data() const { return self_()->data(); }

private:
    friend Base;

    T& storage_(std::size_t idx) { return data()[idx]; }

    const T& storage_(std::size_t idx) const { return data()[idx]; }

    const Impl* self_() const { return static_cast<const Impl*>(this); }

    Impl* self_() { return static_cast<Impl*>(this); }
};

}  // namespace ads::lin

#endif  // ADS_LIN_TENSOR_BASE_HPP
