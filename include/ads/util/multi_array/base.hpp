// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MULTI_ARRAY_BASE_HPP
#define ADS_UTIL_MULTI_ARRAY_BASE_HPP

#include <cstddef>
#include <type_traits>

#include "ads/util/meta.hpp"
#include "ads/util/multi_array/ordering/standard.hpp"

namespace ads {

template <                                                  //
    typename T,                                             // element type
    std::size_t Rank,                                       // number of indices
    typename Array,                                         // self type (CRTP)
    template <std::size_t> class Order = standard_ordering  // memory layout
    >
struct multi_array_base : private Order<Rank> {
public:
    using Ordering = Order<Rank>;
    using size_array = typename Ordering::size_array;

    explicit multi_array_base(const size_array& sizes)
    : Ordering{sizes} { }

    template <typename... Indices>
    const T& operator()(Indices... indices) const {
        check_indices_(indices...);
        const auto lin_idx = Ordering::linear_index(indices...);
        return static_cast<const Array*>(this)->storage_(lin_idx);
    }

    template <typename... Indices>
    T& operator()(Indices... indices) {
        check_indices_(indices...);
        const auto lin_idx = Ordering::linear_index(indices...);
        return static_cast<Array*>(this)->storage_(lin_idx);
    }

    int size(int dim) const { return Ordering::size(dim); }

    int size() const { return Ordering::size(); }

    size_array sizes() const { return Ordering::sizes(); }

private:
    template <typename... Indices>
    void check_indices_(Indices...) const {
        static_assert(sizeof...(Indices) == Rank, "Invalid number of indices");
        static_assert(util::all_<std::is_integral, Indices...>::value,
                      "Indices need to be of integral type");
    }
};

}  // namespace ads

#endif  // ADS_UTIL_MULTI_ARRAY_BASE_HPP
