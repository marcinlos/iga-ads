// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP
#define ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP

#include <cassert>
#include <cstddef>

#include "ads/util/meta.hpp"

namespace ads {

namespace impl {

template <std::size_t I, std::size_t Rank>
struct reverse_ordering_indexer_ : reverse_ordering_indexer_<I - 1, Rank> {
    static_assert(I > 0, "Negative index");

    using Base = reverse_ordering_indexer_<I - 1, Rank>;
    using size_array = typename Base::size_array;

    explicit reverse_ordering_indexer_(size_array sizes)
    : Base{sizes} { }

    template <typename... Indices>
    int linearize(int idx, Indices... indices) const {
        assert(idx < n() && "Index out of bounds");
        return idx + n() * Base::linearize(indices...);
    }

    int size() const { return n() * Base::size(); }

private:
    int n() const { return std::get<Rank - I>(Base::sizes); }
};

template <std::size_t Rank>
struct reverse_ordering_indexer_<0, Rank> {
    using size_array = std::array<int, Rank>;

    size_array sizes;

    explicit reverse_ordering_indexer_(size_array sizes)
    : sizes(sizes) { }

    int linearize() const { return 0; }

    int size() const { return 1; }
};

}  // namespace impl

template <std::size_t Rank>
struct reverse_ordering : private impl::reverse_ordering_indexer_<Rank, Rank> {
    using Indexer = impl::reverse_ordering_indexer_<Rank, Rank>;
    using size_array = typename Indexer::size_array;

    explicit reverse_ordering(size_array sizes)
    : Indexer{sizes} { }

    template <typename... Indices>
    int linear_index(Indices... indices) const {
        static_assert(util::all_<std::is_integral, Indices...>::value,
                      "Indices need to be of integral type");
        return Indexer::linearize(indices...);
    }

    int size(int dim) const {
        assert(dim < Rank && "Index larger than rank");
        return Indexer::sizes[dim];
    }

    int size() const { return Indexer::size(); }

    size_array sizes() const { return Indexer::sizes; }
};

}  // namespace ads

#endif  // ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP
