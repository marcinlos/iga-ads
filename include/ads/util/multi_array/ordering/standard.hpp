// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP
#define ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP

#include <array>
#include <cassert>
#include <cstddef>

#include "ads/util.hpp"
#include "ads/util/meta.hpp"

namespace ads {

namespace impl {

template <std::size_t I, std::size_t Rank>
struct standard_ordering_indexer_ : standard_ordering_indexer_<I + 1, Rank> {
    static_assert(I >= 0, "Negative index");
    static_assert(I < Rank, "Index larger than rank");

    using Base = standard_ordering_indexer_<I + 1, Rank>;
    using size_array = typename Base::size_array;

    explicit standard_ordering_indexer_(size_array sizes)
    : Base{sizes} { }

    template <typename... Indices>
    int linearize(int base, int idx, Indices... indices) const {
        assert(idx < n() && "Index out of bounds");
        return Base::linearize(base * n() + idx, indices...);
    }

    int size() const { return n() * Base::size(); }

private:
    int n() const { return std::get<I>(Base::sizes); }
};

template <std::size_t Rank>
struct standard_ordering_indexer_<Rank, Rank> {
    using size_array = std::array<int, Rank>;

    size_array sizes;

    explicit standard_ordering_indexer_(size_array sizes)
    : sizes(sizes) { }

    int linearize(int base) const { return base; }

    int size() const { return 1; }
};

}  // namespace impl

template <std::size_t Rank>
struct standard_ordering : private impl::standard_ordering_indexer_<0, Rank> {
    using Self = standard_ordering<Rank>;
    using Indexer = impl::standard_ordering_indexer_<0, Rank>;
    using size_array = typename Indexer::size_array;

    explicit standard_ordering(size_array sizes)
    : Indexer{sizes} { }

    template <typename... Indices>
    int linear_index(Indices... indices) const {
        static_assert(util::all_<std::is_integral, Indices...>::value,
                      "Indices need to be of integral type");
        return Indexer::linearize(0, indices...);
    }

    int size(int dim) const {
        assert(dim < as_signed(Rank) && "Index larger than rank");
        return Indexer::sizes[dim];
    }

    int size() const { return Indexer::size(); }

    size_array sizes() const { return Indexer::sizes(); }
};

}  // namespace ads

#endif  // ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP
