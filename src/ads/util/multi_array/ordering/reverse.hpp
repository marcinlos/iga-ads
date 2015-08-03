#ifndef ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP_
#define ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP_

#include <cstddef>
#include <cassert>
#include "ads/util/meta.hpp"


namespace ads {

namespace impl {

template <std::size_t Rank>
struct reverse_ordering_indexer_: reverse_ordering_indexer_<Rank - 1> {

    static_assert(Rank > 0, "Negative index");

    using Base = reverse_ordering_indexer_<Rank - 1>;

    std::size_t n;

    template <typename... Sizes>
    reverse_ordering_indexer_(std::size_t n, Sizes... extents)
    : Base(extents...)
    , n(n)
    { }

    template <typename... Indices>
    std::size_t linearize(std::size_t idx, Indices... indices) const {
        assert(idx < n && "Index out of bounds");
        return idx + n * Base::linearize(indices...);
    }

    std::size_t size(std::size_t dim) const {
        return Rank == dim ? n : Base::size(dim);
    }

    std::size_t size() const {
        return n * Base::size();
    }
};


template <>
struct reverse_ordering_indexer_<0> {

    std::size_t linearize() const {
        return 0;
    }

    std::size_t size(std::size_t) const {
        assert(false && "Impossible");
        return 0;
    }

    std::size_t size() const {
        return 1;
    }
};

}


template <std::size_t Rank>
struct reverse_ordering : private impl::reverse_ordering_indexer_<Rank> {

    using Indexer = impl::reverse_ordering_indexer_<Rank>;

    template <typename... Sizes>
    reverse_ordering(Sizes... sizes)
    : Indexer(sizes...)
    { }

    template <typename... Indices>
    std::size_t linear_index(Indices... indices) const {
        static_assert(util::all_<std::is_integral, Indices...>::value, "Indices need to be of integral type");
        return Indexer::linearize(indices...);
    }

    std::size_t size(std::size_t dim) const {
        assert(dim < Rank && "Index larger than rank");
        return Indexer::size(Rank - dim);
    }
};

}


#endif /* ADS_UTIL_MULTI_ARRAY_ORDERING_REVERSE_HPP_ */
