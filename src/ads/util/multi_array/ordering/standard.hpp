#ifndef ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP_
#define ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP_

#include <cstddef>
#include <cassert>
#include "ads/util/meta.hpp"


namespace ads {

namespace impl {

template <std::size_t I, std::size_t Rank>
struct standard_ordering_indexer_ : standard_ordering_indexer_<I + 1, Rank> {

    static_assert(I >= 0, "Negative index");
    static_assert(I < Rank, "Index larger than rank");

    using Base = standard_ordering_indexer_<I + 1, Rank>;

    std::size_t n;

    template <typename... Sizes>
    standard_ordering_indexer_(std::size_t n, Sizes... extents)
    : Base(extents...)
    , n(n)
    { }

    template <typename... Indices>
    std::size_t linearize(std::size_t base, std::size_t idx, Indices... indices) const {
        assert(idx < n && "Index out of bounds");
        return Base::linearize(base * n + idx, indices...);
    }

    std::size_t size(std::size_t dim) const {
        return I == dim ? n : Base::size(dim);
    }

    std::size_t size() const {
        return n * Base::size();
    }
};

template <std::size_t Rank>
struct standard_ordering_indexer_<Rank, Rank> {

    std::size_t linearize(std::size_t base) const {
        return base;
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
struct standard_ordering : private impl::standard_ordering_indexer_<0, Rank> {

    using Indexer = impl::standard_ordering_indexer_<0, Rank>;

    template <typename... Sizes>
    standard_ordering(Sizes... sizes)
    : Indexer(sizes...)
    { }

    template <typename... Indices>
    std::size_t linear_index(Indices... indices) const {
        static_assert(util::all_<std::is_integral, Indices...>::value, "Indices need to be of integral type");
        return Indexer::linearize(0, indices...);
    }

    std::size_t size(std::size_t dim) const {
        assert(dim < Rank && "Index larger than rank");
        return Indexer::size(dim);
    }
};


}


#endif /* ADS_UTIL_MULTI_ARRAY_ORDERING_STANDARD_HPP_ */
