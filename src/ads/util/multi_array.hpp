#ifndef ADS_UTIL_MULTI_ARRAY_HPP_
#define ADS_UTIL_MULTI_ARRAY_HPP_

#include <ads/util/meta.hpp>
#include <cassert>

namespace ads {

//-------------------------------------------------------------------------------------------------
// Standard ordering, linear indexing
//-------------------------------------------------------------------------------------------------

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

template <std::size_t Rank>
struct standard_ordering : private standard_ordering_indexer_<0, Rank> {

    using Indexer = standard_ordering_indexer_<0, Rank>;

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


//-------------------------------------------------------------------------------------------------
// Reverse ordering (first index varies most rapidly), linear indexing
//-------------------------------------------------------------------------------------------------

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

template <std::size_t Rank>
struct reverse_ordering : private reverse_ordering_indexer_<Rank> {

    using Indexer = reverse_ordering_indexer_<Rank>;

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


//-------------------------------------------------------------------------------------------------
// Multidimensional array base class
//-------------------------------------------------------------------------------------------------

template <
    typename T,
    std::size_t Rank,
    typename Array,
    template <std::size_t> class Order = standard_ordering
>
struct multi_array_base: private Order<Rank> {
public:
    using Ordering = Order<Rank>;
    using Self = multi_array_base<T, Rank, Array>;

    template <typename... Sizes>
    multi_array_base(Sizes... sizes)
    : Ordering(sizes...)
    {
        static_assert(sizeof...(sizes) == Rank, "Invalid number of dimension sizes passed");
        static_assert(util::all_<std::is_integral, Sizes...>::value, "Sizes need to be of integral type");
    }

    template <typename... Indices>
    const T& operator ()(Indices... indices) const {
        check_indices_(indices...);
        std::size_t lin_idx = Ordering::linear_index(indices...);
        return static_cast<const Array*>(this)->storage_(lin_idx);
    }

    template <typename... Indices>
    T& operator ()(Indices... indices) {
        check_indices_(indices...);
        std::size_t lin_idx = Ordering::linear_index(indices...);
        return static_cast<Array*>(this)->storage_(lin_idx);
    }

    std::size_t size(std::size_t dim) const {
        return Ordering::size(dim);
    }

private:
    template <typename... Indices>
    void check_indices_(Indices...) const {
        static_assert(sizeof...(Indices) == Rank, "Invalid number of indices");
        static_assert(util::all_<std::is_integral, Indices...>::value, "Indices need to be of integral type");
    }
};


template <
    typename T,
    std::size_t Rank,
    typename Buffer,
    template <std::size_t> class Order = standard_ordering
>
struct multi_array_wrapper : multi_array_base<T, Rank, multi_array_wrapper<T, Rank, Buffer, Order>> {

    using Self = multi_array_wrapper<T, Rank, Buffer, Order>;
    using Base = multi_array_base<T, Rank, Self, Order>;

    Buffer data;

    template <typename... Sizes>
    multi_array_wrapper(Buffer s, Sizes... sizes)
    : Base(sizes...)
    , data(s)
    { }

private:
    friend Base;

    T& storage_(std::size_t idx) {
        return data[idx];
    }

    const T& storage_(std::size_t idx) const {
        return data[idx];
    }
};

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
    return { a.data, sizes... };
}


}



#endif /* ADS_UTIL_MULTI_ARRAY_HPP_ */
