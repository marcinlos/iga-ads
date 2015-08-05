#ifndef ADS_UTIL_MULTI_ARRAY_WRAPPER_HPP_
#define ADS_UTIL_MULTI_ARRAY_WRAPPER_HPP_

#include "ads/util/multi_array/ordering/standard.hpp"
#include "ads/util/multi_array/base.hpp"


namespace ads {


template <
    typename T,
    std::size_t Rank,
    typename Buffer,
    template <std::size_t> class Order = standard_ordering
>
struct multi_array_wrapper : multi_array_base<T, Rank, multi_array_wrapper<T, Rank, Buffer, Order>> {

    using Self = multi_array_wrapper<T, Rank, Buffer, Order>;
    using Base = multi_array_base<T, Rank, Self, Order>;
    using size_array = typename Base::size_array;

    Buffer data;

    multi_array_wrapper(Buffer s, const size_array& sizes)
    : Base { sizes }
    , data { s }
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

}


#endif /* ADS_UTIL_MULTI_ARRAY_WRAPPER_HPP_ */
