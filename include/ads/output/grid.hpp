// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_GRID_HPP
#define ADS_OUTPUT_GRID_HPP

#include <tuple>

#include "ads/output/range.hpp"

namespace ads::output {

template <typename... RangeIters>
struct grid {
    static constexpr std::size_t dim = sizeof...(RangeIters);

    std::tuple<range<RangeIters>...> ranges;

    explicit grid(range<RangeIters>... ranges)
    : ranges{ranges...} { }

    template <std::size_t I>
    auto begin() const {
        return std::get<I>(ranges).begin;
    }

    template <std::size_t I>
    auto end() const {
        return std::get<I>(ranges).end;
    }
};

template <std::size_t I, typename... RangeIters>
auto get_range(const grid<RangeIters...>& g) {
    return std::get<I>(g.ranges);
}

template <std::size_t I, typename... RangeIters>
int size(const grid<RangeIters...>& g) {
    return std::get<I>(g.ranges).size();
}

template <std::size_t I, typename... RangeIters>
auto begin(const grid<RangeIters...>& g) {
    return begin(std::get<I>(g.ranges));
}

template <std::size_t I, typename... RangeIters>
auto end(const grid<RangeIters...>& g) {
    return end(std::get<I>(g.ranges));
}

namespace impl {

template <typename... Conts>
struct grid_from_containers {
    template <typename Cont>
    using range_ = decltype(from_container(std::declval<Cont>()));

    using type = grid<typename range_<Conts>::iterator...>;
};

template <std::size_t I, std::size_t Rank>
struct size_array_helper {
    using Next = size_array_helper<I + 1, Rank>;

    template <typename... RangeIters>
    static void fill_dimension(std::array<std::size_t, sizeof...(RangeIters)>& dims,
                               const grid<RangeIters...>& grid) {
        dims[I] = size<I>(grid);
        Next::fill_dimension(dims, grid);
    }
};

template <std::size_t Rank>
struct size_array_helper<Rank, Rank> {
    template <typename... RangeIters>
    static void fill_dimension(std::array<std::size_t, sizeof...(RangeIters)>&,
                               const grid<RangeIters...>&) {
        // empty
    }
};

}  // namespace impl

template <typename... RangeIters>
grid<RangeIters...> make_grid(const range<RangeIters>&... ranges) {
    return grid<RangeIters...>{ranges...};
}

template <typename... Conts>
auto grid_from_containers(const Conts&... conts) ->
    typename impl::grid_from_containers<Conts...>::type {
    return {from_container(conts)...};
}

template <typename... RangeIters>
auto dims(const grid<RangeIters...>& grid) {
    std::array<std::size_t, sizeof...(RangeIters)> dims;
    impl::size_array_helper<0, sizeof...(RangeIters)>::fill_dimension(dims, grid);
    return dims;
}

}  // namespace ads::output

#endif  // ADS_OUTPUT_GRID_HPP
