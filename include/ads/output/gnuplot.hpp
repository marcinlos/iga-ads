// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_GNUPLOT_HPP
#define ADS_OUTPUT_GNUPLOT_HPP

#include <iterator>

#include "ads/output/grid.hpp"
#include "ads/output/output_base.hpp"
#include "ads/output/raw.hpp"
#include "ads/util.hpp"

namespace ads::output {

template <std::size_t Dim>
struct gnuplot_printer;

template <>
struct gnuplot_printer<1> : output_base {
    static constexpr std::size_t Dim = 1;

    explicit gnuplot_printer(const output_format& format)
    : output_base{format} { }

    template <typename RangeIter, typename... Values>
    void print(std::ostream& os, const range<RangeIter>& xrange, const Values&... values) {
        print(os, make_grid(xrange), values...);
    }

    template <typename Iter, typename... Values>
    void print(std::ostream& os, const grid<Iter>& grid, const Values&... values) {
        raw_printer printer{format};
        printer.print(os, get_range<0>(grid).size(), values...);
    }
};

template <>
struct gnuplot_printer<2> : output_base {
    static constexpr std::size_t Dim = 2;

    explicit gnuplot_printer(const output_format& format)
    : output_base{format} { }

    template <typename Iter1, typename Iter2, typename... Values>
    void print(std::ostream& os, const grid<Iter1, Iter2>& grid, const Values&... values) {
        ads::util::stream_state_saver guard(os);
        prepare_stream(os);
        const auto& xs = get_range<0>(grid);
        const auto& ys = get_range<1>(grid);
        for (int i = 0; i < as_signed(xs.size()); ++i) {
            for (int j = 0; j < as_signed(ys.size()); ++j) {
                print_row(os, xs[i], ys[j], values(i, j)...);
            }
        }
    }
};

}  // namespace ads::output

#endif  // ADS_OUTPUT_GNUPLOT_HPP
