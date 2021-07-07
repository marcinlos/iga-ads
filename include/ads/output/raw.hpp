// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_RAW_HPP
#define ADS_OUTPUT_RAW_HPP

#include <ostream>

#include "ads/output/output_base.hpp"
#include "ads/output/output_format.hpp"
#include "ads/util/io.hpp"

namespace ads::output {

struct raw_printer : output_base {
    using Base = output_base;

    explicit raw_printer(const output_format& format)
    : Base{format} { }

    template <typename... Data>
    void print(std::ostream& os, std::size_t count, const Data&... data) const {
        ads::util::stream_state_saver guard(os);
        prepare_stream(os);
        for (std::size_t i = 0; i < count; ++i) {
            print_row(os, data(i)...);
        }
    }
};

}  // namespace ads::output

#endif  // ADS_OUTPUT_RAW_HPP
