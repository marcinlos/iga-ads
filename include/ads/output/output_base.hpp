// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_OUTPUT_BASE_HPP
#define ADS_OUTPUT_OUTPUT_BASE_HPP

#include "ads/output/output_format.hpp"


namespace ads::output {

struct output_base {

    output_format format;

    explicit output_base(const output_format& format)
    : format { format }
    { }

    void prepare_stream(std::ostream& os) const {
        os.precision(format.precision());
        os.setf(format.flags(), format.mask());
    }

    template <typename T>
    void print_one(T value, std::ostream& os) const {
        os.width(format.width());
        os << value;
    }

    template <typename Value, typename... Values>
    void print_row(std::ostream& os, Value x, Values... values) const {
        print_one(x, os);
        print_row(os, values...);
    }

    void print_row(std::ostream& os) const {
        os << '\n';
    }
};

}

#endif // ADS_OUTPUT_OUTPUT_BASE_HPP
