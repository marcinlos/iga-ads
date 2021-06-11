#ifndef ADS_OUTPUT_RAW_HPP_
#define ADS_OUTPUT_RAW_HPP_

#include <ostream>
#include "ads/util/io.hpp"
#include "ads/output/output_format.hpp"
#include "ads/output/output_base.hpp"

namespace ads::output {

struct raw_printer : output_base {

    using Base = output_base;

    explicit raw_printer(const output_format& format)
    : Base { format }
    { }

    template <typename... Data>
    void print(std::ostream& os, std::size_t count, const Data&... data) const {
        ads::util::stream_state_saver guard(os);
        prepare_stream(os);
        for (std::size_t i = 0; i < count; ++ i) {
            print_row(os, data(i)...);
        }
    }

};

}

#endif /* ADS_OUTPUT_RAW_HPP_ */
