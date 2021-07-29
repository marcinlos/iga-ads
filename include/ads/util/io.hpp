// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_IO_HPP
#define ADS_UTIL_IO_HPP

#include <ostream>

namespace ads::util {

struct stream_state_saver {
    std::ios_base& os;
    std::streamsize orig_precision;
    std::ios::fmtflags orig_flags;

    explicit stream_state_saver(std::ostream& os)
    : os(os)
    , orig_precision(os.precision())
    , orig_flags(os.flags()) { }

    ~stream_state_saver() {
        os.precision(orig_precision);
        os.flags(orig_flags);
    }

    stream_state_saver(const stream_state_saver&) = delete;
    stream_state_saver& operator=(const stream_state_saver&) = delete;
    stream_state_saver(stream_state_saver&&) = delete;
    stream_state_saver& operator=(stream_state_saver&&) = delete;
};

}  // namespace ads::util

#endif  // ADS_UTIL_IO_HPP
