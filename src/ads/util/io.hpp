#ifndef ADS_UTIL_IO_HPP_
#define ADS_UTIL_IO_HPP_

#include <ostream>

namespace ads::util {

struct stream_state_saver {
    std::ios_base& os;
    std::streamsize orig_precision;
    std::ios::fmtflags orig_flags;

    explicit stream_state_saver(std::ostream& os)
    : os(os)
    , orig_precision(os.precision())
    , orig_flags(os.flags())
    { }

    ~stream_state_saver() {
        os.precision(orig_precision);
        os.flags(orig_flags);
    }
};

}


#endif /* ADS_UTIL_IO_HPP_ */
