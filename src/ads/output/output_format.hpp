
#ifndef ADS_OUTPUT_OUTPUT_FORMAT_HPP_
#define ADS_OUTPUT_OUTPUT_FORMAT_HPP_

namespace ads {
namespace output {

struct output_format {
private:
    std::streamsize precision_;
    std::streamsize width_;
    std::ios_base::fmtflags flags_ = static_cast<std::ios_base::fmtflags>(0);
    std::ios_base::fmtflags mask_ = static_cast<std::ios_base::fmtflags>(0);

public:
    output_format(std::streamsize precision, std::streamsize width)
    : precision_ { precision }
    , width_ { width }
    { }

    output_format& precision(std::streamsize precision) {
        this->precision_ = precision;
        return *this;
    }

    std::streamsize precision() const {
        return precision_;
    }

    output_format& width(std::streamsize width) {
        this->width_ = width;
        return *this;
    }

    std::streamsize width() const {
        return width_;
    }

    output_format& flags(std::ios_base::fmtflags flags, std::ios_base::fmtflags mask) {
        this->mask_ |= mask;
        flags_ &= ~mask;
        flags_ |= flags;
        return *this;
    }

    std::ios_base::fmtflags flags() const {
        return flags_;
    }

    std::ios_base::fmtflags mask() const {
        return mask_;
    }
};

inline output_format fixed_format(std::streamsize precision, std::streamsize width) {
    output_format fmt { precision, width };
    fmt.flags(std::ios_base::fixed, std::ios_base::floatfield);
    return fmt;
}

}
}

#endif /* ADS_OUTPUT_OUTPUT_FORMAT_HPP_ */
