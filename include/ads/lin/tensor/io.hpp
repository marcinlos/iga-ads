// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_IO_HPP
#define ADS_LIN_TENSOR_IO_HPP

#include <cstddef>
#include <iomanip>
#include <ostream>

#include "ads/lin/tensor/base.hpp"


namespace ads::lin {

namespace impl {

template <typename T, std::size_t Rank>
struct tensor_printer {
    static_assert(Rank <= 3, "Too high tensor rank");
};


template <typename T>
struct tensor_printer<T, 1> {

    template <typename Impl>
    static void print(std::ostream& os, const tensor_base<T, 1, Impl>& t) {
        for (std::size_t i = 0; i < t.size(0); ++ i) {
            os << std::setw(12) << t(i) << ' ';
        }
    }
};


template <typename T>
struct tensor_printer<T, 2> {

    template <typename Impl>
    static void print(std::ostream& os, const tensor_base<T, 2, Impl>& t) {
        for (std::size_t i = 0; i < t.size(0); ++ i) {
            for (std::size_t j = 0; j < t.size(1); ++ j) {
                os << std::setw(12) << t(i, j) << ' ';
            }
            os << std::endl;
        }
    }
};

template <typename T>
struct tensor_printer<T, 3> {

    template <typename Impl>
    static void print(std::ostream& os, const tensor_base<T, 3, Impl>& t) {
        for (std::size_t i = 0; i < t.size(0); ++ i) {
            for (std::size_t j = 0; j < t.size(1); ++ j) {
                for (std::size_t k = 0; k < t.size(2); ++ k) {
                    os << std::setw(12) << t(i, j, k) << ' ';
                }
                os << std::endl;
            }
            if (i < t.size(0) - 1) {
                os << std::endl;
            }
        }
    }
};

}


template <typename T, std::size_t Rank, typename Impl>
std::ostream& operator <<(std::ostream& os, const tensor_base<T, Rank, Impl>& t) {
    impl::tensor_printer<T, Rank>::print(os, t);
    return os;
}

}

#endif // ADS_LIN_TENSOR_IO_HPP
