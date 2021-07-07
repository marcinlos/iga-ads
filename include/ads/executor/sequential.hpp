// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_EXECUTOR_SEQUENTIAL_HPP
#define ADS_EXECUTOR_SEQUENTIAL_HPP

#include <iterator>
#include <utility>

namespace ads {

class sequential_executor {
public:
    template <typename Fun>
    void synchronized(Fun fun) const {
        fun();
    }

    template <typename Range, typename Fun>
    void for_each(Range range, Fun&& fun) const {
        using std::begin;
        using std::end;

        std::for_each(begin(range), end(range), std::forward<Fun>(fun));
    }
};

}  // namespace ads

#endif  // ADS_EXECUTOR_SEQUENTIAL_HPP
