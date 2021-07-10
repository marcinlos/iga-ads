// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_EXECUTOR_GALOIS_HPP
#define ADS_EXECUTOR_GALOIS_HPP

#include "ads/config.hpp"

#ifdef ADS_USE_GALOIS

#    include <iterator>
#    include <utility>

#    include <galois/Galois.h>
#    include <galois/substrate/SimpleLock.h>

namespace ads {

class galois_executor {
private:
    using lock_type = galois::substrate::SimpleLock;
    lock_type lock;

    galois::SharedMemSys G;

public:
    template <typename Fun>
    void synchronized(Fun fun) const {
        lock.lock();
        fun();
        lock.unlock();
    }

    template <typename Range, typename Fun>
    void for_each(Range range, Fun&& fun) const {
        using std::begin;
        using std::end;

        auto a = begin(range);
        auto b = end(range);
        auto items = galois::iterate(a, b);

        galois::do_all(items, std::forward<Fun>(fun), galois::no_stats{});
    }

    explicit galois_executor(int threads);

    void thread_count(int threads);
};

}  // namespace ads

#endif  // defined(ADS_USE_GALOIS)

#endif  // ADS_EXECUTOR_GALOIS_HPP
