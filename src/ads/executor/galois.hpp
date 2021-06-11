#ifndef ADS_EXECUTOR_GALOIS_HPP_
#define ADS_EXECUTOR_GALOIS_HPP_

#include <iterator>
#include <galois/Galois.h>
#include <galois/substrate/SimpleLock.h>

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
    void for_each(Range range, Fun fun) const {
        using std::begin;
        using std::end;

        auto a = begin(range);
        auto b = end(range);
        auto items = galois::iterate(a, b);

        galois::do_all(items, [&fun](auto&& item) {
            fun(item);
        }, galois::no_stats());
    }

    explicit galois_executor(int threads);

    void thread_count(int threads);
};

}

#endif /* ADS_EXECUTOR_GALOIS_HPP_*/
