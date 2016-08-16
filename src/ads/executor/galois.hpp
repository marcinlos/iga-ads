#ifndef ADS_EXECUTOR_GALOIS_HPP_
#define ADS_EXECUTOR_GALOIS_HPP_

#include <iterator>
#include <Galois/Galois.h>
#include <Galois/Threads.h>
#include <Galois/Runtime/ll/SimpleLock.h>

namespace ads {

class galois_executor {
private:
    using lock_type = Galois::Runtime::LL::SimpleLock<true>;
    lock_type lock;

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

        Galois::for_each(a, b, [&fun](auto&& item, auto&& /*ctx*/) {
            fun(item);
        });
    }

    galois_executor(int threads);

    void thread_count(int threads);
};

}

#endif /* ADS_EXECUTOR_GALOIS_HPP_*/
