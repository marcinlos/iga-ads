#ifndef ADS_EXECUTOR_SEQUENTIAL_HPP_
#define ADS_EXECUTOR_SEQUENTIAL_HPP_

#include <algorithm>


namespace ads {

class sequential_executor {
public:
    template <typename Fun>
    void synchronized(Fun fun) const {
        fun();
    }

    template <typename Range, typename Fun>
    void for_each(Range range, Fun fun) const {
        using std::begin;
        using std::end;

        std::for_each(begin(range), end(range), [&fun](auto&& item) {
            fun(item);
        });
    }
};

}

#endif /* ADS_EXECUTOR_SEQUENTIAL_HPP_ */
