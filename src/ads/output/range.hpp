#ifndef ADS_OUTPUT_RANGE_HPP_
#define ADS_OUTPUT_RANGE_HPP_

#include <iterator>
#include <type_traits>

namespace ads {
namespace output {

template <typename RangeIter>
struct range {

    using iterator = RangeIter;
    using iter_traits = std::iterator_traits<iterator>;
    using value_type = typename iter_traits::value_type;
    using reference = typename iter_traits::reference;

    RangeIter begin;
    RangeIter end;

    std::size_t size() const {
        using std::distance;
        return distance(begin, end);
    }

    value_type operator [](std::size_t i) const {
        return begin[i];
    }

    reference operator [](std::size_t i) {
        return begin[i];
    }

private:
    using iterator_category = typename iter_traits::iterator_category;

    static_assert(std::is_same<iterator_category, std::random_access_iterator_tag>::value,
            "Range can be constructed only from random access iterators");
};

template <typename RangeIter>
RangeIter begin(const range<RangeIter>& r) {
    return r.begin;
}

template <typename RangeIter>
RangeIter end(const range<RangeIter>& r) {
    return r.end;
}


template <typename Cont>
auto from_container(const Cont& cont) -> range<decltype(begin(cont))> {
    using std::begin;
    using std::end;
    return range<decltype(begin(cont))> { begin(cont), end(cont) };
}


}
}

#endif /* ADS_OUTPUT_RANGE_HPP_ */
