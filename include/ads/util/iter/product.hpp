// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_ITER_PRODUCT_HPP
#define ADS_UTIL_ITER_PRODUCT_HPP

#include <iterator>
#include <tuple>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range.hpp>


namespace ads::util {

namespace impl {

template <typename Iter>
using val_type = typename std::iterator_traits<Iter>::value_type;

template <typename... Iter>
using iter_tuple = std::tuple<val_type<Iter>...>;

}

template <typename Iter, typename Out = impl::iter_tuple<Iter, Iter>>
class iter_product2
  : public boost::iterator_facade<
        iter_product2<Iter, Out>,
        Out,
        boost::forward_traversal_tag,
        Out
    >
{
private:
    Iter iter1, iter2;
    boost::iterator_range<Iter> range2;

public:
    iter_product2(Iter iter1, Iter iter2, boost::iterator_range<Iter> range2)
    : iter1{ iter1 }
    , iter2{ iter2 }
    , range2{ range2 }
    { }

private:
    friend class boost::iterator_core_access;

    void increment() {
        using boost::begin;
        using boost::end;

        ++ iter2;
        if (iter2 == end(range2)) {
            iter2 = begin(range2);
            ++ iter1;
        }
    }

    Out dereference() const {
        return Out{ *iter1, *iter2 };
    }

    bool equal(const iter_product2<Iter, Out>& other) const {
        return iter1 == other.iter1 && iter2 == other.iter2 && range2 == other.range2;
    }
};

template <typename Iter, typename Out = impl::iter_tuple<Iter, Iter>>
iter_product2<Iter, Out>
product_iter(Iter iter1, Iter iter2, boost::iterator_range<Iter> range2) {
    return { iter1, iter2, range2 };
}


template <typename Out, typename Iter>
boost::iterator_range<iter_product2<Iter, Out>>
product_range(boost::iterator_range<Iter> rx, boost::iterator_range<Iter> ry) {
    using boost::begin;
    using boost::end;

    auto it_begin = product_iter<Iter, Out>(begin(rx), begin(ry), ry);
    auto it_end   = product_iter<Iter, Out>(end(rx),   begin(ry), ry);
    return boost::make_iterator_range(it_begin, it_end);
}



template <typename Iter, typename Out = impl::iter_tuple<Iter, Iter, Iter>>
class iter_product3
  : public boost::iterator_facade<
        iter_product3<Iter, Out>,
        Out,
        boost::forward_traversal_tag,
        Out
    >
{
private:
    using iter_range = boost::iterator_range<Iter>;

    Iter iter1, iter2, iter3;
    iter_range range2, range3;

public:
    iter_product3(Iter iter1, Iter iter2, Iter iter3, iter_range range2, iter_range range3)
    : iter1{ iter1 }
    , iter2{ iter2 }
    , iter3{ iter3 }
    , range2{ range2 }
    , range3{ range3 }
    { }

private:
    friend class boost::iterator_core_access;

    void increment() {
        using boost::begin;
        using boost::end;

        ++ iter3;
        if (iter3 == end(range3)) {
            iter3 = begin(range3);
            increment_level2();
        }
    }

    void increment_level2() {
        using boost::begin;
        using boost::end;

        ++ iter2;
        if (iter2 == end(range2)) {
            iter2 = begin(range2);
            ++ iter1;
        }
    }

    Out dereference() const {
        return Out{ *iter1, *iter2, *iter3 };
    }

    bool equal(const iter_product3<Iter, Out>& other) const {
        return iter1 == other.iter1 && iter2 == other.iter2 && iter3 == other.iter3
                && range2 == other.range2 && range3 == other.range3;
    }
};

template <typename Iter, typename Out = impl::iter_tuple<Iter, Iter, Iter>>
iter_product3<Iter, Out>
product_iter(Iter iter1, Iter iter2, Iter iter3, boost::iterator_range<Iter> range2, boost::iterator_range<Iter> range3) {
    return { iter1, iter2, iter3, range2, range3 };
}

template <typename Out, typename Iter>
boost::iterator_range<iter_product3<Iter, Out>>
product_range(boost::iterator_range<Iter> rx, boost::iterator_range<Iter> ry, boost::iterator_range<Iter> rz) {
    using boost::begin;
    using boost::end;

    auto it_begin = product_iter<Iter, Out>(begin(rx), begin(ry), begin(rz), ry, rz);
    auto it_end   = product_iter<Iter, Out>(end(rx),   begin(ry), begin(rz), ry, rz);
    return boost::make_iterator_range(it_begin, it_end);
}

}

#endif // ADS_UTIL_ITER_PRODUCT_HPP
