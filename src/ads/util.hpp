#ifndef ADS_UTIL_HPP_
#define ADS_UTIL_HPP_


namespace ads {


template <typename Num>
Num lerp(Num t, Num a, Num b) {
    return (1 - t) * a + t * b;
}

template <typename Int1, typename Int2, typename Num>
Num lerp(Int1 i, Int2 n, Num a, Num b) {
    Num t = static_cast<Num>(i) / n;
    return lerp(t, a, b);
}

}

#endif /* ADS_UTIL_HPP_ */
