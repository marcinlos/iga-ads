#ifndef ADS_UTIL_META_HPP_
#define ADS_UTIL_META_HPP_

#include <type_traits>


namespace ads::util {

template <typename...>
struct and_;

template <>
struct and_<> : std::true_type
{ };

template <typename T>
struct and_<T> : T
{ };

template <typename T, typename... Ts>
struct and_<T, Ts...> : std::conditional<T::value, and_<Ts...>, T>::type
{ };

template <template <typename...> class Pred, typename... Types>
struct all_ : and_<Pred<Types>...>
{ };

}

#endif /* ADS_UTIL_META_HPP_ */
