#ifndef ADS_LIN_TENSOR_TENSOR_HPP_
#define ADS_LIN_TENSOR_TENSOR_HPP_

#include "ads/lin/tensor/base.hpp"


namespace ads {
namespace lin {

namespace impl {

constexpr inline std::size_t product() {
    return 1;
}

template <typename Value, typename... Values>
constexpr Value product(Value a, Values... as) {
    return a * product(as...);
}

}

template <typename T, std::size_t Rank>
struct tensor : tensor_base<T, Rank, tensor<T, Rank>> {
private:

    using Self = tensor<T, Rank>;
    using Base = tensor_base<T, Rank, Self>;

    std::vector<T> buffer_;

public:
    template <typename... Sizes>
    tensor(Sizes... sizes)
    : Base { sizes... }
    , buffer_(impl::product(sizes...))
    { }

    T* data() {
        return buffer_.data();
    }

    const T* data() const {
        return buffer_.data();
    }
};

}
}


#endif /* ADS_LIN_TENSOR_TENSOR_HPP_ */
