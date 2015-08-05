#ifndef ADS_LIN_TENSOR_TENSOR_HPP_
#define ADS_LIN_TENSOR_TENSOR_HPP_

#include "ads/lin/tensor/base.hpp"


namespace ads {
namespace lin {

namespace impl {

template <std::size_t N>
inline std::size_t product(const std::array<std::size_t, N>& ns) {
    std::size_t p = 1;
    for (auto n: ns) {
        p *= n;
    }
    return p;
}

}

template <typename T, std::size_t Rank>
struct tensor : tensor_base<T, Rank, tensor<T, Rank>> {
private:

    using Self = tensor<T, Rank>;
    using Base = tensor_base<T, Rank, Self>;
    using size_array = typename Base::size_array;

    std::vector<T> buffer_;

public:
    tensor(const size_array& sizes)
    : Base { sizes }
    , buffer_(impl::product(sizes))
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
