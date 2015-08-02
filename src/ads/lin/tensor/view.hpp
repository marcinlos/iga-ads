#ifndef ADS_LIN_TENSOR_VIEW_HPP_
#define ADS_LIN_TENSOR_VIEW_HPP_

#include "ads/lin/tensor/base.hpp"


namespace ads {
namespace lin {

template <typename T, std::size_t Rank>
struct tensor_view : tensor_base<T, Rank, tensor_view<T, Rank>> {
private:

    using Self = tensor_view<T, Rank>;
    using Base = tensor_base<T, Rank, Self>;

    T* data_;

public:
    template <typename... Sizes>
    tensor_view(T* data, Sizes... sizes)
    : Base { sizes... }
    , data_ { data }
    { }

    T* data() {
        return data_;
    }

    const T* data() const {
        return data_;
    }
};


}
}


#endif /* ADS_LIN_TENSOR_VIEW_HPP_ */
