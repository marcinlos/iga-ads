#ifndef ADS_PROBLEMS_TUMOR_STATE_HPP_
#define ADS_PROBLEMS_TUMOR_STATE_HPP_

#include "ads/lin/tensor.hpp"


namespace ads {
namespace tumor {


struct state {
    static constexpr std::size_t Dim = 2;

    using field = lin::tensor<double, Dim>;

    field b;
    field c;

    field n, f, m;

    field M, A;

    state(const std::array<std::size_t, Dim>& shape)
    : b{ shape }
    , c{ shape }
    , n{ shape }, f{ shape }, m{ shape }
    , M{ shape }, A{ shape }
    { }

    void clear() {
        for (field* x : { &b, &c, &n, &f, &m, &M, &A }) {
            zero(*x);
        }
    }
};


}
}



#endif /* ADS_PROBLEMS_TUMOR_STATE_HPP_ */
