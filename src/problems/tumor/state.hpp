#ifndef PROBLEMS_TUMOR_STATE_HPP_
#define PROBLEMS_TUMOR_STATE_HPP_

#include "ads/lin/tensor.hpp"


namespace ads {
namespace tumor {


template <std::size_t Dim>
struct state {

    using field = lin::tensor<double, Dim>;

    field b;
    field c;

    field n, f, m;

    field M, A;

    state(std::array<std::size_t, Dim> shape)
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



#endif /* PROBLEMS_TUMOR_STATE_HPP_ */
