#ifndef PROBLEMS_TUMOR_STATE_HPP_
#define PROBLEMS_TUMOR_STATE_HPP_

#include "ads/lin/tensor.hpp"


namespace tumor {

template <std::size_t Dim>
struct state {

    using field = ads::lin::tensor<double, Dim>;

    field b;
    field c;

    field M, A;

    state(std::array<std::size_t, Dim> shape)
    : b{ shape }
    , c{ shape }
    , M{ shape }, A{ shape }
    { }

    void clear() {
        for (field* x : { &b, &c, &M, &A }) {
            zero(*x);
        }
    }
};

}

#endif /* PROBLEMS_TUMOR_STATE_HPP_ */
