// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_STATE_HPP
#define TUMOR_STATE_HPP

#include "ads/lin/tensor.hpp"


namespace tumor {

template <std::size_t Dim>
struct state {

    using field = ads::lin::tensor<double, Dim>;

    field b;
    field c;
    field o;

    field M, A;

    state(std::array<std::size_t, Dim> shape)
    : b{ shape }
    , c{ shape }
    , o{ shape }
    , M{ shape }, A{ shape }
    { }

    void clear() {
        for (field* x : { &b, &c, &o, &M, &A }) {
            zero(*x);
        }
    }
};

}

#endif // TUMOR_STATE_HPP
