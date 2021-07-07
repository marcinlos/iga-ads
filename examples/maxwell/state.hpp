// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_STATE_HPP
#define MAXWELL_STATE_HPP

#include "ads/lin/tensor.hpp"

namespace ads {

struct state {
    using field = ads::lin::tensor<double, 3>;

    field E1;
    field E2;
    field E3;

    field H1;
    field H2;
    field H3;

    state(const std::array<std::size_t, 3> shape)
    : E1{shape}
    , E2{shape}
    , E3{shape}
    , H1{shape}
    , H2{shape}
    , H3{shape} { }

    void clear() {
        zero(E1);
        zero(E2);
        zero(E3);

        zero(H1);
        zero(H2);
        zero(H3);
    }
};

}  // namespace ads

#endif  // MAXWELL_STATE_HPP
