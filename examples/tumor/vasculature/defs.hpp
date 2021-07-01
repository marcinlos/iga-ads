// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_

#include <cstddef>

#include "ads/lin/tensor.hpp"
#include "ads/util/math/vec.hpp"


namespace tumor::vasc {

constexpr std::size_t Dim = 2;

using vector = ads::math::vec<Dim>;
using val_array = ads::lin::tensor<double, Dim>;

}


#endif /* PROBLEMS_TUMOR_VASCULATURE_DEFS_HPP_ */
