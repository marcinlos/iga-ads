// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_VASCULATURE_DEFS_HPP
#define TUMOR_VASCULATURE_DEFS_HPP

#include <cstddef>

#include "ads/lin/tensor.hpp"
#include "ads/util/math/vec.hpp"


namespace tumor::vasc {

constexpr std::size_t Dim = 2;

using vector = ads::math::vec<Dim>;
using val_array = ads::lin::tensor<double, Dim>;

}


#endif // TUMOR_VASCULATURE_DEFS_HPP
