// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_TENSOR_DEFS_HPP
#define ADS_LIN_TENSOR_DEFS_HPP

#include "ads/lin/tensor/tensor.hpp"


namespace ads::lin {

using vector = tensor<double, 1>;
using matrix = tensor<double, 2>;

}

#endif // ADS_LIN_TENSOR_DEFS_HPP
