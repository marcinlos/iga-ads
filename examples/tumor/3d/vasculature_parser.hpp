// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_3D_VASCULATURE_PARSER_HPP
#define TUMOR_3D_VASCULATURE_PARSER_HPP

#include <iostream>
#include <vector>

#include "vasculature.hpp"

namespace tumor {

void normalize_positions(std::vector<vessels::point_type>& points);

vessels parse_vessels(std::istream& is);

}  // namespace tumor

#endif  // TUMOR_3D_VASCULATURE_PARSER_HPP
