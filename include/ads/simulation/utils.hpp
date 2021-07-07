// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SIMULATION_UTILS_HPP
#define ADS_SIMULATION_UTILS_HPP

#include "ads/simulation/dimension.hpp"

namespace ads {

inline double min_element_size(const dimension& U) {
    return 2 * *std::min_element(U.basis.J, U.basis.J + U.elements);
}

inline double max_element_size(const dimension& U) {
    return 2 * *std::max_element(U.basis.J, U.basis.J + U.elements);
}

}  // namespace ads

#endif  // ADS_SIMULATION_UTILS_HPP
