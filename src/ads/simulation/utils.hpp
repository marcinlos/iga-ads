#ifndef ADS_SIMULATION_UTILS_HPP_
#define ADS_SIMULATION_UTILS_HPP_

#include "ads/simulation/dimension.hpp"

namespace ads {

    inline double min_element_size(const dimension& U) {
        return 2 * *std::min_element(U.basis.J, U.basis.J + U.elements);
    }

    inline double max_element_size(const dimension& U) {
        return 2 * *std::max_element(U.basis.J, U.basis.J + U.elements);
    }

}


#endif // ADS_SIMULATION_UTILS_HPP_
