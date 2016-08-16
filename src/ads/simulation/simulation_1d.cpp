#include "ads/simulation/simulation_1d.hpp"

namespace ads {

    simulation_1d::simulation_1d(const config_1d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    { }

}
