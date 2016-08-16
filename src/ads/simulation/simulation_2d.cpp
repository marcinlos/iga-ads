#include "ads/simulation/simulation_2d.hpp"

namespace ads {

    simulation_2d::simulation_2d(const config_2d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , buffer{shape()}
    { }

}
