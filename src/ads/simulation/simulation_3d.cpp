#include "ads/simulation/simulation_3d.hpp"

namespace ads {

    simulation_3d::simulation_3d(const config_3d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , z{config.z, config.derivatives}
    , buffer{shape()}
    { }

}
