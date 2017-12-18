#include "ads/simulation/simulation_2d.hpp"

namespace ads {

    simulation_2d::simulation_2d(const config_2d& config)
    : simulation_2d{
        dimension{config.x, config.derivatives},
        dimension{config.y, config.derivatives},
        config.steps}
    { }

    simulation_2d::simulation_2d(dimension x, dimension y, const timesteps_config& steps)
    : simulation_base{steps}
    , x{std::move(x)}
    , y{std::move(y)}
    , buffer{shape()}
    { }


}
