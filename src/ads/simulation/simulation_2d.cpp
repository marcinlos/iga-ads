// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/simulation/simulation_2d.hpp"

namespace ads {

simulation_2d::simulation_2d(const config_2d& config)
: simulation_2d{
    dimension{config.x, config.derivatives},
    dimension{config.y, config.derivatives},
    config.steps,
} { }

simulation_2d::simulation_2d(const dimension& x, const dimension& y, const timesteps_config& steps)
: simulation_base{steps}
, x{x}
, y{y}
, buffer{shape()} { }

}  // namespace ads
