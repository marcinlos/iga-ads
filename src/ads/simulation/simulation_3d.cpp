// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/simulation/simulation_3d.hpp"

namespace ads {

simulation_3d::simulation_3d(const config_3d& config)
: simulation_3d{
    dimension{config.x, config.derivatives},
    dimension{config.y, config.derivatives},
    dimension{config.z, config.derivatives},
    config.steps,
} { }

simulation_3d::simulation_3d(const dimension& x, const dimension& y, const dimension& z,
                             const timesteps_config& steps)
: simulation_base{steps}
, x{x}
, y{y}
, z{z}
, buffer{shape()} { }

}  // namespace ads
