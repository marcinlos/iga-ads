// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/simulation/simulation_1d.hpp"

namespace ads {

simulation_1d::simulation_1d(const config_1d& config)
: simulation_1d{dimension{config.x, config.derivatives}, config.steps} { }

simulation_1d::simulation_1d(const dimension& x, const timesteps_config& steps)
: simulation_base{steps}
, x{x} { }

}  // namespace ads
