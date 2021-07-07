// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "flow.hpp"

using namespace ads;
using namespace ads::problems;

int main() {
    dim_config dim{2, 20};
    timesteps_config steps{10000, 1e-7};
    int ders = 1;

    config_3d c{dim, dim, dim, steps, ders};
    flow sim{c};
    sim.run();
}
