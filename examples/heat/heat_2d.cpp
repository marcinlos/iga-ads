// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_2d.hpp"

using namespace ads;
using namespace ads::problems;

int main() {
    dim_config dim{2, 40};
    timesteps_config steps{10000, 1e-5};
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    heat_2d sim{c};
    sim.run();
}
