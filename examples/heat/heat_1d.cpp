// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_1d.hpp"


using namespace ads;
using namespace ads::problems;

int main() {
    dim_config dim{ 2, 16 };
    timesteps_config steps{ 10, 1e-5 };
    int ders = 1;

    config_1d c{dim, steps, ders};
    heat_1d sim{c};
    sim.run();
}
