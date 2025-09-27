// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "fire.hpp"

int main() {
    ads::dim_config dim{2, 200, 0.0, 100.0};
    ads::timesteps_config steps{10000, 1e-3};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::fire sim{c};
    sim.run();
}
