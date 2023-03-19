// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_2d.hpp"

int main() {
    ads::dim_config dim{2, 40};
    ads::timesteps_config steps{10000, 1e-5};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::heat_2d sim{c};
    sim.run();
}
