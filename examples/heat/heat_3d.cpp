// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_3d.hpp"

int main() {
    ads::dim_config dim{2, 12};
    ads::timesteps_config steps{5000, 1e-7};
    int ders = 1;

    ads::config_3d c{dim, dim, dim, steps, ders};
    ads::problems::heat_3d sim{c};
    sim.run();
}
