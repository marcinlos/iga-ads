// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "heat_1d.hpp"

int main() {
    ads::dim_config dim{2, 16};
    ads::timesteps_config steps{10, 1e-5};
    int ders = 1;

    ads::config_1d c{dim, steps, ders};
    ads::problems::heat_1d sim{c};
    sim.run();
}
