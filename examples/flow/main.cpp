// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "flow.hpp"

int main() {
    ads::dim_config dim{2, 20};
    ads::timesteps_config steps{10000, 1e-7};
    int ders = 1;

    ads::config_3d c{dim, dim, dim, steps, ders};
    ads::problems::flow sim{c};
    sim.run();
}
