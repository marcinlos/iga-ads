// SPDX-FileCopyrightText: 2015 - 2025 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "poisson.hpp"

int main() {
    ads::dim_config dim{2, 40};
    ads::timesteps_config steps{1, 0};
    int ders = 1;

    auto conf = ads::config_2d{dim, dim, steps, ders};
    auto sim = ads::problems::poisson_2d{conf};
    sim.run();
}
