// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "coupled.hpp"

#include <cstdlib>

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: coupled <p> <n> <steps> <dt>" << std::endl;
        std::exit(-1);
    }
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    int nsteps = std::atoi(argv[3]);
    double dt = std::atof(argv[4]);

    // dim_config dim{ 2, 80 };
    ads::dim_config dim{p, n};

    // timesteps_config steps{ 100, 1e-2 };
    ads::timesteps_config steps{nsteps, dt};

    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::coupled sim{c};
    sim.run();
}
