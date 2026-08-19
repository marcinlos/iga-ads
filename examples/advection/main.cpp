// SPDX-FileCopyrightText: 2026 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "advection.hpp"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: advection <p> <n> <nsteps>" << std::endl;
        return 0;
    }
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    int nsteps = std::atoi(argv[3]);

    if (n <= 0) {
        std::cerr << "Invalid value of n: " << argv[1] << std::endl;
    }

    double T = 1.0;
    double dt = T / nsteps;

    ads::dim_config dim{p, n};
    ads::timesteps_config steps{nsteps, dt};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::advection sim{c};
    sim.run();
}
