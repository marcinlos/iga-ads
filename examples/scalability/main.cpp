// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "test2d.hpp"
#include "test3d.hpp"

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: scalability <dim> <threads> <p> <n> <steps>" << std::endl;
        return 0;
    }
    int D = std::atoi(argv[1]);
    int threads = std::atoi(argv[2]);
    int p = std::atoi(argv[3]);
    int n = std::atoi(argv[4]);
    int ts = std::atoi(argv[5]);

    ads::dim_config dim{p, n};
    ads::timesteps_config steps{ts, 1e-6};
    int ders = 1;

    if (D == 2) {
        ads::config_2d c{dim, dim, steps, ders};
        ads::problems::scalability_2d sim{c, threads};
        sim.run();
    } else if (D == 3) {
        ads::config_3d c{dim, dim, dim, steps, ders};
        ads::problems::scalability_3d sim{c, threads};
        sim.run();
    } else {
        std::cerr << "Invalid dimension: " << D << std::endl;
    }
}
