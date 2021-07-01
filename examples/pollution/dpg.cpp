// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>

#include "pollution_dpg_v2_2d.hpp"


using namespace ads;


int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: pollution_dpg <N> <p> <k> <sep>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int p = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    int sep = std::atoi(argv[4]);

    // dim_config dim{ p, n, 0, 5000, p + k + 1 };
    // dim_config dim{ p, n, 0, 5000, p + 1, k }; // obie przestrzenie z separatorami
    // dim_config dim{ p, n, 0, 1 }; // enrichment separatorami
    dim_config dim{ p, n, 0.0, 1.0, sep};


    timesteps_config steps{ 3000, 0.5*1e-2 };

    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    pollution_dpg_v2_2d sim{c, k};
    // pollution_dpg_2d sim{c, k};

    sim.run();
}
