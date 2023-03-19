// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>

#include "pollution_dpg_v2_2d.hpp"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: pollution_dpg <N> <p> <k> <sep>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int p = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    int sep = std::atoi(argv[4]);

    // ads::dim_config dim{ p, n, 0, 5000, p + k + 1 };
    // ads::dim_config dim{ p, n, 0, 5000, p + 1, k }; // obie przestrzenie z separatorami
    // ads::dim_config dim{ p, n, 0, 1 }; // enrichment separatorami
    ads::dim_config dim{p, n, 0.0, 1.0, sep};

    ads::timesteps_config steps{3000, 0.5 * 1e-2};

    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::pollution_dpg_v2_2d sim{c, k};
    // ads::pollution_dpg_2d sim{c, k};

    sim.run();
}
