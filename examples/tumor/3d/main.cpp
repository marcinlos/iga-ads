// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>

#include "../vasculature.hpp"
#include "ads/simulation.hpp"
#include "tumor_3d.hpp"
#include "vasculature_parser.hpp"


using namespace ads;

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: tumor_3d threads p n vasc_size steps" << std::endl;
        std::exit(1);
    }
    int threads = std::atoi(argv[1]);
    int p = std::atoi(argv[2]);
    int n = std::atoi(argv[3]);
    int vasc_size = std::atoi(argv[4]); // 300; //16 * 3;
    int nsteps = std::atoi(argv[5]);

    auto vessels = tumor::parse_vessels(std::cin);
    // auto vessels = tumor::vessels{};

    auto vasc = tumor::vasculature{ vasc_size, vasc_size, vasc_size, std::move(vessels) };

    dim_config dim { p, n, 0, 5000.0 };
    dim_config dimz { p, n, 0, 5000.0 };

    int ders = 1;

    timesteps_config steps { nsteps, 0.1 }; // 200h
    config_3d c { dim, dim, dimz, steps, ders };

    tumor::params par;

    tumor::tumor_3d sim { c, par, std::move(vasc), threads };
    sim.run();

}
