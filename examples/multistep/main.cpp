// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <exception>
#include <iostream>

#include "multistep2d.hpp"
#include "multistep3d.hpp"
#include "scheme.hpp"

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: multistep <dim> <p> <n> <scheme> <order> <steps> <dt>" << std::endl;
        std::cerr << "Scheme format: \"s | a(s-1) ... a(0) | b(s) b(s-1) ... b(0) \"" << std::endl;
        return 0;
    }
    auto D = std::atoi(argv[1]);
    auto p = std::atoi(argv[2]);
    auto n = std::atoi(argv[3]);
    auto scheme_name = std::string{argv[4]};
    auto order = std::atoi(argv[5]);
    auto nsteps = std::atoi(argv[6]);
    auto dt = std::atof(argv[7]);

    ads::dim_config dim{p, n};
    ads::timesteps_config steps{nsteps, dt};
    int ders = 1;

    try {
        auto scm = ads::get_scheme(scheme_name);
        // std::cout << "Scheme: " << scm << std::endl;

        if (D == 2) {
            ads::config_2d c{dim, dim, steps, ders};
            ads::problems::multistep2d sim{c, scm, order};
            sim.run();
        } else if (D == 3) {
            ads::config_3d c{dim, dim, dim, steps, ders};
            ads::problems::multistep3d sim{c, scm, order};
            sim.run();
        }
    } catch (std::exception const& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }
}
