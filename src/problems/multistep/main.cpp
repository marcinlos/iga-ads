#include "problems/multistep/multistep2d.hpp"
#include "problems/multistep/multistep3d.hpp"
#include "problems/multistep/scheme.hpp"
#include <iostream>

using namespace ads;
using namespace ads::problems;


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

    dim_config dim{ p, n };
    timesteps_config steps{ nsteps, dt };
    int ders = 1;

    auto scm = ads::get_scheme(scheme_name);
    // std::cout << "Scheme: " << scm << std::endl;

    if (D == 2) {
        config_2d c{dim, dim, steps, ders};
        multistep2d sim{c, scm, order};
        sim.run();
    } else if (D == 3) {
        config_3d c{dim, dim, dim, steps, ders};
        multistep3d sim{c, scm, order};
        sim.run();
    }
}
