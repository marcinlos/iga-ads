#include "problems/stokes/stokes.hpp"

using namespace ads;
using namespace ads::problems;

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: stokes <p> <n> <steps> <dt>" << std::endl;
        std::exit(-1);
    }
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    int nsteps = std::atoi(argv[3]);
    double dt = std::atof(argv[4]);

    dim_config dim{ p, n };
    timesteps_config steps{ nsteps, dt };

    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    stokes sim{c};
    sim.run();
}
