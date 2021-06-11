#include "coupled.hpp"

using namespace ads;


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
    dim_config dim{ p, n };

    // timesteps_config steps{ 100, 1e-2 };
    timesteps_config steps{ nsteps, dt };

    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    coupled sim{c};
    sim.run();
}
