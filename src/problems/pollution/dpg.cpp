#include "problems/pollution/pollution_dpg_2d.hpp"

using namespace ads;


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: pollution_dpg <N> <p> <k>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int p = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);

    dim_config dim{ p, n, 0, 5000, p + k + 1 };
    timesteps_config steps{ 100, 10 };

    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    pollution_dpg_2d sim{c, k};
    sim.run();
}
