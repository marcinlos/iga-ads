#include "problems/validation/validation.hpp"

using namespace ads;
using namespace ads::problems;


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: validation <n>" << std::endl;
        return 0;
    }
    int n = std::atoi(argv[1]);
    if (n <= 0) {
        std::cerr << "Invalid value of n: " << argv[1] << std::endl;
    }

    dim_config dim{ 2, n };
    timesteps_config steps{ 100000, 1e-6 };
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    validation sim{c};
    sim.run();
}
