#include "test2d.hpp"
#include "test3d.hpp"


using namespace ads;
using namespace ads::problems;


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

    dim_config dim{ p, n };
    timesteps_config steps{ ts, 1e-6 };
    int ders = 1;

    if (D == 2) {
        config_2d c{dim, dim, steps, ders};
        scalability_2d sim{c, threads};
        sim.run();
    } else if (D == 3) {
        config_3d c{dim, dim, dim, steps, ders};
        scalability_3d sim{c, threads};
        sim.run();
    } else {
        std::cerr << "Invalid dimension: " << D << std::endl;
    }
}
