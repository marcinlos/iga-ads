#include "problems/pollution/pollution_3d.hpp"

using namespace ads;


int main() {
    dim_config dim{ 2, 100, 0, 5000 };
    timesteps_config steps{ 10000, 10 };
    int ders = 1;

    config_3d c{dim, dim, dim, steps, ders};
    pollution_3d sim{c};
    sim.run();
}
