#include "problems/pollution/pollution_2d.hpp"

using namespace ads;


int main() {
    dim_config dim{ 2, 400, 0, 5000 };
    timesteps_config steps{ 10000, 10 };
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    pollution_2d sim{c};
    sim.run();
}
