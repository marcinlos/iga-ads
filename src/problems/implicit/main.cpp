#include "problems/implicit/implicit.hpp"

using namespace ads;


int main() {
    dim_config dim{ 2, 80 };
    timesteps_config steps{ 100, 1e-1 };
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    implicit_2d sim{c};
    sim.run();
}
