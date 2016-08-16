#include <iostream>

#include "ads/problems/heat/heat_2d.hpp"

using namespace ads;
using namespace ads::problems;


int main() {
    dim_config dim{ 2, 12 };
    timesteps_config steps{ 5000, 1e-5 };
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    heat_2d sim{c};
    sim.run();
}
