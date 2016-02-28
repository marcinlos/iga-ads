#include <iostream>

#include "ads/util/iter/product.hpp"
#include "ads/problems/heat/heat_1d.hpp"
#include "ads/problems/heat/heat_3d.hpp"
#include "ads/simulation.hpp"

using namespace ads;
using namespace ads::problems;


int main() {
//    dim_config x{ 2, 20 };
//    timesteps_config steps{ 20000, 1e-5 };
//    int ders = 1;
//
//    config_1d c{x, steps, ders};
//    heat_1d sim{c};
//    sim.run();

    dim_config dim{ 2, 12 };
    timesteps_config steps{ 100, 1e-7 };
    int ders = 1;

    config_3d c{dim, dim, dim, steps, ders};
    heat_3d sim{c};
    sim.run();
}
