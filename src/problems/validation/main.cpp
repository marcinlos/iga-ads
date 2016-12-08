#include "problems/validation/validation.hpp"

using namespace ads;
using namespace ads::problems;


int main() {
    dim_config dim{ 2, 20 };
    timesteps_config steps{ 100000, 1e-6 };
    int ders = 1;

    config_2d c{dim, dim, steps, ders};
    validation sim{c};
    sim.run();
}
