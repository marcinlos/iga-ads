#include "ads/simulation.hpp"
#include "ads/executor/sequential.hpp"
#include "ads/executor/galois.hpp"
#include "problems/elasticity/elasticity_victor.hpp"



using namespace ads;

int main() {
    dim_config dim { 2, 20 };
    int ders = 1;

    timesteps_config steps { 4000, 2.7e-2 };
    config_3d c { dim, dim, dim, steps, ders };

    problems::elasticity_victor sim{c};

    sim.run();
}
