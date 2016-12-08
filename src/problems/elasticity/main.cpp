#include "ads/simulation.hpp"
#include "ads/executor/sequential.hpp"
#include "ads/executor/galois.hpp"
#include "problems/elasticity/elasticity.hpp"


using namespace ads;

int main() {
    dim_config dim { 2, 12 };
    int ders = 1;

    timesteps_config steps { 40000, 1e-4 };
    config_3d c { dim, dim, dim, steps, ders };

    problems::linear_elasticity sim{c};
    sim.run();
}
