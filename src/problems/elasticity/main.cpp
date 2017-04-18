#include "ads/simulation.hpp"
#include "ads/executor/sequential.hpp"
#include "ads/executor/galois.hpp"
#include "problems/elasticity/elasticity.hpp"
#include "problems/elasticity/implicit.hpp"



using namespace ads;

int main() {
    dim_config dim { 2, 12 };
    int ders = 1;

    timesteps_config steps { 4000, 1e-3 };
    config_3d c { dim, dim, dim, steps, ders };

    // problems::implicit_elasticity sim{c};
    problems::linear_elasticity sim{c};

    sim.run();
}
