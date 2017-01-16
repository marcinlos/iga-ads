#include "ads/simulation.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/3d/tumor_3d.hpp"


using namespace ads;

int main() {
    dim_config dim { 2, 30, 0, 3000.0 };
    int ders = 1;

    timesteps_config steps { 1000, 0.001 };
    config_3d c { dim, dim, dim, steps, ders };

    tumor::params p;

    tumor::tumor_3d sim { c, p };
    sim.run();
}
