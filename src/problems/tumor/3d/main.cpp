#include "ads/simulation.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/3d/tumor_3d.hpp"


using namespace ads;

int main() {
    int n = 20;
    dim_config dim { 2, n, 0, 5000.0 };
    dim_config dimz { 2, n, 0, 3000.0 };

    int ders = 1;

    timesteps_config steps { 2000, 0.1 }; // 200h
    config_3d c { dim, dim, dimz, steps, ders };

    tumor::params p;

    tumor::tumor_3d sim { c, p };
    sim.run();
}
