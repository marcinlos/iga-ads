#include "ads/simulation.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/tumor.hpp"


using namespace ads;

int main() {
    dim_config dim { 2, 80, 0, 3000.0 };
    int ders = 1;

    timesteps_config steps { 10000, 0.1 };
    config_2d c { dim, dim, steps, ders };

    tumor::params p;
    tumor::vasc::random_vasculature rand_vasc { 0 };

    tumor::tumor_2d sim { c, p, rand_vasc() };
    sim.run();
}
