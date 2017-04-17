#include "ads/simulation.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/3d/tumor_3d.hpp"


#include "problems/tumor/3d/vasculature_parser.hpp"

using namespace ads;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: tumor_3d threads n steps" << std::endl;
        std::exit(1);
    }
    int threads = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    int nsteps = std::atoi(argv[3]);

    // auto vessels = tumor::parse_vessels(std::cin);
    auto vessels = tumor::vessels{};

    constexpr int vasc_size = 100;
    auto vasc = tumor::vasculature{ vasc_size, vasc_size, vasc_size, std::move(vessels) };

    dim_config dim { 2, n, 0, 5000.0 };
    dim_config dimz { 2, n, 0, 3000.0 };

    int ders = 1;

    timesteps_config steps { nsteps, 0.1 }; // 200h
    config_3d c { dim, dim, dimz, steps, ders };

    tumor::params p;

    tumor::tumor_3d sim { c, p, std::move(vasc), threads };
    sim.run();

}
