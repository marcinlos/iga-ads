#include <cstdlib>
#include "problems/victor/victor.hpp"


using namespace ads;


int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Usage: victor <N> <p_trial> <C_trial> <p_test> <C_test> <steps>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int p_trial = std::atoi(argv[2]);
    int C_trial = std::atoi(argv[3]);
    int p_test = std::atoi(argv[4]);
    int C_test = std::atoi(argv[5]);
    int nsteps = std::atoi(argv[6]);

    int quad = std::max(p_trial, p_test) + 1;
    dim_config trial{ p_trial, n, 0.0, 1.0, quad, p_trial - 1 - C_trial};
    dim_config test { p_test,  n, 0.0, 1.0, quad, p_test  - 1 - C_test };

    auto trial_dim = dimension(trial, 0).B.dofs();
    auto test_dim = dimension(test, 0).B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    }

    timesteps_config steps{ nsteps, 0.5*1e-2 };

    int ders = 1;
    config_2d c{trial, trial, steps, ders};
    victor sim{c, test, test};

    sim.run();
}
