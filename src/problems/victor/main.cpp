#include <cstdlib>
#include "problems/victor/victor.hpp"


using namespace ads;

bspline::basis create_basis(double a, double b, int p, int elements, int repeated_nodes) {
    int points = elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    bspline::knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }

    auto x0 = 0.5;
    auto y0 = 0.9;

    for (int i = 1; i < points - 1; ++i) {
        auto t = lerp(i, elements, 0.0, 1.0);

        auto s = t < x0 ? t / x0 * y0 : (t - x0) / (1 - x0) * (1 - y0) + y0;
        for (int j = 0; j < r; ++ j) {
            knot[p + 1 + (i - 1) * r + j] = lerp(s, a, b);
        }
    }
    return {std::move(knot), p};
}


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
    timesteps_config steps{ nsteps, 0.5*1e-2 };
    int ders = 1;

    auto trial_basis = create_basis(0, 1, p_trial, n, p_trial - 1 - C_trial);
    auto dtrial = dimension{ trial_basis, quad, ders };

    auto test_basis = create_basis(0, 1, p_test, n, p_test - 1 - C_test);
    auto dtest = dimension{ test_basis, quad, ders };

    auto trial_dim = dtrial.B.dofs();
    auto test_dim = dtest.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    }

    victor sim{dtrial, dtrial, dtest, dtest, steps};
    sim.run();
}
