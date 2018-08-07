#include <cstdlib>
#include <iostream>
#include <string>

#include "ads/bspline/bspline.hpp"
#include "problems/stokes/stokes.hpp"

using namespace ads;

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: stokes <N> <p_trial> <C_trial> <p_test> <C_test>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int subdivision = 1;

    int p_trial = std::atoi(argv[2]);
    int C_trial = std::atoi(argv[3]);
    int p_test = std::atoi(argv[4]);
    int C_test = std::atoi(argv[5]);

    int quad = std::max(p_trial, p_test) + 1;

    timesteps_config steps{ 1, 0 };
    int ders = 2;
    int rep_trial = p_trial - 1 - C_trial;
    int rep_test = p_test - 1 - C_test;

    auto trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_x = dimension{ trial_basis_x, quad, ders, subdivision };

    auto trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_y = dimension{ trial_basis_y, quad, ders, subdivision };

    auto test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_x = dimension{ test_basis_x, quad, ders, 1 };

    auto test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_y = dimension{ test_basis_y, quad, ders, 1 };

    auto trial_dim = dtrial_x.B.dofs();
    auto test_dim = dtest_x.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
    }

    auto sim = stokes{dtrial_x, dtrial_y, dtest_x, dtest_y, steps};
    sim.run();
}
