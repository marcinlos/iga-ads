#include <cstdlib>
#include <iostream>
#include <string>

#include "ads/bspline/bspline.hpp"
#include "problems/stokes/stokes.hpp"
#include "problems/stokes/stokes_conforming.hpp"
#include "problems/stokes/stokes_constrained.hpp"

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

    int quad = std::max(p_trial, p_test) + 1 + 1; // to integrate velocity

    timesteps_config steps{ 1, 0 };
    int ders = 2;
    int rep_trial = p_trial - 1 - C_trial;
    int rep_test = p_test - 1 - C_test;

    // Pressure spaces
    auto trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_x = dimension{ trial_basis_x, quad, ders, subdivision };

    auto trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_y = dimension{ trial_basis_y, quad, ders, subdivision };

    auto test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_x = dimension{ test_basis_x, quad, ders, 1 };

    auto test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_y = dimension{ test_basis_y, quad, ders, 1 };
    // Velocity spaces
    auto U1_trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U1_dtrial_x = dimension{ U1_trial_basis_x, quad, ders, subdivision };

    auto U1_trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U1_dtrial_y = dimension{ U1_trial_basis_y, quad, ders, subdivision };

    auto U1_test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U1_dtest_x = dimension{ U1_test_basis_x, quad, ders, 1 };

    auto U1_test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U1_dtest_y = dimension{ U1_test_basis_y, quad, ders, 1 };


    auto U2_trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U2_dtrial_x = dimension{ U2_trial_basis_x, quad, ders, subdivision };

    auto U2_trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U2_dtrial_y = dimension{ U2_trial_basis_y, quad, ders, subdivision };

    auto U2_test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U2_dtest_x = dimension{ U2_test_basis_x, quad, ders, 1 };

    auto U2_test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U2_dtest_y = dimension{ U2_test_basis_y, quad, ders, 1 };



    // Sanity check
    auto trial_dim = dtrial_x.B.dofs();
    auto test_dim = dtest_x.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
        auto dofs = trial_dim * trial_dim + test_dim * test_dim;
        std::cout << "dofs = " << dofs << std::endl;
    }

    auto trial = space_set{
        U1_dtrial_x, U1_dtrial_y,
        U2_dtrial_x, U2_dtrial_y,
        dtrial_x, dtrial_y
    };

    auto test = space_set{
        U1_dtest_x, U1_dtest_y,
        U2_dtest_x, U2_dtest_y,
        dtest_x, dtest_y
    };

    // auto sim = stokes_conforming{trial, test, steps};
    // auto sim = stokes_constrained{trial, test, steps};
    auto sim = stokes{dtrial_x, dtrial_y, dtest_x, dtest_y, steps};

    sim.run();
}
