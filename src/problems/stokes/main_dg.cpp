#include <cstdlib>
#include <iostream>
#include <string>
#include "problems/stokes/stokes_dg.hpp"


using namespace ads;


int main(int argc, char* argv[]) {
    if (argc != 26) {
        std::cerr << "Usage: stokes <N> <p_trial> <C_trial> <p_test> <C_test>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int subdivision = 1;

    int idx = 2;
    int vxp_trial_x = std::atoi(argv[idx ++]);
    int vxc_trial_x = std::atoi(argv[idx ++]);
    int vxp_trial_y = std::atoi(argv[idx ++]);
    int vxc_trial_y = std::atoi(argv[idx ++]);
    int vyp_trial_x = std::atoi(argv[idx ++]);
    int vyc_trial_x = std::atoi(argv[idx ++]);
    int vyp_trial_y = std::atoi(argv[idx ++]);
    int vyc_trial_y = std::atoi(argv[idx ++]);
    int pp_trial_x = std::atoi(argv[idx ++]);
    int pc_trial_x = std::atoi(argv[idx ++]);
    int pp_trial_y = std::atoi(argv[idx ++]);
    int pc_trial_y = std::atoi(argv[idx ++]);

    int vxp_test_x = std::atoi(argv[idx ++]);
    int vxc_test_x = std::atoi(argv[idx ++]);
    int vxp_test_y = std::atoi(argv[idx ++]);
    int vxc_test_y = std::atoi(argv[idx ++]);
    int vyp_test_x = std::atoi(argv[idx ++]);
    int vyc_test_x = std::atoi(argv[idx ++]);
    int vyp_test_y = std::atoi(argv[idx ++]);
    int vyc_test_y = std::atoi(argv[idx ++]);
    int pp_test_x = std::atoi(argv[idx ++]);
    int pc_test_x = std::atoi(argv[idx ++]);
    int pp_test_y = std::atoi(argv[idx ++]);
    int pc_test_y = std::atoi(argv[idx ++]);

    int pmax_trial_x = std::max(std::max(vxp_trial_x, vyp_trial_x), pp_trial_x);
    int pmax_trial_y = std::max(std::max(vxp_trial_y, vyp_trial_y), pp_trial_y);
    int pmax_test_x = std::max(std::max(vxp_test_x, vyp_test_x), pp_test_x);
    int pmax_test_y = std::max(std::max(vxp_test_y, vyp_test_y), pp_test_y);

    int p_max = std::max(std::max(pmax_trial_x, pmax_trial_y), std::max(pmax_test_x, pmax_test_y));

    // int quad = p_max + 1; // to integrate velocity
    int quad = 7; // due to exponential

    timesteps_config steps{ 1, 0 };
    int ders = 2;
    int vx_rep_trial_x = vxp_trial_x - 1 - vxc_trial_x;
    int vx_rep_trial_y = vxp_trial_y - 1 - vxc_trial_y;
    int vy_rep_trial_x = vyp_trial_x - 1 - vyc_trial_x;
    int vy_rep_trial_y = vyp_trial_y - 1 - vyc_trial_y;
    int p_rep_trial_x = pp_trial_x - 1 - pc_trial_x;
    int p_rep_trial_y = pp_trial_y - 1 - pc_trial_y;

    int vx_rep_test_x = vxp_test_x - 1 - vxc_test_x;
    int vx_rep_test_y = vxp_test_y - 1 - vxc_test_y;
    int vy_rep_test_x = vyp_test_x - 1 - vyc_test_x;
    int vy_rep_test_y = vyp_test_y - 1 - vyc_test_y;
    int p_rep_test_x = pp_test_x - 1 - pc_test_x;
    int p_rep_test_y = pp_test_y - 1 - pc_test_y;


    // Pressure spaces
    auto trial_basis_x = bspline::create_basis(0, 1, pp_trial_x, n, p_rep_trial_x);
    auto dtrial_x = dimension{ trial_basis_x, quad, ders, subdivision };

    auto trial_basis_y = bspline::create_basis(0, 1, pp_trial_y, n, p_rep_trial_y);
    auto dtrial_y = dimension{ trial_basis_y, quad, ders, subdivision };

    auto test_basis_x = bspline::create_basis(0, 1, pp_test_x, subdivision*n, p_rep_test_x);
    auto dtest_x = dimension{ test_basis_x, quad, ders, 1 };

    auto test_basis_y = bspline::create_basis(0, 1, pp_test_y, subdivision*n, p_rep_test_y);
    auto dtest_y = dimension{ test_basis_y, quad, ders, 1 };

    // Velocity spaces
    auto U1_trial_basis_x = bspline::create_basis(0, 1, vxp_trial_x, n, vx_rep_trial_x);
    auto U1_dtrial_x = dimension{ U1_trial_basis_x, quad, ders, subdivision };

    auto U1_trial_basis_y = bspline::create_basis(0, 1, vxp_trial_y, n, vx_rep_trial_y);
    auto U1_dtrial_y = dimension{ U1_trial_basis_y, quad, ders, subdivision };

    auto U1_test_basis_x = bspline::create_basis(0, 1, vxp_test_x, subdivision*n, vx_rep_test_x);
    auto U1_dtest_x = dimension{ U1_test_basis_x, quad, ders, 1 };

    auto U1_test_basis_y = bspline::create_basis(0, 1, vxp_test_y, subdivision*n, vx_rep_test_y);
    auto U1_dtest_y = dimension{ U1_test_basis_y, quad, ders, 1 };


    auto U2_trial_basis_x = bspline::create_basis(0, 1, vyp_trial_x, n, vy_rep_trial_x);
    auto U2_dtrial_x = dimension{ U2_trial_basis_x, quad, ders, subdivision };

    auto U2_trial_basis_y = bspline::create_basis(0, 1, vyp_trial_y, n, vy_rep_trial_y);
    auto U2_dtrial_y = dimension{ U2_trial_basis_y, quad, ders, subdivision };

    auto U2_test_basis_x = bspline::create_basis(0, 1, vyp_test_x, subdivision*n, vy_rep_test_x);
    auto U2_dtest_x = dimension{ U2_test_basis_x, quad, ders, 1 };

    auto U2_test_basis_y = bspline::create_basis(0, 1, vyp_test_y, subdivision*n, vy_rep_test_y);
    auto U2_dtest_y = dimension{ U2_test_basis_y, quad, ders, 1 };


    // Sanity check
    auto trial_v_dim = U1_dtrial_x.dofs() * U1_dtrial_y.dofs() + U2_dtrial_x.dofs() * U2_dtrial_y.dofs();
    auto trial_p_dim = dtrial_x.dofs() * dtrial_y.dofs();
    auto trial_dim = trial_v_dim + trial_p_dim;

    auto test_v_dim = U1_dtest_x.dofs() * U1_dtest_y.dofs() + U2_dtest_x.dofs() * U2_dtest_y.dofs();
    auto test_p_dim = dtest_x.dofs() * dtest_y.dofs();
    auto test_dim = test_v_dim + test_p_dim;

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::cerr << "dim Ux = " << U1_dtrial_x.dofs() * U1_dtrial_y.dofs()
                  << " = " << U1_dtrial_x.dofs() << " x " << U1_dtrial_y.dofs() << std::endl;
        std::cerr << "dim Uy = " << U2_dtrial_x.dofs() * U2_dtrial_y.dofs()
                  << " = " << U2_dtrial_x.dofs() << " x " << U2_dtrial_y.dofs() << std::endl;
        std::cerr << "dim P =  " << dtrial_x.dofs() * dtrial_y.dofs()
                  << " = " << dtrial_x.dofs() << " x " << dtrial_y.dofs() << std::endl;

        std::cerr << "dim Vx = " << U1_dtest_x.dofs() * U1_dtest_y.dofs()
                  << " = " << U1_dtest_x.dofs() << " x " << U1_dtest_y.dofs() << std::endl;
        std::cerr << "dim Vy = " << U2_dtest_x.dofs() * U2_dtest_y.dofs()
                  << " = " << U2_dtest_x.dofs() << " x " << U2_dtest_y.dofs() << std::endl;
        std::cerr << "dim Q =  " << dtest_x.dofs() * dtest_y.dofs()
                  << " = " << dtest_x.dofs() << " x " << dtest_y.dofs() << std::endl;
        std::exit(1);
    } else {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
        std::cout << "dofs = " << trial_dim + test_dim << std::endl;
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

    auto sim = stokes_dg{trial, test, steps};
    sim.run();
}
