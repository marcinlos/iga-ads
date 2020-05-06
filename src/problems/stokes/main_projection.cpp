#include <cstdlib>
#include <iostream>
#include <string>

#include "ads/bspline/bspline.hpp"
#include "problems/stokes/stokes_projection.hpp"


using namespace ads;

int total_dimension(const space_set& s) {
    auto u_dim = s.U1x.dofs() * s.U1y.dofs() + s.U2x.dofs() * s.U2y.dofs();
    auto p_dim = s.Px.dofs() * s.Py.dofs();
    return u_dim + p_dim;
}

dimension make_dimension(int p, int c, int n, int quad, int ders) {
    int rep = p - 1 - c;
    auto basis = bspline::create_basis(0, 1, p, n, rep);
    return dimension{ basis, quad, ders, 1 };
}

int main(int argc, char* argv[]) {
    if (argc != 29) {
        std::cerr << "Usage: stokes_projection <N> <steps> <Navier-Stokes> <Re> <spaces...>" << std::endl;
        std::exit(1);
    }
    int idx = 1;
    int n = std::atoi(argv[idx ++]);
    int nsteps = std::atoi(argv[idx ++]);
    bool NS = std::atoi(argv[idx ++]);
    auto Re = std::atof(argv[idx ++]);

    int vxp_trial_x = std::atoi(argv[idx ++]);
    int vxc_trial_x = std::atoi(argv[idx ++]);
    int vxp_trial_y = std::atoi(argv[idx ++]);
    int vxc_trial_y = std::atoi(argv[idx ++]);
    int vyp_trial_x = std::atoi(argv[idx ++]);
    int vyc_trial_x = std::atoi(argv[idx ++]);
    int vyp_trial_y = std::atoi(argv[idx ++]);
    int vyc_trial_y = std::atoi(argv[idx ++]);
    int pp_trial_x  = std::atoi(argv[idx ++]);
    int pc_trial_x  = std::atoi(argv[idx ++]);
    int pp_trial_y  = std::atoi(argv[idx ++]);
    int pc_trial_y  = std::atoi(argv[idx ++]);

    int vxp_test_x = std::atoi(argv[idx ++]);
    int vxc_test_x = std::atoi(argv[idx ++]);
    int vxp_test_y = std::atoi(argv[idx ++]);
    int vxc_test_y = std::atoi(argv[idx ++]);
    int vyp_test_x = std::atoi(argv[idx ++]);
    int vyc_test_x = std::atoi(argv[idx ++]);
    int vyp_test_y = std::atoi(argv[idx ++]);
    int vyc_test_y = std::atoi(argv[idx ++]);
    int pp_test_x  = std::atoi(argv[idx ++]);
    int pc_test_x  = std::atoi(argv[idx ++]);
    int pp_test_y  = std::atoi(argv[idx ++]);
    int pc_test_y  = std::atoi(argv[idx ++]);

    using std::max;
    int pmax_trial = max({vxp_trial_x, vxp_trial_y, vyp_trial_x, vyp_trial_y, pp_trial_x, pp_trial_y});
    int pmax_test = max({vxp_test_x, vxp_test_y, vyp_test_x, vyp_test_y, pp_test_x, pp_test_y});
    int p_max = max(pmax_trial, pmax_test);

    int quad = p_max + 1;

    double dt = 2.0 / nsteps;
    timesteps_config steps{ nsteps, dt };
    int ders = 2;

    // Trial
    auto U1_trial_x = make_dimension(vxp_trial_x, vxc_trial_x, n, quad, ders);
    auto U1_trial_y = make_dimension(vxp_trial_y, vxc_trial_y, n, quad, ders);

    auto U2_trial_x = make_dimension(vyp_trial_x, vyc_trial_x, n, quad, ders);
    auto U2_trial_y = make_dimension(vyp_trial_y, vyc_trial_y, n, quad, ders);

    auto trial_x = make_dimension(pp_trial_x, pc_trial_x, n, quad, ders);
    auto trial_y = make_dimension(pp_trial_y, pc_trial_y, n, quad, ders);

    // Test
    auto U1_test_x = make_dimension(vxp_test_x, vxc_test_x, n, quad, ders);
    auto U1_test_y = make_dimension(vxp_test_y, vxc_test_y, n, quad, ders);

    auto U2_test_x = make_dimension(vyp_test_x, vyc_test_x, n, quad, ders);
    auto U2_test_y = make_dimension(vyp_test_y, vyc_test_y, n, quad, ders);

    auto test_x = make_dimension(pp_test_x, pc_test_x, n, quad, ders);
    auto test_y = make_dimension(pp_test_y, pc_test_y, n, quad, ders);

    auto trial = space_set{
        U1_trial_x, U1_trial_y,
        U2_trial_x, U2_trial_y,
        trial_x, trial_y
    };

    auto test = space_set{
        U1_test_x, U1_test_y,
        U2_test_x, U2_test_y,
        test_x, test_y
    };

    // Sanity check
    auto trial_dim = total_dimension(trial);
    auto test_dim = total_dimension(test);

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
    }

    using Prob = prob_manufactured_NS_nonpoly;
    // using Prob = prob_cavity_flow;
    Prob problem{Re};
    auto sim = stokes_projection<Prob>{trial, test, steps, problem};
    sim.run();
}
