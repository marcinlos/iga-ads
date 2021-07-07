// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>

#include "erikkson_mumps_split.hpp"

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
    for (int i = 1; i < points - 1; ++i) {
        auto t = lerp(i, elements, 0.0, 1.0);

        for (int j = 0; j < r; ++j) {
            knot[p + 1 + (i - 1) * r + j] = lerp(t, a, b);
        }
    }

    return {std::move(knot), p};
}

int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: erikkson_nonstationary <type> <Nx> <Ny> <p_trial> <C_trial> <p_test> "
                     "<C_test> <steps>"
                  << std::endl;
        std::exit(1);
    }
    const std::string type{argv[1]};
    int nx = std::atoi(argv[2]);
    int ny = std::atoi(argv[3]);

    int p_trial = std::atoi(argv[4]);
    int C_trial = std::atoi(argv[5]);
    int p_test = std::atoi(argv[6]);
    int C_test = std::atoi(argv[7]);
    // auto dt = std::atof(argv[8]);
    int nsteps = std::atoi(argv[8]);

    double dt = 1.0 / nsteps;
    nsteps /= 2;
    nsteps += 1;

    // std::cout << "trial (" << p_trial << ", " << C_trial << "), "
    //           << "test (" << p_test << ", " << C_test << ")" << std::endl;

    double S = 1.0;
    int quad = std::max(p_trial, p_test) + 1;

    timesteps_config steps{nsteps, dt};
    int ders = 2;

    auto trial_basis_x = create_basis(0, S, p_trial, nx, p_trial - 1 - C_trial);
    auto dtrial_x = dimension{trial_basis_x, quad, ders, 1};

    auto trial_basis_y = create_basis(0, S, p_trial, ny, p_trial - 1 - C_trial);
    auto dtrial_y = dimension{trial_basis_y, quad, ders, 1};

    auto test_basis_x = create_basis(0, S, p_test, nx, p_test - 1 - C_test);
    auto dtest_x = dimension{test_basis_x, quad, ders, 1};

    auto test_basis_y = create_basis(0, S, p_test, ny, p_test - 1 - C_test);
    auto dtest_y = dimension{test_basis_y, quad, ders, 1};

    auto trial_dim = dtrial_x.B.dofs();
    auto test_dim = dtest_x.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space (" << trial_dim
                  << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
        // std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
    }

    // int dofsU = trial_dim * trial_dim;
    // int dofsV = test_dim * test_dim;
    // int dofs = dofsU + dofsV;
    // std::cout << "DOFs : " << dofs << std::endl;

    // std::cout << "/\\t = " << dt << std::endl;

    // if (type == "igrm-mumps") erikkson_mumps{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();
    // if (type == "igrm-cg") erikkson_CG{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();
    // if (type == "igrm-cg-weak") erikkson_CG_weak{dtrial_x, dtrial_y, dtest_x, dtest_y,
    // steps}.run(); if (type == "supg") erikkson_supg {dtrial_x, dtrial_y, steps}.run(); if (type
    // == "supg-weak") erikkson_supg_weak{dtrial_x, dtrial_y, steps}.run(); if (type == "pollution")
    // pollution_rotation{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run(); if (type == "split")
    // erikkson_mumps_split{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();
    scheme method;
    if (type == "BE")
        method = scheme::BE;
    else if (type == "CN")
        method = scheme::CN;
    else if (type == "PR")
        method = scheme::peaceman_rachford;
    else if (type == "strang-BE")
        method = scheme::strang_BE;
    else if (type == "strang-CN")
        method = scheme::strang_CN;

    else {
        std::cerr << "Unknown scheme: " << type << std::endl;
        std::exit(1);
    }

    erikkson_mumps_split sim{dtrial_x, dtrial_y, dtest_x, dtest_y, method, steps};
    sim.run();
}
