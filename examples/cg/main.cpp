// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <lyra/lyra.hpp>

#include "ads/util.hpp"
#include "advection.hpp"
#include "problems.hpp"
#include "shishkin.hpp"

ads::bspline::basis basis_from_points(const std::vector<double>& points, int p,
                                      int repeated_nodes) {
    auto elems = ads::narrow_cast<int>(points.size()) - 1;
    int r = repeated_nodes + 1;
    int size = (elems - 1) * r + 2 * (p + 1);
    auto knot = ads::bspline::knot_vector(size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = points[0];
        knot[size - i - 1] = points[elems];
    }

    for (int i = 1; i < elems; ++i) {
        for (int j = 0; j < r; ++j) {
            knot[p + 1 + (i - 1) * r + j] = points[i];
        }
    }
    return {std::move(knot), p};
}

std::vector<double> make_points(int n) {
    std::vector<double> points{};

    points.push_back(0);

    double x = 0;
    double h = 1;

    for (int i = 0; i < n; ++i) {
        h /= 2;
        x += h;
        points.push_back(x);
    }

    points.push_back(1);
    return points;
}

void validate_dimensions(const ads::dimension& trial, const ads::dimension& test, bool print_dim) {
    auto trial_dim = trial.B.dofs();
    auto test_dim = test.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space (" << trial_dim
                  << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else if (print_dim) {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
    }
}

void print_dofs(const ads::dimension& trial, const ads::dimension& test) {
    auto trial_dim = trial.B.dofs();
    auto test_dim = test.B.dofs();

    int dofs_trial = trial_dim * trial_dim;
    int dofs_test = test_dim * test_dim;
    int dofs = dofs_trial + dofs_test;
    std::cout << "DOFs : " << dofs << std::endl;
}

auto make_config_parser(ads::advection_config& cfg) {
    return lyra::opt(cfg.tol_outer, "eps")["--tol-outer"]  //
           ("outer iterations tolerance")
         | lyra::opt(cfg.tol_inner, "eps")["--tol-inner"]  //
           ("inner iterations tolerance")
         | lyra::opt(cfg.max_outer_iters, "N")["--max-outer-iters"]  //
           ("maximum # of outer iterations")
         | lyra::opt(cfg.max_inner_iters, "N")["--max-inner-iters"]  //
           ("maximum # of inner iterations")
         | lyra::opt(cfg.use_cg)["--cg"]  //
           ("use CG solver?")
         | lyra::opt(cfg.weak_bc)["--weak-bc"]  //
           ("impose Dirichlet BC weakly?")
         | lyra::opt(cfg.print_inner)["--print-inner"]  //
           ("print errors for inner iterations?")
         | lyra::opt(cfg.print_inner_count)["--print-inner-count"]  //
           ("print number of inner iterations?")
         | lyra::opt(cfg.print_outer)["--print-outer"]  //
           ("print errors for outer iterations?")
         | lyra::opt(cfg.print_outer_count)["--print-outer-count"]  //
           ("print number of outer iterations?")
         | lyra::opt(cfg.print_inner_total)["--print-inner-total"]  //
           ("print total number of inner iterations?")
         | lyra::opt(cfg.print_times)["--print-times"]  //
           ("print integration and solver times?")
         | lyra::opt(cfg.print_errors)["--print-errors"]  //
           ("print solution errors?")
         | lyra::opt(cfg.plot)["--plots"]  //
           ("save data for plots?")
         | lyra::opt(cfg.threads, "N")["--threads"]  //
           ("number of threads");
}

template <typename Fun>
void with_problem(const std::string& name, double peclet, Fun&& fun) {
    auto epsilon = 1 / peclet;

    if (name == "erikkson") {
        fun(erikkson{epsilon});
    } else if (name == "problem1") {
        fun(problem1{epsilon});
    } else {
        std::cerr << "Unknown problem: " << name << std::endl;
        std::exit(-1);
    }
}

int main(int argc, char* argv[]) {
    std::string problem;
    int nx;
    int ny;
    int p_trial;
    int C_trial;
    int p_test;
    int C_test;

    double peclet = 1e2;
    double eta = 0.1;
    bool adapt_x = false;
    bool adapt_y = false;

    bool print_dof_count = false;
    bool print_dims = false;

    ads::advection_config cfg;

    bool show_help = false;
    auto const cli =                                                                 //
        lyra::help(show_help)                                                        //
        | lyra::arg(problem, "problem")("problem to solve").required()               //
        | lyra::arg(nx, "Nx")("mesh resolution in x direction").required()           //
        | lyra::arg(ny, "Ny")("mesh resolution in y direction").required()           //
        | lyra::arg(p_trial, "p trial")("B-spline order of trial space").required()  //
        | lyra::arg(C_trial, "C trial")("continuity of trial space").required()      //
        | lyra::arg(p_test, "p test")("B-spline order of test space").required()     //
        | lyra::arg(C_test, "C test")("continuity of test space").required()         //
        | lyra::opt(adapt_x)["--adaptx"]("adapt in x direction")                     //
        | lyra::opt(adapt_y)["--adapty"]("adapt in y direction")                     //
        | lyra::opt(peclet, "val")["--Pe"]("Peclet number")                          //
        | lyra::opt(eta, "val")["--eta"]("eta (solver parameter)")                   //
        | lyra::opt(print_dof_count)["--dofs"]("print # of DOFs")                    //
        | lyra::opt(print_dims)["--dims"]("print dimensions of spacs")               //
        | make_config_parser(cfg);

    auto const result = cli.parse({argc, argv});

    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }

    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }

    int quad = std::max(p_trial, p_test) + 1;

    int rep_trial = p_trial - 1 - C_trial;
    int rep_test = p_test - 1 - C_test;

    auto bd_layer = shishkin_const(nx, 1 / peclet);

    auto make_dim = [&](int p, int n, int rep, bool adapt) {
        int ders = 1;
        auto basis = create_basis(0.0, 1.0, p, n, rep, adapt, bd_layer);
        return ads::dimension{basis, quad, ders, 1};
    };

    auto points_x = make_points(nx);

    auto trial_basis_x = basis_from_points(points_x, p_trial, rep_trial);
    auto dtrial_x = ads::dimension{trial_basis_x, quad, 1, 1};
    // auto dtrial_x = make_dim(p_trial, nx, rep_trial, adapt_x);
    auto dtrial_y = make_dim(p_trial, ny, rep_trial, adapt_y);

    auto test_basis_x = basis_from_points(points_x, p_test, rep_test);
    auto dtest_x = ads::dimension{test_basis_x, quad, 1, 1};
    // auto dtest_x = make_dim(p_test, nx, rep_test, adapt_x);
    auto dtest_y = make_dim(p_test, ny, rep_test, adapt_y);

    validate_dimensions(dtrial_x, dtest_x, print_dims);

    if (print_dof_count)
        print_dofs(dtrial_x, dtest_x);

    with_problem(problem, peclet, [&](auto prob) {
        using Problem = decltype(prob);
        auto sim =
            ads::advection<Problem>{dtrial_x, dtrial_y, dtest_x, dtest_y, peclet, eta, cfg, prob};
        sim.run();
    });
}
