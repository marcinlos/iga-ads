#include "shishkin.hpp"
#include "advection.hpp"
#include "problems.hpp"
#include <clara.hpp>

using namespace ads;
using namespace clara;

void validate_dimensions(const ads::dimension& trial, const ads::dimension& test, bool print_dim) {
    auto trial_dim = trial.B.dofs();
    auto test_dim = test.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
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

auto make_config_parser(advection_config& cfg) {
    return
        Opt(cfg.tol_outer, "outer iterations tolerance")["--tol-outer"] |
        Opt(cfg.tol_inner, "inner iterations tolerance")["--tol-inner"] |
        Opt(cfg.max_outer_iters, "maximum # of outer iterations")["--max-outer-iters"] |
        Opt(cfg.max_inner_iters, "maximum # of inner iterations")["--max-inner-iters"] |
        Opt(cfg.use_cg, "use CG solver?")["--cg"] |
        Opt(cfg.weak_bc, "impose Dirichlet BC weakly?")["--weak-bc"] |
        Opt(cfg.print_inner, "print errors for inner iterations?")["--print-inner"] |
        Opt(cfg.print_inner_count, "print number of inner iterations?")["--print-inner-count"] |
        Opt(cfg.print_outer, "print errors for outer iterations?")["--print-outer"] |
        Opt(cfg.print_outer_count, "print number of outer iterations?")["--print-outer-count"] |
        Opt(cfg.print_inner_total, "print total number of inner iterations?")["--print-inner-total"] |
        Opt(cfg.print_times, "print integration and solver times?")["--print-times"] |
        Opt(cfg.print_errors, "print solution errors?")["--print-errors"] |
        Opt(cfg.plot, "save data for plots?")["--plots"] |
        Opt(cfg.threads, "number of threads")["--threads"]
        ;
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

    advection_config cfg;

    bool help = false;
    auto cli = Help(help)
        | Arg(problem, "problem").required()
        | Arg(nx, "Nx").required()
        | Arg(ny, "Ny").required()
        | Arg(p_trial, "p trial").required()
        | Arg(C_trial, "C trial").required()
        | Arg(p_test, "p test").required()
        | Arg(C_test, "C test").required()
        | Opt(adapt_x, "adapt in x direction")["--adaptx"]
        | Opt(adapt_y, "adapt in y direction")["--adapty"]
        | Opt(peclet, "Peclet number")["--Pe"]
        | Opt(eta, "eta (solver parameter)")["--eta"]
        | Opt(print_dof_count, "print # of DOFs")["--dofs"]
        | Opt(print_dims, "print dimensions of spacs")["--dims"]
        | make_config_parser(cfg);

    auto result = cli.parse(Args(argc, argv));

    if (! result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::exit(1);
    }
    if (help) {
        cli.writeToStream(std::cout);
        return 0;
    }
    if (argc < 8) {
        std::cerr << "Usage: cg <problem> <Nx> <Ny> <p trial> <C trial> <p test> <C test>" << std::endl;
        std::exit(1);
    }

    int quad = std::max(p_trial, p_test) + 1;

    int rep_trial = p_trial - 1 - C_trial;
    int rep_test  = p_test - 1 - C_test;

    auto bd_layer = shishkin_const(nx, 1 / peclet);

    auto make_dim = [&](int p, int n, int rep, bool adapt) {
        int ders = 1;
        auto basis = create_basis(0.0, 1.0, p, n, rep, adapt, bd_layer);
        return dimension{basis, quad, ders, 1};
    };

    auto dtrial_x = make_dim(p_trial, nx, rep_trial, adapt_x);
    auto dtrial_y = make_dim(p_trial, ny, rep_trial, adapt_y);

    auto dtest_x = make_dim(p_test, nx, rep_test, adapt_x);
    auto dtest_y = make_dim(p_test, ny, rep_test, adapt_y);

    validate_dimensions(dtrial_x, dtest_x, print_dims);

    if (print_dof_count) print_dofs(dtrial_x, dtest_x);

    with_problem(problem, peclet, [&](auto prob) {
        using Problem = decltype(prob);
        auto sim = advection<Problem>{dtrial_x, dtrial_y, dtest_x, dtest_y, peclet, eta, cfg, prob};
        sim.run();
    });
}
