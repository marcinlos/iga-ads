#include "shishkin.hpp"
#include "advection.hpp"

using namespace ads;

void validate_dimensions(const ads::dimension& trial, const ads::dimension& test) {
    auto trial_dim = trial.B.dofs();
    auto test_dim = test.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
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


int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: cg <Nx> <Ny> <adapt> <p_trial> <C_trial> <p_test> <C_test>" << std::endl;
        std::exit(1);
    }
    int nx = std::atoi(argv[1]);
    int ny = std::atoi(argv[2]);

    bool adapt = std::atoi(argv[3]);

    int p_trial = std::atoi(argv[4]);
    int C_trial = std::atoi(argv[5]);
    int p_test = std::atoi(argv[6]);
    int C_test = std::atoi(argv[7]);

    double peclet = 1e2;

    std::cout << "trial (" << p_trial << ", " << C_trial << "), "
              << "test (" << p_test << ", " << C_test << ")" << std::endl;
    std::cout << "adaptations: " << std::boolalpha << adapt << std::endl;

    int quad = std::max(p_trial, p_test) + 1;

    bool adapt_x = true  && adapt;
    bool adapt_y = false && adapt;

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

    validate_dimensions(dtrial_x, dtest_x);
    print_dofs(dtrial_x, dtest_x);

    auto sim = advection{dtrial_x, dtrial_y, dtest_x, dtest_y, peclet};
    sim.run();
}
