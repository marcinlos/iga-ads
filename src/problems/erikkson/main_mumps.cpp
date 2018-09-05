#include "problems/erikkson/erikkson_mumps.hpp"
#include "problems/erikkson/erikkson_supg.hpp"
#include "problems/erikkson/erikkson_cg.hpp"
#include "problems/erikkson/pollution_rotation.hpp"

// #include "problems/erikkson/erikkson_quanling.hpp"
// #include "problems/erikkson/pollution_cg.hpp"
// #include "problems/erikkson/erikkson_mumps_split.hpp"


using namespace ads;

double shishkin_const(int n, double eps) {
    return std::log(n) * eps;
}

bspline::basis create_basis(double a, double b, int p, int elements, int repeated_nodes, bool adapt, double d) {
    int points = elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    bspline::knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }

    auto x0 = 0.5;
    auto y0 = 1 - d;

    for (int i = 1; i < points - 1; ++i) {
        auto t = lerp(i, elements, 0.0, 1.0);

        auto s = adapt ? (t < x0 ? t / x0 * y0 : (t - x0) / (1 - x0) * (1 - y0) + y0) : t;
        for (int j = 0; j < r; ++ j) {
            knot[p + 1 + (i - 1) * r + j] = lerp(s, a, b);
        }
    }

    return {std::move(knot), p};
}

bspline::basis create_checkboard_basis(double a, double b, int p, int elements, int repeated_nodes, bool adapt) {
    int total_elements = 2 * elements + 2;
    int points = total_elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    bspline::knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }

    double eps = 1e-2;
    double x0 = 0;
    double x1 = 0.5 - eps;
    double x2 = 0.5;
    double x3 = 1 - eps;

    for (int i = 1; i < points - 1; ++i) {

        double s;
        if (i <= elements) {
            auto t = lerp(i, elements, 0.0, 1.0);
            s = lerp(t, x0, x1);
        } else {
            auto t = lerp(i - elements - 1, elements, 0.0, 1.0);
            s = lerp(t, x2, x3);
        }
        std::cout << s << std::endl;

        for (int j = 0; j < r; ++ j) {
            knot[p + 1 + (i - 1) * r + j] = lerp(s, a, b);
        }
    }

    return {std::move(knot), p};
}


int main(int argc, char* argv[]) {
    if (argc != 10) {
        std::cerr << "Usage: erikkson_mumps <type> <N> <subdivision> <adapt> <p_trial> <C_trial> <p_test> <C_test> <steps>" << std::endl;
        std::exit(1);
    }
    const std::string type{ argv[1] };
    int n = std::atoi(argv[2]);
    int subdivision = std::atoi(argv[3]);
    bool adapt = std::atoi(argv[4]);

    int p_trial = std::atoi(argv[5]);
    int C_trial = std::atoi(argv[6]);
    int p_test = std::atoi(argv[7]);
    int C_test = std::atoi(argv[8]);
    int nsteps = std::atoi(argv[9]);

    std::cout << "trial (" << p_trial << ", " << C_trial << "), "
              << "test (" << p_test << ", " << C_test << ")" << std::endl;

    // double S = 5000.0;
    double S = 1.0;

    int quad = std::max(p_trial, p_test) + 1;

    std::cout << "adaptations: " << std::boolalpha << adapt << std::endl;

    timesteps_config steps{ nsteps, 1e-1 };
    int ders = 2;

    bool adapt_x = true && adapt;
    bool adapt_y = true && adapt;

    // auto d = 1e-2;
    auto d = shishkin_const(n, 1e-2);

    auto trial_basis_x = create_basis(-S, S, p_trial, n, p_trial - 1 - C_trial, adapt_x, d);
    // auto trial_basis_x = create_checkboard_basis(0, S, p_trial, n, p_trial - 1 - C_trial, adapt_x);
    auto dtrial_x = dimension{ trial_basis_x, quad, ders, subdivision };

    auto trial_basis_y = create_basis(-S, S, p_trial, n, p_trial - 1 - C_trial, adapt_y, d);
    // auto trial_basis_y = create_checkboard_basis(0, S, p_trial, n, p_trial - 1 - C_trial, adapt_y);
    auto dtrial_y = dimension{ trial_basis_y, quad, ders, subdivision };

    auto test_basis_x = create_basis(-S, S, p_test, subdivision*n, p_test - 1 - C_test, adapt_x, d);
    // auto test_basis_x = create_checkboard_basis(0, S, p_test, subdivision*n, p_test - 1 - C_test, adapt_x);
    auto dtest_x = dimension{ test_basis_x, quad, ders, 1 };

    auto test_basis_y = create_basis(-S, S, p_test, subdivision*n, p_test - 1 - C_test, adapt_y, d);
    // auto test_basis_y = create_checkboard_basis(0, S, p_test, subdivision*n, p_test - 1 - C_test, adapt_y);
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

    if (type == "igrm-mumps") erikkson_mumps{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();
    if (type == "igrm-cg") erikkson_CG{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();
    if (type == "supg") erikkson_supg {dtrial_x, dtrial_y, steps}.run();
    if (type == "pollution") pollution_rotation{dtrial_x, dtrial_y, dtest_x, dtest_y, steps}.run();


    // erikkson_mumps_split sim{dtrial_x, dtrial_y, dtest_x, dtest_y, steps};
    // pollution_CG sim{dtrial_x, dtrial_y, dtest_x, dtest_y, steps};
    // erikkson_quanling sim{dtrial_x, dtrial_y, dtest_x, dtest_y, steps};

    // sim.run();
}
