#include <algorithm>
#include <libunittest/all.hpp>
#include "ads/lin/tensor.hpp"
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"


using namespace unittest::assertions;

namespace ads {
namespace lin {


struct banded_solver_test : unittest::testcase<> {

    static void run() {
        UNITTEST_CLASS(banded_solver_test)
        UNITTEST_RUN(test_example_matrix)
    }

    void test_example_matrix() {
        int kl = 1;
        int ku = 2;
        int n = 6;
        int d = 4;

        band_matrix m(kl, ku, n);
        matrix b(n, d);

        for (int i = 0; i < n; ++ i) {
            for (int j = std::max(0, i - kl); j < std::min(n, i + ku + 1); ++ j) {
                m(i, j) = (i + 1) * 10 + j + 1;
            }
            for (int j = 0; j < d; ++ j) {
                b(i, j) = (j + 1) * (i + 1);
            }
        }

        double solution[] = { 0.230377, -0.126052, -0.0016554, -0.00111222, 0.203603, -0.109609 };
        matrix x(n, d);
        for (int i = 0; i < n; ++ i) {
            for (int j = 0; j < d; ++ j) {
                x(i, j) = (j + 1) * solution[i];
            }
        }

        solver_ctx ctx(m);
        solve(m, b, ctx);

        assert_true(approx_equal(x, b, 1e-5));
    }
};

REGISTER(banded_solver_test)

}
}
