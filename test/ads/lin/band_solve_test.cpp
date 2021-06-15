#include "ads/lin/band_solve.hpp"

#include <algorithm>

#include <catch2/catch.hpp>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/tensor.hpp"


using namespace ads::lin;

TEST_CASE("Banded matrix") {

    SECTION("Solver") {
        int kl = 1;
        int ku = 2;
        int n = 6;
        int d = 4;

        band_matrix m(kl, ku, n);
        matrix b({n, d});

        for (int i = 0; i < n; ++ i) {
            for (int j = std::max(0, i - kl); j < std::min(n, i + ku + 1); ++ j) {
                m(i, j) = (i + 1) * 10 + j + 1;
            }
            for (int j = 0; j < d; ++ j) {
                b(i, j) = (j + 1) * (i + 1);
            }
        }

        double solution[] = { 0.230377, -0.126052, -0.0016554, -0.00111222, 0.203603, -0.109609 };
        matrix x({n, d});
        for (int i = 0; i < n; ++ i) {
            for (int j = 0; j < d; ++ j) {
                x(i, j) = (j + 1) * solution[i];
            }
        }

        solver_ctx ctx(m);
        solve(m, b, ctx);

        CHECK(approx_equal(x, b, 1e-5));
    }

    SECTION("Vector multiplication") {
        band_matrix A{ 1, 1, 4, 3 };
        A(0, 0) = A(1, 1) = A(2, 2) = A(3, 2) = 1;
        std::vector<double> x{1, 2, 3, 4, 5, 6}, y(4*2);

        multiply(A, x, y, 2);

        auto expected = std::vector<double>{1, 2, 3, 3, 4, 5, 6, 6};
        CHECK_THAT(y, Catch::Matchers::Approx(expected));
    }
}

