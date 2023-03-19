// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/lin/dense_solve.hpp"

#include <algorithm>

#include <catch2/catch.hpp>

#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/tensor.hpp"

namespace lin = ads::lin;

TEST_CASE("Dense matrix") {
    SECTION("Solver") {
        int kl = 1;
        int ku = 2;
        int n = 6;
        int d = 4;

        lin::dense_matrix m(n, n);
        lin::matrix b({n, d});

        for (int i = 0; i < n; ++i) {
            for (int j = std::max(0, i - kl); j < std::min(n, i + ku + 1); ++j) {
                m(i, j) = (i + 1) * 10 + j + 1;
            }
            for (int j = 0; j < d; ++j) {
                b(i, j) = (j + 1) * (i + 1);
            }
        }

        double solution[] = {0.230377, -0.126052, -0.0016554, -0.00111222, 0.203603, -0.109609};
        lin::matrix x({n, d});
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < d; ++j) {
                x(i, j) = (j + 1) * solution[i];
            }
        }

        lin::solver_ctx ctx(m);
        lin::factorize(m, ctx);
        lin::solve_with_factorized(m, b, ctx);

        CHECK(approx_equal(x, b, 1e-5));
    }
}
