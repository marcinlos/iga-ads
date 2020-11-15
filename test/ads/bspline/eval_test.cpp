#include <catch2/catch.hpp>
#include "ads/bspline/bspline.hpp"
#include "ads/util.hpp"
#include <iomanip>
#include <numeric>


using namespace ads::bspline;

TEST_CASE("B-spline evaluation", "[splines]") {
    const int p = 2;
    const int elements = 5;
    const int dofs = elements + p;

    double a = 0, b = 1;
    basis basis = create_basis(a, b, p, elements);
    eval_ctx ctx(p);

    SECTION("Partition of unity property holds") {
        const int N = 100;
        double y[N + 1][dofs] = {0};

        for (int i = 0; i <= N; ++ i) {
            double x = ads::lerp(i, N, a, b);
            int span = find_span(x, basis);

            int offset = span - p;
            eval_basis(span, x, basis, y[i] + offset, ctx);
        }

        for (int i = 0; i <= N; ++ i) {
            double sum = std::accumulate(std::begin(y[i]), std::end(y[i]), 0.0);
            INFO("i = " << i << ", x = " << ads::lerp(i, N, a, b));
            REQUIRE(sum == Approx(1.0));
        }
    }

    SECTION("Computing derivatives") {
        const int d = 2;
        const int N = 100;
        double y[N + 1][d + 1][dofs] = {0};

        for (int i = 0; i <= N; ++ i) {
            double x = ads::lerp(i, N, a, b);
            int span = find_span(x, basis);

            int offset = span - p;

            std::vector<double*> ys(d + 1);
            for (int j = 0; j <= d; ++ j) {
                ys[j] = y[i][j] + offset;
            }
            eval_basis_with_derivatives(span, x, basis, ys.data(), d, ctx);
        }

        SECTION("Partition of unity property holds") {
            for (int i = 0; i <= N; ++ i) {
                double sum = std::accumulate(std::begin(y[i][0]), std::end(y[i][0]), 0.0);
                INFO("i = " << i << ", x = " << ads::lerp(i, N, a, b));
                REQUIRE(sum == Approx(1.0));
            }
        }

        SECTION("Derivatives of sum of basis functions are zero") {
            for (int i = 0; i <= N; ++ i) {
                for (int j = 1; j <= d; ++ j) {
                    double sum = std::accumulate(std::begin(y[i][j]), std::end(y[i][j]), 0.0);
                    INFO("order: " << j << ", i = " << i << ", x = " << ads::lerp(i, N, a, b));
                    REQUIRE(sum == Approx(0.0).margin(1e-7));
                }
            }
        }
    }
}

