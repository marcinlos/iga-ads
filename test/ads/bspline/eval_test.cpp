#include <libunittest/all.hpp>
#include "ads/bspline/bspline.hpp"
#include "ads/util.hpp"
#include <iomanip>
#include <numeric>

using namespace unittest::assertions;
using namespace ads::bspline;

struct bspline_eval_test: unittest::testcase<> {

    static void run() {
        UNITTEST_CLASS(bspline_eval_test)
        UNITTEST_RUN(eval_test)
        UNITTEST_RUN(eval_derivatives_test)
    }

    void eval_test() {
        const int p = 2;
        const int elements = 5;
        const int dofs = elements + p;

        double a = 0, b = 1;
        basis basis = create_basis(a, b, p, elements);
        eval_ctx ctx(p);

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
            assert_approx_equal(1, sum, 1e-5, "partition of unity property violated");
        }
    }

    void eval_derivatives_test() {
        const int p = 4;
        const int elements = 5;
        const int dofs = elements + p;
        const int d = 2;

        double a = 0, b = 1;
        basis basis = create_basis(a, b, p, elements);
        eval_ctx ctx(p);

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

        for (int i = 0; i <= N; ++ i) {
            double sum = std::accumulate(std::begin(y[i][0]), std::end(y[i][0]), 0.0);
            assert_approx_equal(1, sum, 1e-5, "partition of unity property violated");

            for (int j = 1; j <= d; ++ j) {
                double sum = std::accumulate(std::begin(y[i][j]), std::end(y[i][j]), 0.0);
                assert_approx_equal(0, sum, 1e-5, "derivative of sum should be zero");
            }
        }
    }

};

REGISTER(bspline_eval_test)
