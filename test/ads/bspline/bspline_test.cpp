#include <libunittest/all.hpp>
#include "ads/bspline/bspline.hpp"

using namespace unittest::assertions;
using namespace ads::bspline;

struct bspline_test: unittest::testcase<> {

    static void run() {
        UNITTEST_CLASS(bspline_test)
        UNITTEST_RUN(create_basis_test)
        UNITTEST_RUN(find_span_test)
    }

    void create_basis_test() {
        basis b = create_basis(0, 1, 2, 4);
        double expected[] = { 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1 };
        assert_equal_containers(b.knot, expected, "quadratic basis");
    }

    void find_span_test() {
        basis b = create_basis(0, 1, 2, 4);

        assert_equal(2, find_span(-1, b), "x = -1"); // degree = 2 is the lowest possible result
        assert_equal(2, find_span(0, b), "x = 0");
        assert_equal(2, find_span(0.1, b), "x = 0.1");
        assert_equal(3, find_span(0.25, b), "x = 0.25");
        assert_equal(3, find_span(0.3, b), "x = 0.3");
        assert_equal(4, find_span(0.7, b), "x = 0.7");
        assert_equal(5, find_span(0.75, b), "x = 0.75");
        assert_equal(5, find_span(0.9, b), "x = 0.9");
        assert_equal(5, find_span(1, b), "x = 1");
        assert_equal(5, find_span(2, b), "x = 2");
    }
};

REGISTER(bspline_test)
