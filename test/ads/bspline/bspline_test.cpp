#include <catch2/catch.hpp>
#include "ads/bspline/bspline.hpp"


using namespace ads::bspline;
using namespace Catch::Matchers;

TEST_CASE("B-spline basis", "[splines]") {

    basis b = create_basis(0.0, 1.0, 2, 4);

    SECTION("Knot vector is correct after creation") {
        REQUIRE_THAT(b.knot, Equals<double>({0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1}));
    }

    SECTION("Finding span") {
        SECTION("point outside domain") {
            // degree = 2 is the lowest possible result
            CHECK(find_span(-1, b) == 2);
            CHECK(find_span(2.0, b) == 5);
        }
        SECTION("point on domain boundary") {
            CHECK(find_span(0.0, b) == 2);
            CHECK(find_span(1.0, b) == 5);
        }
        SECTION("point on element boundary") {
            CHECK(find_span(0.25, b) == 3);
            CHECK(find_span(0.75, b) == 5);
        }
        SECTION("point in element interior") {
            CHECK(find_span(0.1, b) == 2);
            CHECK(find_span(0.3, b) == 3);
            CHECK(find_span(0.7, b) == 4);
            CHECK(find_span(0.9, b) == 5);
        }
    }

    SECTION("First non-zero dofs") {
        REQUIRE_THAT(first_nonzero_dofs(b), Equals<int>({0, 1, 2, 3}));
    }
}

