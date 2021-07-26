// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/basis_data.hpp"

#include <catch2/catch.hpp>

#include "ads/bspline/bspline.hpp"

TEST_CASE("basis_data memory management") {
    auto basis = ads::bspline::create_basis(0.0, 1.0, 2, 5);
    auto data = ads::basis_data{basis, 1};

    SECTION("can be copied") {
        auto copy = data;

        REQUIRE(copy.degree == data.degree);
        REQUIRE(copy.elements == data.elements);
        REQUIRE(copy.dofs == data.dofs);
        REQUIRE(copy.derivatives == data.derivatives);
        REQUIRE(copy.quad_order == data.quad_order);
        REQUIRE(copy.elem_division == data.elem_division);
    }
}
