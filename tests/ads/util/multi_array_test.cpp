// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/util/multi_array.hpp"

#include <catch2/catch_all.hpp>

using Catch::Matchers::Equals;

TEST_CASE("Multiarray") {
    int n = 2;
    int m = 5;
    auto buffer = std::vector<int>(n * m);
    ads::multi_array_wrapper<int, 2, int*> a(buffer.data(), {n, m});

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a(i, j) = 10 * (i + 1) + j + 1;
        }
    }

    SECTION("Buffer is filled correctly") {
        auto expected = std::vector<int>({11, 12, 13, 14, 15, 21, 22, 23, 24, 25});
        CHECK_THAT(buffer, Equals(expected));
    }

    SECTION("Reshaping") {
        ads::multi_array_wrapper<int, 1, int*> a1(buffer.data(), {n * m});
        a1(7) = 7;
        ads::multi_array_wrapper<int, 2, int*> a2 = ads::reshape<2>(a, n, m);
        CHECK(a2(1, 2) == 7);
    }

    SECTION("Access through const ref") {
        const auto& ca = a;
        CHECK(ca(0, 1) == 12);
        // ca(0, 1) = 3; // doesn't compile
    }

    SECTION("Size") {
        const int p = 2;
        const int q = 5;
        const int r = 3;
        const int s = 2;
        auto large_buffer = std::vector<int>(p * q * r * s);
        ads::multi_array_wrapper<int, 4, int*> aa(large_buffer.data(), {p, q, r, s});

        CHECK(aa.size(0) == p);
        CHECK(aa.size(1) == q);
        CHECK(aa.size(2) == r);
        CHECK(aa.size(3) == s);
    }
}
