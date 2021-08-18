// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/solver.hpp"

#include <cstddef>
#include <vector>

#include <catch2/catch.hpp>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"
#include "ads/lin/solver_ctx.hpp"

using matrix = ads::lin::band_matrix;

template <std::size_t N>
using tensor = ads::lin::tensor<double, N>;

template <std::size_t N>
auto to_vector(const tensor<N>& a) -> std::vector<double> {
    auto const* data = a.data();
    auto const size = a.size();
    return {data, data + size};
}

using Catch::Matchers::Approx;

TEST_CASE("ADS one dimension", "[ads]") {
    auto const n = 4;

    auto M = matrix{1, 1, n};
    M(0, 0) = 1;
    M(0, 1) = 2;
    M(1, 0) = 2;
    M(1, 1) = 3;
    M(1, 2) = 1;
    M(2, 1) = -1;
    M(2, 2) = 4;
    M(3, 2) = 1;
    M(3, 3) = 3;

    auto ctx = ads::lin::solver_ctx{M};
    ads::lin::factorize(M, ctx);

    auto rhs = tensor<1>{{n}};
    rhs(0) = 5;
    rhs(1) = 11;
    rhs(2) = 10;
    rhs(3) = 15;

    auto dim = ads::dim_data{M, ctx};
    auto const expected = std::vector<double>{1, 2, 3, 4};

    SECTION("without buffer") {
        ads_solve(rhs, dim);
        CHECK_THAT(to_vector(rhs), Approx(expected));
    }

    SECTION("with buffer") {
        auto buf = tensor<1>{{n}};
        ads_solve(rhs, buf, dim);
        CHECK_THAT(to_vector(rhs), Approx(expected));
    }
}

TEST_CASE("ADS standard 2D", "[ads]") {
    auto const nx = 4;
    auto const ny = 3;

    auto Mx = matrix{1, 1, nx};
    auto My = matrix{1, 1, ny};

    Mx(0, 0) = 1;
    Mx(0, 1) = 2;
    Mx(1, 0) = 2;
    Mx(1, 1) = 3;
    Mx(1, 2) = 1;
    Mx(2, 1) = -1;
    Mx(2, 2) = 4;
    Mx(3, 2) = 1;
    Mx(3, 3) = 3;

    My(0, 0) = 2;
    My(0, 1) = 3;
    My(1, 0) = 1;
    My(1, 1) = 2;
    My(1, 2) = -3;
    My(2, 1) = 2;
    My(2, 2) = 2;

    auto ctx_x = ads::lin::solver_ctx{Mx};
    ads::lin::factorize(Mx, ctx_x);

    auto ctx_y = ads::lin::solver_ctx{My};
    ads::lin::factorize(My, ctx_y);

    auto rhs = tensor<2>{{nx, ny}};

    rhs(0, 0) = 61;
    rhs(1, 0) = 127;
    rhs(2, 0) = 86;
    rhs(3, 0) = 123;

    rhs(0, 1) = -48;
    rhs(1, 1) = -96;
    rhs(2, 1) = -48;
    rhs(3, 1) = -64;

    rhs(0, 2) = 92;
    rhs(1, 2) = 188;
    rhs(2, 2) = 112;
    rhs(3, 2) = 156;

    auto dim_x = ads::dim_data{Mx, ctx_x};
    auto dim_y = ads::dim_data{My, ctx_y};

    auto const expected = std::vector<double>{
        1, 2,  3,  4,   //
        5, 6,  7,  8,   //
        9, 10, 11, 12,  //
    };

    auto buf = tensor<2>{{nx, ny}};
    ads_solve(rhs, buf, dim_x, dim_y);

    CHECK_THAT(to_vector(rhs), Approx(expected));
}

TEST_CASE("ADS standard 3D", "[ads]") {
    auto const nx = 4;
    auto const ny = 3;
    auto const nz = 2;

    auto Mx = matrix{1, 1, nx};
    auto My = matrix{1, 1, ny};
    auto Mz = matrix{1, 1, nz};

    Mx(0, 0) = 1;
    Mx(0, 1) = 2;
    Mx(1, 0) = 2;
    Mx(1, 1) = 3;
    Mx(1, 2) = 1;
    Mx(2, 1) = -1;
    Mx(2, 2) = 4;
    Mx(3, 2) = 1;
    Mx(3, 3) = 3;

    My(0, 0) = 2;
    My(0, 1) = 3;
    My(1, 0) = 1;
    My(1, 1) = 2;
    My(1, 2) = -3;
    My(2, 1) = 2;
    My(2, 2) = 2;

    Mz(0, 0) = 1;
    Mz(0, 1) = 2;
    Mz(1, 0) = 3;
    Mz(1, 1) = 4;

    auto ctx_x = ads::lin::solver_ctx{Mx};
    ads::lin::factorize(Mx, ctx_x);

    auto ctx_y = ads::lin::solver_ctx{My};
    ads::lin::factorize(My, ctx_y);

    auto ctx_z = ads::lin::solver_ctx{Mz};
    ads::lin::factorize(Mz, ctx_z);

    auto rhs = tensor<3>{{nx, ny, nz}};

    auto dim_x = ads::dim_data{Mx, ctx_x};
    auto dim_y = ads::dim_data{My, ctx_y};
    auto dim_z = ads::dim_data{Mz, ctx_z};

    rhs(0, 0, 0) = 543;
    rhs(1, 0, 0) = 1101;
    rhs(2, 0, 0) = 618;
    rhs(3, 0, 0) = 849;

    rhs(0, 1, 0) = -144;
    rhs(1, 1, 0) = -288;
    rhs(2, 1, 0) = -144;
    rhs(3, 1, 0) = -192;

    rhs(0, 2, 0) = 564;
    rhs(1, 2, 0) = 1140;
    rhs(2, 2, 0) = 624;
    rhs(3, 2, 0) = 852;

    rhs(0, 0, 1) = 1147;
    rhs(1, 0, 1) = 2329;
    rhs(2, 0, 1) = 1322;
    rhs(3, 0, 1) = 1821;

    rhs(0, 1, 1) = -336;
    rhs(1, 1, 1) = -672;
    rhs(2, 1, 1) = -336;
    rhs(3, 1, 1) = -448;

    rhs(0, 2, 1) = 1220;
    rhs(1, 2, 1) = 2468;
    rhs(2, 2, 1) = 1360;
    rhs(3, 2, 1) = 1860;

    auto const expected = std::vector<double>{
        1,  2,  3,  4,   //
        5,  6,  7,  8,   //
        9,  10, 11, 12,  //

        13, 14, 15, 16,  //
        17, 18, 19, 20,  //
        21, 22, 23, 24,  //
    };

    auto buf = tensor<3>{{nx, ny, nz}};
    ads_solve(rhs, buf, dim_x, dim_y, dim_z);

    CHECK_THAT(to_vector(rhs), Approx(expected));
}
