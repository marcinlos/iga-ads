// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/solver.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include <catch2/catch_all.hpp>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"
#include "ads/lin/solver_ctx.hpp"

using matrix = ads::lin::band_matrix;

template <std::size_t N>
using tensor = ads::lin::tensor<double, N>;

using Catch::Matchers::Approx;

template <std::size_t N>
auto to_vector(const tensor<N>& a) -> std::vector<double> {
    auto const* data = a.data();
    auto const size = a.size();
    return {data, data + size};
}

auto make_matrix(int kl, int ku, int n, std::vector<std::vector<double>> const& rows) -> matrix {
    auto M = matrix{kl, ku, n};
    for (int i = 0; i < n; ++i) {
        auto j = std::max(i - kl, 0);
        for (auto val : rows[i]) {
            M(i, j) = val;
            ++j;
        }
    }
    return M;
}

struct factored_matrix {
    using matrix = ads::lin::band_matrix;
    using context = ads::lin::solver_ctx;

    matrix mat;
    context ctx;

    explicit factored_matrix(matrix m)
    : mat(std::move(m))
    , ctx{mat} {
        ads::lin::factorize(mat, ctx);
    }
};

auto make_factored_matrix(int kl, int ku, int n, std::vector<std::vector<double>> const& rows)
    -> factored_matrix {
    auto M = make_matrix(kl, ku, n, rows);
    return factored_matrix{std::move(M)};
}

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

TEST_CASE("ADS generalized 2D", "[ads]") {
    auto const na = 2;
    auto const nb = 3;

    auto A = matrix{1, 1, na};
    auto B1 = matrix{1, 1, nb};
    auto B2 = matrix{1, 1, nb};

    A(0, 0) = 1;
    A(0, 1) = 2;
    A(1, 0) = 3;
    A(1, 1) = 4;

    B1(0, 0) = 2;
    B1(0, 1) = 3;
    B1(1, 0) = 1;
    B1(1, 1) = 2;
    B1(1, 2) = -3;
    B1(2, 1) = 2;
    B1(2, 2) = 2;

    B2(0, 0) = 1;
    B2(0, 1) = 2;
    B2(1, 0) = 3;
    B2(1, 1) = 1;
    B2(1, 2) = 2;
    B2(2, 1) = 2;
    B2(2, 2) = -1;

    auto ctx_A = ads::lin::solver_ctx{A};
    ads::lin::factorize(A, ctx_A);

    auto ctx_B1 = ads::lin::solver_ctx{B1};
    ads::lin::factorize(B1, ctx_B1);

    auto ctx_B2 = ads::lin::solver_ctx{B2};
    ads::lin::factorize(B2, ctx_B2);

    auto dim_A = ads::dim_data{A, ctx_A};

    auto const step_B = [&](auto& rhs) {
        auto* data = rhs.data();
        auto const rhs_size = rhs.size(0);

        solve_with_factorized(B1, data, ctx_B1, 1);
        solve_with_factorized(B2, data + rhs_size, ctx_B2, 1);
    };

    auto const expected = std::vector<double>{
        1, 2,  //
        3, 4,  //
        5, 6,  //
    };

    SECTION("special dim second") {
        auto rhs = tensor<2>{{na, nb}};

        rhs(0, 0) = 43;
        rhs(1, 0) = 61;

        rhs(0, 1) = -24;
        rhs(1, 1) = 136;

        rhs(0, 2) = 56;
        rhs(1, 2) = 11;

        auto buf = tensor<2>{{na, nb}};
        ads_solve(rhs, buf, dim_A, step_B);

        CHECK_THAT(to_vector(rhs), Approx(expected));
    }

    SECTION("special dim first") {
        auto rhs = tensor<2>{{nb, na}};

        rhs(0, 0) = 54;
        rhs(1, 0) = -12;
        rhs(2, 0) = 54;

        rhs(0, 1) = 71;
        rhs(1, 1) = 149;
        rhs(2, 1) = 19;

        auto buf = tensor<2>{{nb, na}};
        ads_solve(rhs, buf, step_B, dim_A);

        CHECK_THAT(to_vector(rhs), Approx(expected));
    }
}

TEST_CASE("ADS generalized 3D", "[ads]") {
    auto const na = 3;
    auto const nb = 2;
    auto const nc = 4;

    auto A = matrix{1, 1, na};
    auto B = matrix{1, 1, nb};

    A(0, 0) = 2;
    A(0, 1) = 3;
    A(1, 0) = 1;
    A(1, 1) = 2;
    A(1, 2) = -3;
    A(2, 1) = 2;
    A(2, 2) = 2;

    B(0, 0) = 1;
    B(0, 1) = 2;
    B(1, 0) = 3;
    B(1, 1) = 4;

    auto C11 = make_factored_matrix(1, 1, nc, {{1, 2}, {2, 3, 1}, {1, -4, 0}, {1, 3}});
    auto C21 = make_factored_matrix(1, 1, nc, {{1, 3}, {2, -2, 1}, {3, -2, 1}, {1, 1}});
    auto C31 = make_factored_matrix(1, 1, nc, {{3, 0}, {1, 2, 1}, {2, -1, 3}, {3, 1}});
    auto C12 = make_factored_matrix(1, 1, nc, {{2, 3}, {-1, 3, 2}, {2, 4, 1}, {-2, 1}});
    auto C22 = make_factored_matrix(1, 1, nc, {{-1, -2}, {3, 2, 1}, {2, -1, 1}, {1, 1}});
    auto C32 = make_factored_matrix(1, 1, nc, {{2, 1}, {-1, 3, 2}, {4, 2, 1}, {-3, 2}});

    auto ctx_A = ads::lin::solver_ctx{A};
    ads::lin::factorize(A, ctx_A);

    auto ctx_B = ads::lin::solver_ctx{B};
    ads::lin::factorize(B, ctx_B);

    auto dim_A = ads::dim_data{A, ctx_A};
    auto dim_B = ads::dim_data{B, ctx_B};

    auto const step_C = [&](auto& rhs) {
        auto* data = rhs.data();
        auto const rhs_size = rhs.size(0);

        solve_with_factorized(C11.mat, data + 0 * rhs_size, C11.ctx, 1);
        solve_with_factorized(C21.mat, data + 1 * rhs_size, C21.ctx, 1);
        solve_with_factorized(C31.mat, data + 2 * rhs_size, C31.ctx, 1);
        solve_with_factorized(C12.mat, data + 3 * rhs_size, C12.ctx, 1);
        solve_with_factorized(C22.mat, data + 4 * rhs_size, C22.ctx, 1);
        solve_with_factorized(C32.mat, data + 5 * rhs_size, C32.ctx, 1);
    };

    auto const expected = std::vector<double>{
        1,  2,  3,  //
        4,  5,  6,  //

        7,  8,  9,   //
        10, 11, 12,  //

        13, 14, 15,  //
        16, 17, 18,  //

        19, 20, 21,  //
        22, 23, 24,  //
    };

    SECTION("special dim first") {
        auto rhs = tensor<3>{{nc, na, nb}};

        rhs(0, 0, 0) = 543;
        rhs(1, 0, 0) = 1101;
        rhs(2, 0, 0) = -618;
        rhs(3, 0, 0) = 849;

        rhs(0, 1, 0) = -192;
        rhs(1, 1, 0) = -48;
        rhs(2, 1, 0) = -96;
        rhs(3, 1, 0) = -96;

        rhs(0, 2, 0) = 540;
        rhs(1, 2, 0) = 768;
        rhs(2, 2, 0) = 828;
        rhs(3, 2, 0) = 828;

        rhs(0, 0, 1) = 1900;
        rhs(1, 0, 1) = 1681;
        rhs(2, 0, 1) = 2968;
        rhs(3, 0, 1) = -394;

        rhs(0, 1, 1) = 336;
        rhs(1, 1, 1) = -672;
        rhs(2, 1, 1) = -224;
        rhs(3, 1, 1) = -224;

        rhs(0, 2, 1) = 1192;
        rhs(1, 2, 1) = 1748;
        rhs(2, 2, 1) = 3024;
        rhs(3, 2, 1) = -388;

        auto buf = tensor<3>{rhs.sizes()};
        ads_solve(rhs, buf, step_C, dim_A, dim_B);

        CHECK_THAT(to_vector(rhs), Approx(expected));
    }

    SECTION("special dim second") {
        auto rhs = tensor<3>{{nb, nc, na}};

        rhs(0, 0, 0) = 351;
        rhs(1, 0, 0) = 1325;

        rhs(0, 1, 0) = 732;
        rhs(1, 1, 0) = 1382;

        rhs(0, 2, 0) = -501;
        rhs(1, 2, 0) = 2471;

        rhs(0, 3, 0) = 718;
        rhs(1, 3, 0) = -293;

        rhs(0, 0, 1) = -384;
        rhs(1, 0, 1) = 672;

        rhs(0, 1, 1) = -96;
        rhs(1, 1, 1) = -1344;

        rhs(0, 2, 1) = -192;
        rhs(1, 2, 1) = -448;

        rhs(0, 3, 1) = -192;
        rhs(1, 3, 1) = -448;

        rhs(0, 0, 2) = 492;
        rhs(1, 0, 2) = 1196;

        rhs(0, 1, 2) = 752;
        rhs(1, 1, 2) = 1912;

        rhs(0, 2, 2) = 872;
        rhs(1, 2, 2) = 3276;

        rhs(0, 3, 2) = 872;
        rhs(1, 3, 2) = -380;

        auto buf = tensor<3>{rhs.sizes()};
        ads_solve(rhs, buf, dim_B, step_C, dim_A);

        CHECK_THAT(to_vector(rhs), Approx(expected));
    }

    SECTION("special dim third") {
        auto rhs = tensor<3>{{na, nb, nc}};

        rhs(0, 0, 0) = 342;
        rhs(1, 0, 0) = -48;
        rhs(2, 0, 0) = 162;

        rhs(0, 1, 0) = 1210;
        rhs(1, 1, 0) = 84;
        rhs(2, 1, 0) = 522;

        rhs(0, 0, 1) = 774;
        rhs(1, 0, 1) = -12;
        rhs(2, 0, 1) = 504;

        rhs(0, 1, 1) = 1934;
        rhs(1, 1, 1) = -168;
        rhs(2, 1, 1) = 1648;

        rhs(0, 0, 2) = -792;
        rhs(1, 0, 2) = -24;
        rhs(2, 0, 2) = 864;

        rhs(0, 1, 2) = 3542;
        rhs(1, 1, 2) = -56;
        rhs(2, 1, 2) = 2674;

        rhs(0, 0, 3) = 1206;
        rhs(1, 0, 3) = -24;
        rhs(2, 0, 3) = 864;

        rhs(0, 1, 3) = -326;
        rhs(1, 1, 3) = -56;
        rhs(2, 1, 3) = -118;

        auto buf = tensor<3>{rhs.sizes()};
        ads_solve(rhs, buf, dim_A, dim_B, step_C);

        CHECK_THAT(to_vector(rhs), Approx(expected));
    }
}
