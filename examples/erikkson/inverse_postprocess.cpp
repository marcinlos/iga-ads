// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include <lyra/lyra.hpp>

#include "ads/bspline/eval.hpp"
#include "ads/executor/galois.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/simulation/simulation_2d.hpp"
#include "ads/util/math/vec.hpp"
#include "inverse.hpp"

struct space_desc {
    int n;
    int p;
    int c;
    std::string path;
};

auto parse_args(int argc, char* argv[]) {
    struct {
        space_desc coarse;
        space_desc fine;
    } args{};

    bool show_help = false;

    auto const cli = lyra::help(show_help)                                                    //
                   | lyra::arg(args.coarse.path, "path1")("coarse solution path").required()  //
                   | lyra::arg(args.coarse.n, "N1")("coarse mesh resolution").required()      //
                   | lyra::arg(args.coarse.p, "p1")("coarse p").required()                    //
                   | lyra::arg(args.coarse.c, "c1")("coarse continuity").required()           //
                   | lyra::arg(args.fine.path, "path2")("fine solution path").required()      //
                   | lyra::arg(args.fine.n, "N2")("fine mesh resolution").required()          //
                   | lyra::arg(args.fine.p, "p2")("fine p").required()                        //
                   | lyra::arg(args.fine.c, "c2")("fine continuity").required()               //
        ;

    auto const result = cli.parse({argc, argv});

    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }

    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }

    return args;
}

auto read_solution(std::string const& path, ads::dimension const& x, ads::dimension const& y)
    -> ads::lin::tensor<double, 2> {
    auto data = ads::lin::tensor<double, 2>{{x.dofs(), y.dofs()}};
    auto input = std::ifstream{path};
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    for (int i = 0; i < data.size(); ++i) {
        int ix, iy;
        double val;
        input >> ix >> iy >> val;
        data(ix, iy) = val;
    }
    return data;
}

class postprocess : public ads::simulation_2d {
private:
    ads::dimension dcoarse_x;
    ads::dimension dcoarse_y;
    ads::dimension dfine_x;
    ads::dimension dfine_y;

    vector_type u_coarse, u_fine;

    ads::bspline::eval_ders_ctx ctx_coarse_x;
    ads::bspline::eval_ders_ctx ctx_coarse_y;
    ads::bspline::eval_ders_ctx ctx_fine_x;
    ads::bspline::eval_ders_ctx ctx_fine_y;

    ads::galois_executor executor{4};

public:
    postprocess(vector_type u_coarse,                                              //
                ads::dimension const& dcoarse_x, ads::dimension const& dcoarse_y,  //
                vector_type u_fine,                                                //
                ads::dimension const& dfine_x, ads::dimension const& dfine_y)
    : simulation_2d{dcoarse_x, dcoarse_y, ads::timesteps_config{1, 1.0}}
    , dcoarse_x{dcoarse_x}
    , dcoarse_y{dcoarse_y}
    , dfine_x{dfine_x}
    , dfine_y{dfine_y}
    , u_coarse{std::move(u_coarse)}
    , u_fine{std::move(u_fine)}
    , ctx_coarse_x{dcoarse_x.p, 1}
    , ctx_coarse_y{dcoarse_y.p, 1}
    , ctx_fine_x{dfine_x.p, 1}
    , ctx_fine_y{dfine_y.p, 1} { }

    auto eval_coarse(double x, double y) -> value_type {
        return ads::bspline::eval_ders(x, y, u_coarse, dcoarse_x.B, dcoarse_y.B, ctx_coarse_x,
                                       ctx_coarse_y);
    }

    auto eval_fine(double x, double y) -> value_type {
        return ads::bspline::eval_ders(x, y, u_fine, dfine_x.B, dfine_y.B, ctx_fine_x, ctx_fine_y);
    }

    auto before() -> void override {
        auto const diffL2 =
            errorL2(u_fine, dfine_x, dfine_y, [&](auto x) { return eval_coarse(x[0], x[1]); });
        auto const diffH1 =
            errorH1(u_fine, dfine_x, dfine_y, [&](auto x) { return eval_coarse(x[0], x[1]); });

        std::cout << diffL2 << " " << diffH1 << std::endl;

        auto const points = std::vector<ads::math::vec<2>>{
            {0.25, 0.25}, {0.25, 0.50}, {0.25, 0.75},  //
            {0.50, 0.25}, {0.50, 0.50}, {0.50, 0.75},  //
            {0.75, 0.25}, {0.75, 0.50}, {0.75, 0.75},  //
        };

        for (auto const [x, y] : points) {
            const auto val = eval_fine(x, y).val;
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
};

auto main(int argc, char* argv[]) -> int {
    auto const args = parse_args(argc, argv);
    auto const quad_order = std::max(args.coarse.p, args.fine.p) + 1;
    auto const ders = 1;

    auto const eps = 1e-2;
    auto const d_coarse = shishkin_const(args.coarse.n, eps);
    auto const d_fine = shishkin_const(args.coarse.n, eps);

    auto coarse_basis_x = basis_from_points(make_points(args.coarse.n, false, true, d_coarse),
                                            args.coarse.p, args.coarse.c);
    auto coarse_basis_y = basis_from_points(make_points(args.coarse.n, true, true, d_coarse),
                                            args.coarse.p, args.coarse.c);

    auto fine_basis_x =
        basis_from_points(make_points(args.fine.n, false, true, d_fine), args.fine.p, args.fine.c);
    auto fine_basis_y =
        basis_from_points(make_points(args.fine.n, true, true, d_fine), args.fine.p, args.fine.c);

    auto dcoarse_x = ads::dimension{coarse_basis_x, quad_order, ders};
    auto dcoarse_y = ads::dimension{coarse_basis_y, quad_order, ders};

    auto dfine_x = ads::dimension{fine_basis_x, quad_order, ders};
    auto dfine_y = ads::dimension{fine_basis_y, quad_order, ders};

    try {
        auto const u_coarse = read_solution(args.coarse.path, dcoarse_x, dcoarse_y);
        auto const u_fine = read_solution(args.fine.path, dfine_x, dfine_y);
        postprocess(u_coarse, dcoarse_x, dcoarse_y, u_fine, dfine_x, dfine_y).run();
    } catch (std::exception const& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }
}
