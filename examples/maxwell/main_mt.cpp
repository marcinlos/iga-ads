// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <iostream>

#include <lyra/lyra.hpp>

#include "ads/experimental/all.hpp"
#include "maxwell_absorbing_bc.hpp"
#include "maxwell_base.hpp"
#include "mt_problem.hpp"
#include "mt_source.hpp"

auto parse_args(int argc, char* argv[]) {
    struct {
        int n, p, c, step_count;
        double T;
    } args{};

    bool show_help = false;

    auto const* const desc =
        "Solver for non-stationary Maxwell equations for magnetotelluric problem \n"
        "with absorbing BC using ADS";

    auto const cli = lyra::help(show_help).description(desc)  //
                   | common_arg_parser(args)                  //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

// Example invocation:
// <prog> 8 2 1 10 0.1
int main(int argc, char* argv[]) {
    auto const args = parse_args(argc, argv);
    auto const dt = args.T / args.step_count;
    auto const steps = ads::timesteps_config{args.step_count, dt};

    auto const n = args.n;
    auto const nx = n * 3;
    auto const ny = n * 3;
    auto const nz = n;
    auto const end_x = 450e3;
    auto const end_y = 450e3;
    auto const end_z = 150e3;

    auto const rep = args.p - args.c - 1;
    auto const dim_x = ads::dim_config{args.p, nx, 0.0, end_x, args.p + 1, rep};
    auto const dim_y = ads::dim_config{args.p, ny, 0.0, end_y, args.p + 1, rep};
    auto const dim_z = ads::dim_config{args.p, nz, 0.0, end_z, args.p + 1, rep};
    auto const cfg = ads::config_3d{dim_x, dim_y, dim_z, steps, 1};

    // new interface
    auto xs = ads::evenly_spaced(0.0, end_x, nx);
    auto ys = ads::evenly_spaced(0.0, end_y, ny);
    auto zs = ads::evenly_spaced(0.0, end_z, nz);

    auto bx = ads::make_bspline_basis(xs, args.p, args.c);
    auto by = ads::make_bspline_basis(ys, args.p, args.c);
    auto bz = ads::make_bspline_basis(zs, args.p, args.c);

    auto mesh = ads::regular_mesh3{xs, ys, zs};
    auto quad = ads::quadrature3{&mesh, std::max(args.p + 1, 2)};

    using Problem = mt_problem;
    using Source = mt_source;
    auto problem = Problem{};
    auto source = Source{};

    auto space = ads::space3{&mesh, bx, by, bz};

    auto sim = maxwell_absorbing_bc<Problem, Source>{
        cfg, mesh, quad, space, std::move(problem), std::move(source)};
    sim.run();
}
