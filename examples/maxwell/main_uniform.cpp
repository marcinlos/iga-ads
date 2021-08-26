// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <iostream>

#include <lyra/lyra.hpp>

#include "maxwell_base.hpp"
#include "maxwell_uniform.hpp"

auto parse_args(int argc, char* argv[]) {
    struct {
        int n, p, c, step_count;
        double T;
    } args{};

    bool show_help = false;

    auto const* const desc =
        "Solver for non-stationary Maxwell equations with uniform material data\n"
        "using ADS";

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

    auto const dim = ads::dim_config{args.p, args.n};
    auto const cfg = ads::config_3d{dim, dim, dim, steps, 1};

    auto sim = maxwell_uniform{cfg};
    sim.run();
}
