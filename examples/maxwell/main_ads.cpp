// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <iostream>
#include <string>

#include <lyra/lyra.hpp>

#include "maxwell_ads.hpp"
#include "maxwell_base.hpp"

auto parse_args(int argc, char* argv[]) {
    struct {
        int n, p, c, step_count;
        double T;
        std::string data_file;
    } args{};

    bool show_help = false;

    auto const desc = "Solver for non-stationary Maxwell equations with non-uniform material\n"
                      "data using generalized ADS";

    auto const cli = lyra::help(show_help).description(desc)                                  //
                   | common_arg_parser(args)                                                  //
                   | lyra::arg(args.data_file, "file")("file with material data").required()  //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

// Example invocation:
// <prog> 8 2 1 10 0.1 mri.dat
int main(int argc, char* argv[]) {
    auto const args = parse_args(argc, argv);
    auto const dt = args.T / args.step_count;
    auto const steps = ads::timesteps_config{args.step_count, dt};

    auto const dim = ads::dim_config{args.p, args.n};
    auto const cfg = ads::config_3d{dim, dim, dim, steps, 1};

    auto sim = maxwell_ads{cfg, args.data_file};
    sim.run();
}
