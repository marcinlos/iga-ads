// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <iostream>

#include <clara.hpp>

#include "maxwell_galerkin.hpp"

auto parse_args(int argc, char* argv[]) {
    struct {
        int n, p, c, step_count;
        double T;
    } args;

    bool help = false;

    using clara::Help, clara::Arg;
    auto const cli = Help(help)                                //
                   | Arg(args.n, "N").required()               //
                   | Arg(args.p, "p").required()               //
                   | Arg(args.c, "c").required()               //
                   | Arg(args.step_count, "steps").required()  //
                   | Arg(args.T, "T").required()               //
        ;

    auto const result = cli.parse({argc, argv});

    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::exit(1);
    }
    if (help) {
        cli.writeToStream(std::cout);
        std::exit(0);
    }
    if (argc < 6) {
        cli.writeToStream(std::cout);
        std::exit(1);
    }

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

    auto sim = maxwell_galerkin{cfg};
    sim.run();
}
