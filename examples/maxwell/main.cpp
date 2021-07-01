// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <clara.hpp>

#include "maxwell_galerkin.hpp"


using namespace ads;
using namespace clara;

int main(int argc, char* argv[]) {
    int n;
    int p;
    int c;
    int step_count;

    bool help = false;
    auto cli = Help(help)
        | Arg(n, "N").required()
        | Arg(p, "p").required()
        | Arg(c, "c").required()
        | Arg(step_count, "steps").required();

    auto result = cli.parse(Args(argc, argv));

    if (! result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::exit(1);
    }

    if (help) {
        cli.writeToStream(std::cout);
        return 0;
    }

    if (argc < 5) {
        cli.writeToStream(std::cout);
        std::exit(1);
    }

    auto T = 1.0;
    auto dt = T / step_count;
    auto steps = timesteps_config{step_count, dt};

    auto dim = dim_config{p, n};
    auto cfg = config_3d{dim, dim, dim, steps, 1};

    auto sim = maxwell_galerkin{cfg};
    sim.run();
}
