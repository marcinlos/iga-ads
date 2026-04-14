// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "fire.hpp"

#include <lyra/lyra.hpp>

int main(int argc, char* argv[]) {
    int n = 200;
    int p = 2;
    int threads = 8;
    int timesteps = 10'000;
    double dt = 1e-3;
    int plot_every = 10;

    bool show_help = false;
    auto const cli  //
        = lyra::cli() | lyra::help(show_help)
        | lyra::opt(n, "mesh_size")["-n"]["--elems"]  //
          ("Number of elements in each direction")
        | lyra::opt(p, "order")["-p"]  //
          ("Polynomial order")
        | lyra::opt(threads, "threads")["--threads"]  //
          ("How many threads to use")
        | lyra::opt(timesteps, "steps")["--steps"]  //
          ("Number of timesteps")
        | lyra::opt(dt, "step_size")["--dt"]  //
          ("Timestep size")
        | lyra::opt(plot_every, "every")["--plot-every"]  //
          ("How often to save plot data")
        //
        ;

    auto const result = cli.parse({argc, argv});

    if (!result) {
        std::cerr << "Error: " << result.message() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }

    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }

    ads::dim_config dim{p, n, 0.0, 100.0};
    ads::timesteps_config steps{timesteps, dt};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::fire sim{c, threads, plot_every};
    sim.run();
}
