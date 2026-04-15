// SPDX-FileCopyrightText: 2026 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "fire_implicit.hpp"

#include <cstdlib>
#include <optional>

#include <lyra/lyra.hpp>

#include "params.hpp"

namespace {

ads::bspline::basis create_basis(double a, double b, int p, int elements, int repeated_nodes) {
    int points = elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    ads::bspline::knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }
    for (int i = 1; i < points - 1; ++i) {
        auto t = ads::lerp(i, elements, 0.0, 1.0);

        for (int j = 0; j < r; ++j) {
            knot[p + 1 + (i - 1) * r + j] = ads::lerp(t, a, b);
        }
    }

    return {std::move(knot), p};
}

auto parse_scheme(std::string_view name) -> std::optional<scheme> {
    if (name == "FE") {
        return scheme::FE;
    } else if (name == "BE") {
        return scheme::BE;
    } else if (name == "CN") {
        return scheme::CN;
    } else if (name == "PR") {
        return scheme::peaceman_rachford;
    } else if (name == "strang-BE") {
        return scheme::strang_BE;
    } else if (name == "strang-CN") {
        return scheme::strang_CN;
    } else {
        return std::nullopt;
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    int n = 200;
    int p = 2;
    std::string scheme_name;
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
        | lyra::opt(scheme_name, "scheme")["--scheme"]  //
          ("Time discretization scheme")
              .required()
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

    auto const scheme = parse_scheme(scheme_name);
    if (!scheme) {
        std::cerr << "Unknown scheme: " << scheme_name << std::endl;
        std::exit(1);
    }

    auto const steps = ads::timesteps_config{timesteps, dt};

    auto const ders = 2;
    auto const quad = p + 1;

    auto const S = 100.0;
    auto trial_basis_x = create_basis(0, S, p, n, 0);
    auto dtrial_x = ads::dimension{trial_basis_x, quad, ders};

    auto trial_basis_y = create_basis(0, S, p, n, 0);
    auto dtrial_y = ads::dimension{trial_basis_y, quad, ders};

    auto test_basis_x = create_basis(0, S, p, n, 0);
    auto dtest_x = ads::dimension{test_basis_x, quad, ders};

    auto test_basis_y = create_basis(0, S, p, n, 0);
    auto dtest_y = ads::dimension{test_basis_y, quad, ders};

    auto const params = fire_params{};
    auto sim = fire_implicit{
        dtrial_x, dtrial_y, dtest_x, dtest_y, params, threads, *scheme, steps,
    };
    sim.run();
}
