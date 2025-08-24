// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "co2_sequestration_2d_trim.hpp"

#include <iostream>

#include <lyra/lyra.hpp>

int main(int argc, char* argv[]) {
    double mesh_x;
    double mesh_y;
    std::array<int, 2> resolution{1, 1};
    int iter;
    double timestep_size;

    bool show_help = false;
    bool verbose = false;

    double mu_w = 1.25;
    double mu_g = 0.2;
    double K = 1;
    double phi = 1;
    double rho_w = 2;
    double rho_g = 1;
    double g = 1;
    double qg_x = -1;
    double qg_y = -1;
    double qg_rate = 1e-6;
    std::string porosity_map = "none";
    std::string permeability_map = "none";

    auto cli =
        lyra::help(show_help)                                                                  //
        | lyra::arg(mesh_x, "mesh_size - x")("mesh resolution in the x direction").required()  //
        | lyra::arg(mesh_y, "mesh_size - y")("mesh resolution in the y direction").required()
        | lyra::arg(resolution, "resolution")(
            "resolution of the PINN grid / S&P data i/o").required()  //
        | lyra::arg(iter, "iter")("the current iteration").required()                 //
        | lyra::arg(timestep_size, "timestep_size")("size of the timestep").required()         //
        | lyra::opt(mu_w, "mu_w")["--mu_w"]("mu_w (parameter) - brine viscosity")              //
        | lyra::opt(mu_g, "mu_g")["--mu_g"]("mu_g (parameter) - gas viscosity")                //
        | lyra::opt(K, "K")["--K"]("K (parameter) - permeability tensor ")                     //
        | lyra::opt(phi, "phi")["--phi"]("phi (parameter) - porosity")                         //
        | lyra::opt(rho_w, "rho_w")["--rho_w"]("rho_w (parameter) - density of the brine")     //
        | lyra::opt(rho_g, "rho_g")["--rho_g"]("rho_g (parameter) - density of the gas")       //
        | lyra::opt(g, "g")["--g"]("g (parameter) - gravitational acceleration")               //
        | lyra::opt(qg_x, "qg_x")["--qg_x"](
            "qg_x (parameter) - gas injection location in x direction")  //
        | lyra::opt(qg_y, "qg_y")["--qg_y"](
            "qg_y (parameter) - gas injection location in y direction")                           //
        | lyra::opt(qg_rate, "qg_rate")["--qg_rate"]("qg_rate (parameter) - gas injection rate")  //
        | lyra::opt(porosity_map, "porosity_map")["--porosity_map"]("porosity map file path")     //
        | lyra::opt(permeability_map, "permeability_map")["--permeability_map"](
            "permeability map file path")  //
        | lyra::opt(verbose)["--verbose"];

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

    // if the user does not specify gas injection location - set it to the middle of the mesh
    if (qg_x == -1) {
        qg_x = mesh_x / 2;
    }

    if (qg_y == -1) {
        qg_y = mesh_y / 2;
    }

    if (verbose) {
        std::cout << "Argument parsing is complete" << std::endl;
        // print the values:
        std::cout << "mesh_size x: " << mesh_x << std::endl;
        std::cout << "mesh_size y: " << mesh_y << std::endl;
        std::cout << "resolution: " << resolution[0] << " x " << resolution[1] << std::endl;
        std::cout << "iter: " << iter << std::endl;
        std::cout << "timestep_size: " << timestep_size << std::endl;
        std::cout << "mu_w: " << mu_w << std::endl;
        std::cout << "mu_g: " << mu_g << std::endl;
        std::cout << "K: " << K << std::endl;
        std::cout << "phi: " << phi << std::endl;
        std::cout << "rho_w: " << rho_w << std::endl;
        std::cout << "rho_g: " << rho_g << std::endl;
        std::cout << "g: " << g << std::endl;
        std::cout << "qg_x: " << qg_x << std::endl;
        std::cout << "qg_y: " << qg_y << std::endl;
        std::cout << "qg_rate: " << qg_rate << std::endl;
        std::cout << "porosity_map: " << porosity_map << std::endl;
        std::cout << "permeability_map: " << permeability_map << std::endl;
    }

    int n_elem_x = static_cast<int>(mesh_x);
    int n_elem_y = static_cast<int>(mesh_y);

    ads::dim_config dim_x{2, 2 * n_elem_x, 0, mesh_x};
    ads::dim_config dim_y{2, 2 * n_elem_y, 0, mesh_y};
    ads::timesteps_config steps{1, timestep_size};
    int ders = 1;  // order of derivative

    ads::config_2d c{dim_x, dim_y, steps, ders};
    auto sim = ads::problems::co2_sequestration_2d_trim(c, resolution, iter, mu_w, mu_g, K, phi, rho_w, rho_g, g,
                                                   qg_x, qg_y, qg_rate, porosity_map, permeability_map,
                                                   verbose);
    sim.run();
}
