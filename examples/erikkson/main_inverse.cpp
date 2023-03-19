// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <iostream>

#include <lyra/lyra.hpp>

#include "erikkson_inverse.hpp"
#include "inverse.hpp"

auto parse_args(int argc, char* argv[]) {
    struct {
        int n;
        int p_trial, c_trial;
        int p_test, c_test;
        double bx, by;
    } args{};

    bool show_help = false;

    auto const cli = lyra::help(show_help)                                        //
                   | lyra::arg(args.n, "N")("mesh resolution").required()         //
                   | lyra::arg(args.p_trial, "p")("trial p").required()           //
                   | lyra::arg(args.c_trial, "c")("trial continuity").required()  //
                   | lyra::arg(args.p_test, "P")("test p").required()             //
                   | lyra::arg(args.c_test, "C")("test continuity").required()    //
                   | lyra::arg(args.bx, "bx")("x advection").required()           //
                   | lyra::arg(args.by, "by")("y advection").required()           //
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

auto main(int argc, char* argv[]) -> int {
    auto const args = parse_args(argc, argv);

    auto const quad_order = std::max(args.p_trial, args.p_test) + 1;
    auto const ders = 1;

    auto const eps = 1e-2;
    auto const d = shishkin_const(args.n, eps);

    auto trial_basis_x =
        basis_from_points(make_points(args.n, false, true, d), args.p_trial, args.c_trial);
    auto trial_basis_y =
        basis_from_points(make_points(args.n, true, true, d), args.p_trial, args.c_trial);
    auto test_basis_x =
        basis_from_points(make_points(args.n, false, true, d), args.p_test, args.c_test);
    auto test_basis_y =
        basis_from_points(make_points(args.n, true, true, d), args.p_test, args.c_test);

    auto dtrial_x = ads::dimension{trial_basis_x, quad_order, ders};
    auto dtrial_y = ads::dimension{trial_basis_y, quad_order, ders};
    auto dtest_x = ads::dimension{test_basis_x, quad_order, ders};
    auto dtest_y = ads::dimension{test_basis_y, quad_order, ders};

    ads::erikkson_inverse{args.bx, args.by, dtrial_x, dtrial_y, dtest_x, dtest_y}.run();
}
