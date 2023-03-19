// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "results.hpp"

#include <iostream>

auto print_result_info(maxwell_result_info const& res) -> void {
    std::cout << "  |E|     = " << res.E.norm.total().L2 << "  " << res.E.norm.total().H1 << '\n'
              << "  |rot E| = " << res.E.norm.rot << '\n'
              << "  |div E| = " << res.E.norm.div << '\n'
              << "  E err L2 = " << res.E.error.total().L2 << "  " << res.E.error.total().H1 << '\n'
              << "    rel L2 = " << res.E.relative_error().L2 * 100 << "%  "
              << res.E.relative_error().H1 * 100 << "%" << '\n'
              << "  E err rot = " << res.E.error.rot << '\n'
              << "  |H| = " << res.H.norm.total().L2 << "  " << res.H.norm.total().H1 << '\n'
              << "  |rot H| = " << res.H.norm.rot << '\n'
              << "  |div H| = " << res.H.norm.div << '\n'
              << "  H err L2 = " << res.H.error.total().L2 << "  " << res.H.error.total().H1 << '\n'
              << "    rel L2 = " << res.H.relative_error().L2 * 100 << "%  "
              << res.H.relative_error().H1 * 100 << "%" << '\n'
              << "  H err rot = " << res.H.error.rot << std::endl;
}
