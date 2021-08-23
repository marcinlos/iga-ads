// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "head_data.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <system_error>

#include <fmt/core.h>

auto read_density_data(std::istream& input) -> ads::lin::tensor<std::uint8_t, 3> {
    using byte = std::uint8_t;
    int nx;
    int ny;
    int nz;
    input >> nz >> ny >> nx;
    auto data = ads::lin::tensor<byte, 3>{{nx, ny, nz}};

    auto const check_range = [nx, ny, nz](int ix, int iy, int iz, int val) {
        if (ix < 0 || ix >= nx)
            throw std::range_error{fmt::format("x index {} not in [0, {}) range", ix, nx)};
        if (iy < 0 || iy >= ny)
            throw std::range_error{fmt::format("y index {} not in [0, {}) range", iy, ny)};
        if (iz < 0 || iz >= nz)
            throw std::range_error{fmt::format("z index {} not in [0, {}) range", iz, nz)};
        if (val < 0 || val > 255)
            throw std::range_error{fmt::format("value {} not in [0, 255] range", val)};
    };

    for (int i = 0; i < data.size(); ++i) {
        int ix;
        int iy;
        int iz;
        int val;
        input >> iz >> iy >> ix >> val;
        check_range(ix, iy, iz, val);
        data(ix, iy, iz) = static_cast<byte>(val);
    }
    return data;
}

auto read_density_data(std::string_view path) -> ads::lin::tensor<std::uint8_t, 3> {
    std::ifstream input;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        input.open(std::string{path});
        return read_density_data(input);
    } catch (std::system_error const& e) {
        std::cerr << "Error reading " << path << ":\n" << e.code().message() << std::endl;
        std::cerr << "Make sure " << path << " is in the current working directory." << std::endl;
        std::exit(1);
    } catch (std::exception const& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }
}
