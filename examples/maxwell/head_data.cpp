// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "head_data.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <system_error>

auto read_density_data(std::istream& input) -> ads::lin::tensor<std::uint8_t, 3> {
    using byte = std::uint8_t;
    int nx;
    int ny;
    int nz;
    input >> nz >> ny >> nx;
    auto data = ads::lin::tensor<byte, 3>{{nx, ny, nz}};

    for (int i = 0; i < data.size(); ++i) {
        int ix;
        int iy;
        int iz;
        int val;
        input >> iz >> iy >> ix >> val;
        data(ix, iy, iz) = static_cast<byte>(val);
    }
    return data;
}

auto read_density_data(const std::string& path) -> ads::lin::tensor<std::uint8_t, 3> {
    std::ifstream input;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        input.open(path);
        return read_density_data(input);
    } catch (const std::system_error& e) {
        std::cerr << "Error reading " << path << ":\n" << e.code().message() << std::endl;
        std::cerr << "Make sure " << path << " is in the current working directory." << std::endl;
        std::exit(1);
    }
}
