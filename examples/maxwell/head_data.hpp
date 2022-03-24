// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_HEAD_DATA_HPP
#define MAXWELL_HEAD_DATA_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <iosfwd>
#include <string_view>
#include <utility>

#include "ads/lin/tensor/tensor.hpp"

class head_data {
public:
    using byte = std::uint8_t;
    using density_data = ads::lin::tensor<byte, 3>;
    using point_type = std::array<double, 3>;

private:
    enum class material { air, tissue, bone };

    static constexpr std::array<double, 3> eps_vals = {1.0, 45.8, 16.6};
    static constexpr std::array<double, 3> mu_vals = {1.0, 1.0, 1.0};
    static constexpr std::array<double, 3> sigma_vals = {1.0, 15.31, 4.8};

    density_data density_map_;

public:
    explicit head_data(density_data data)
    : density_map_{std::move(data)} { }

    auto eps(point_type x) const -> double { return eps_vals[index_at(x)]; }

    auto mu(point_type x) const -> double { return mu_vals[index_at(x)]; }

    auto sigma(point_type x) const -> double { return sigma_vals[index_at(x)]; }

    auto empty(point_type x) const -> bool { return material_at(x) == material::air; }

private:
    auto density(point_type x) const -> byte {
        const auto ix = discretize(x[0], density_map_.size(0));
        const auto iy = discretize(x[1], density_map_.size(1));
        const auto iz = discretize(x[2], density_map_.size(2));

        return density_map_(ix, iy, iz);
    }

    auto as_index(material m) const noexcept -> int { return static_cast<int>(m); }

    auto material_at(point_type x) const noexcept -> material { return as_material(density(x)); }

    auto index_at(point_type x) const noexcept -> int { return as_index(material_at(x)); }

    static auto discretize(double t, int n) noexcept -> int {
        return std::min(static_cast<int>(t * n), n - 1);
    }

    auto as_material(byte n) const noexcept -> material {
        if (n <= 1) {
            return material::air;
        } else if (n <= 240) {
            return material::tissue;
        } else {
            return material::bone;
        }
    }
};

auto read_density_data(std::istream& input) -> ads::lin::tensor<std::uint8_t, 3>;

auto read_density_data(std::string_view path) -> ads::lin::tensor<std::uint8_t, 3>;

#endif  // MAXWELL_HEAD_DATA_HPP
