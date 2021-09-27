// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "inverse.hpp"

#include <utility>

#include "ads/util.hpp"

auto choose_nodes(bool left, bool right, double width) -> std::array<double, 4> {
    if (left && right) {
        return {0.0, width, 1.0 - width, 1.0};
    } else if (left) {
        return {0.0, width, 0.75, 1.0};
    } else if (right) {
        return {0.0, 0.25, 1.0 - width, 1.0};
    } else {
        return {0.0, 0.25, 0.75, 1.0};
    }
}

auto make_points(int elems, bool left, bool right, double width) -> std::vector<double> {
    auto const nodes = choose_nodes(left, right, width);
    auto const map = [&](double t) {
        if (t < 0.25) {
            auto const s = t / 0.25;
            return ads::lerp(s, nodes[0], nodes[1]);
        } else if (t < 0.75) {
            auto const s = (t - 0.25) / (0.75 - 0.25);
            return ads::lerp(s, nodes[1], nodes[2]);
        } else {
            auto const s = (t - 0.75) / (1.0 - 0.75);
            return ads::lerp(s, nodes[2], nodes[3]);
        }
    };

    auto points = std::vector<double>(elems + 1);
    for (int i = 0; i <= elems; ++i) {
        auto const t = static_cast<double>(i) / elems;
        points[i] = map(t);
    }
    return points;
}

auto basis_from_points(std::vector<double> const& points, int p, int c) -> ads::bspline::basis {
    auto const repeated_nodes = p - 1 - c;
    auto const elems = ads::narrow_cast<int>(points.size() - 1);
    int const r = repeated_nodes + 1;
    int const size = (elems - 1) * r + 2 * (p + 1);
    auto knot = ads::bspline::knot_vector(size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = points[0];
        knot[size - i - 1] = points[elems];
    }

    for (int i = 1; i < elems; ++i) {
        for (int j = 0; j < r; ++j) {
            knot[p + 1 + (i - 1) * r + j] = points[i];
        }
    }
    return {std::move(knot), p};
}
