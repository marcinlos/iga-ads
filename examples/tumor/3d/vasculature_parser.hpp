// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_3D_VASCULATURE_PARSER_HPP
#define TUMOR_3D_VASCULATURE_PARSER_HPP

#include <iostream>
#include <tuple>
#include <vector>

#include "vasculature.hpp"


namespace tumor {

    void normalize_positions(std::vector<vessels::point_type>& points) {
        auto inf = std::numeric_limits<double>::infinity();
        double xmin = inf, xmax = -inf;
        double ymin = inf, ymax = -inf;
        double zmin = inf, zmax = -inf;
        for (const auto& p : points) {
            xmin = std::min(xmin, p.x);
            xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y);
            ymax = std::max(ymax, p.y);
            zmin = std::min(zmin, p.z);
            zmax = std::max(zmax, p.z);
        }
        for (auto& p : points) {
            p.x = (p.x - xmin) / (xmax - xmin);
            p.y = (p.y - ymin) / (ymax - ymin);
            p.z = (p.z - zmin) / (zmax - zmin);
        }
    }

vessels parse_vessels(std::istream& is) {
    vessels vs;

    constexpr auto ALL = std::numeric_limits<std::streamsize>::max();
    auto skip_lines = [&](int n = 1) { for (int i = 0; i < n; ++ i) is.ignore(ALL, '\n'); };
    auto skip_word = [&]{ is.ignore(ALL, ' '); };

    skip_lines(4);

    skip_word();
    int joint_count;
    is >> joint_count;
    skip_lines();

    using node_ptr = vessels::node_ptr;
    using point = vessels::point_type;

    std::vector<point> points;
    points.reserve(joint_count);

    for (int i = 0; i < joint_count; ++ i) {
        double x, y, z;
        is >> x >> y >> z;
        points.push_back(point{ x, y, z });

    }
    normalize_positions(points);

    std::vector<node_ptr> nodes;
    nodes.reserve(joint_count);

    for (const auto& pos : points) {
        auto node = vs.make_node(pos);
        nodes.push_back(node);
    }

    skip_word();
    int line_count;
    is >> line_count;
    skip_lines();

    std::vector<std::tuple<int, int>> lines;
    lines.reserve(line_count);

    for (int i = 0; i < line_count; ++ i) {
        skip_word();
        int j, k;
        is >> j >> k;
        lines.emplace_back(j, k);
    }

    skip_lines(4);

    std::vector<double> diameters;
    diameters.reserve(line_count);

    for (int i = 0; i < line_count; ++ i) {
        double d;
        is >> d;
        diameters.push_back(d);
    }

    skip_lines(2);

    std::vector<vessel_type> types;
    types.reserve(line_count);

    for (int i = 0; i < line_count; ++ i) {
        int type;
        is >> type;
        diameters.push_back(type);
    }

    for (int i = 0; i < line_count; ++ i) {
        int j = std::get<0>(lines[i]);
        int k = std::get<1>(lines[i]);
        double d = diameters[i];
        vs.connect(nodes[j], nodes[k], types[i], d);
    }

    return vs;
}

}

#endif // TUMOR_3D_VASCULATURE_PARSER_HPP
