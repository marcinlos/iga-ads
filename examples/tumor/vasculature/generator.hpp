// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_VASCULATURE_GENERATOR_HPP
#define TUMOR_VASCULATURE_GENERATOR_HPP

#include <cmath>
#include <random>

#include "ads/util.hpp"
#include "defs.hpp"
#include "vasculature.hpp"

namespace tumor::vasc {

class random_vasculature {
private:
    std::mt19937 rng;
    config cfg;

    using node = vasculature::node;
    using node_ptr = vasculature::node_ptr;
    using segment = vasculature::segment;
    using segment_ptr = vasculature::segment_ptr;

public:
    random_vasculature(config cfg, int seed = 0)
    : rng{seed}
    , cfg{cfg} { }

    vasculature operator()() {
        std::vector<vasculature::node_ptr> nodes;
        int n = 5;
        for (int i = 0; i <= n; ++i) {
            vector pos = {0.1,
                          ads::lerp(i, n, 0.1, 0.9)};  // random_point({0.05, 0.1}, {0.1, 0.9});
            vector bias = {1, 0};
            node_ptr root = grow_tree(pos, bias);
            nodes.push_back(root);
        }
        for (int i = 0; i <= n; ++i) {
            vector pos = {0.9, ads::lerp(i, n, 0.1, 0.9)};  // random_point({0.9, 0.1}, {0.85,
                                                            // 0.9});
            vector bias = {-1, 0};
            node_ptr root = grow_tree(pos, bias);
            nodes.push_back(root);
        }
        return vasculature{nodes, cfg};
    }

private:
    std::vector<node_ptr> nodes;

    node_ptr grow_tree(vector p, vector bias) {
        node_ptr root = new node{p};
        nodes.push_back(root);

        double bias_strength = 10;
        vector dir = normalized(random_dir() + bias_strength * bias);
        grow_from(root, dir, 0.03, 60);
        return root;
    }

    node_ptr find_neighbor(node_ptr node, node_ptr prev, double dist) {
        for (node_ptr n : nodes) {
            if (n != node && n != prev) {
                if (norm(node->position - n->position) < dist) {
                    return n;
                }
            }
        }
        return nullptr;
    }

    void grow_from(node_ptr root, vector dir, double segment_length, double expected_length) {
        node_ptr n = root;

        int segments = 0;
        int min_segments = static_cast<int>(expected_length / 2);
        int max_segments = static_cast<int>(expected_length * 2);

        while (true) {
            double length = rand(0.5 * segment_length, 1.1 * segment_length);
            double max_angle = 0.07 * M_PI;
            double dphi = rand(-max_angle, max_angle);

            dir = rotate(dir, dphi);
            vector s = length * dir;
            vector end = n->position + s;
            if (!inside_domain(end)) {
                break;
            }
            node_ptr prev = n;
            n = new node{end};
            connect(prev, n, 1);

            nodes.push_back(n);
            node_ptr neighor = find_neighbor(n, prev, segment_length * 0.5);
            if (neighor != nullptr) {
                connect(neighor, n, 1);
                break;
            }

            if (flip_coin(0.05)) {
                // small branch
                double dev_back = 0.1 * M_PI;
                double dev_fwd = 0.3 * M_PI;
                double angle = rand_sign() * (M_PI / 2 + rand(dev_back, -dev_fwd));
                grow_from(prev, rotate(dir, angle), segment_length, expected_length * 0.5);
            }

            if (flip_coin(0.06)) {
                // split
                double min_angle = 0.2 * M_PI;
                double max_angle = 0.25 * M_PI;
                double angle1 = rand(min_angle, max_angle);
                double angle2 = rand(min_angle, max_angle);
                grow_from(n, rotate(dir, angle1), segment_length, expected_length * 0.6);
                grow_from(n, rotate(dir, -angle2), segment_length, expected_length * 0.6);
                break;
            }

            ++segments;
            if (segments > min_segments) {
                double p = 1 / expected_length;
                if (flip_coin(p)) {
                    break;
                }
            }
            if (segments >= max_segments) {
                break;
            }
        }
    }

    segment_ptr connect(node_ptr a, node_ptr b, double stability) {
        segment_ptr s = new segment{a, b, stability};
        a->segments.push_back(s);
        b->segments.push_back(s);
        return s;
    }

    bool inside_domain(vector v) const { return 0 <= v.x && v.x <= 1 && 0 <= v.y && v.y <= 1; }

    double rand(double a, double b) {
        std::uniform_real_distribution<> dist{a, b};
        return dist(rng);
    }

    int rand_sign() {
        std::uniform_int_distribution<> dist{0, 1};
        return 1 - 2 * dist(rng);
    }

    bool flip_coin(double p) { return rand(0, 1) < p; }

    vector random_dir() { return random_dir(0, 2 * M_PI); }

    vector random_dir(double phi1, double phi2) {
        double phi = rand(phi1, phi2) * M_PI;
        return {std::cos(phi), std::sin(phi)};
    }

    vector random_point(vector a, vector b) { return {rand(a.x, b.x), rand(a.y, b.y)}; }

    vector rotate(vector a, double phi) const {
        double cos_phi = std::cos(phi);
        double sin_phi = std::sin(phi);
        return {
            a.x * cos_phi - a.y * sin_phi,
            a.x * sin_phi + a.y * cos_phi,
        };
    }
};

}  // namespace tumor::vasc

#endif  // TUMOR_VASCULATURE_GENERATOR_HPP
