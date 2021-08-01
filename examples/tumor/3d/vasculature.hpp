// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_3D_VASCULATURE_HPP
#define TUMOR_3D_VASCULATURE_HPP

#include <algorithm>
#include <random>
#include <set>
#include <vector>

#include "../vasculature/config.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/output/vtk.hpp"
#include "ads/output_manager.hpp"
#include "ads/util/function_value.hpp"
#include "ads/util/math/vec.hpp"
#include "rasterizer.hpp"

namespace tumor {

enum class vessel_type {
    vein,
    arthery,
    capillary,
    sprout,
};

class vessels {
public:
    using point_type = ads::math::vec<3>;
    using value_type = ads::function_value_3d;

    struct node;
    struct edge;

    using node_ptr = node*;
    using edge_ptr = edge*;

    struct node {
        point_type pos;
        std::vector<edge_ptr> edges;
    };

    struct edge {
        node_ptr src, dst;
        vessel_type type;
        double stability;
        double radius;
        double inside_tumor;
    };

private:
    std::vector<node_ptr> roots_;
    tumor::vasc::config cfg;
    std::mt19937 rng;

    std::set<node_ptr> nodes_;
    std::set<edge_ptr> edges_;

public:
    vessels() {
        auto xysize = 5000.0;
        auto zsize = 3000.0;
        auto zmin = 300 / zsize;
        auto zmax = 2400 / zsize;

        auto dh = 200;

        for (auto i = dh / 2; i < xysize; i += dh) {
            auto x = i / xysize;
            for (auto j = dh / 2; j < xysize; j += dh) {
                auto y = j / xysize;
                make_line({x, y, zmin}, {x, y, zmax}, 10);
            }
        }
    }

    void make_line(point_type a, point_type b, int steps) {
        auto* na = make_node(a);
        roots_.push_back(na);
        for (int i = 0; i < steps; ++i) {
            double t = static_cast<double>(i) / (steps - 1);
            auto* nb = make_node(a + t * (b - a));
            connect(na, nb, vessel_type::arthery, 1.0);
            na = nb;
        }
    }

    const std::set<edge_ptr>& edges() const { return edges_; }

    template <typename Tumor, typename TAF>
    void update(Tumor&& tumor, TAF&& taf, int iter, double dt) {
        if (iter % 240 == 0) {
            create_sprouts(tumor, taf, dt);
            for (edge_ptr s : edges_) {
                // Vessel collapse
                if (s->stability <= 0) {
                    if (flip_coin(dt / cfg.t_ec_collapse)) {
                        remove(s);
                    }
                }
            }
        }

        for (edge_ptr s : edges_) {
            auto p = center(s);
            double b = tumor(p.x, p.y, p.z);
            // Wall degradation
            if (b > 1) {
                s->stability -= cfg.degeneration * dt;
                s->inside_tumor += dt;
            }
            // Vessel dilatation
            if (s->inside_tumor > cfg.t_ec_switch && s->radius < cfg.r_max) {
                double c = taf(p.x, p.y, p.z).val;
                if (c > cfg.c_switch) {
                    s->radius += dt * cfg.dilatation;
                }
            }
        }
    }

    node_ptr make_node(point_type p) {
        auto* n = new node{p, {}};
        nodes_.insert(n);
        return n;
    }

    edge_ptr connect(node_ptr a, node_ptr b, vessel_type type, double radius) {
        auto* e = new edge{a, b, type, cfg.init_stability, radius, 0};
        edges_.insert(e);
        a->edges.push_back(e);
        b->edges.push_back(e);
        return e;
    }

private:
    template <typename Tumor, typename TAF>
    void create_sprouts(Tumor&&, TAF&& taf, double dt) {
        auto nodes_copy = nodes_;
        for (node_ptr n : nodes_copy) {
            auto p = n->pos;
            value_type c = taf(p.x, p.y, p.z);

            if (c.val > cfg.c_min) {
                // std::cout << "Maybe sprout" << std::endl;
                if (flip_coin(dt / cfg.t_ec_sprout)) {
                    auto dir = normalized(grad(c));
                    auto end = p + cfg.segment_length * dir;
                    if (inside_domain(end)) {
                        // std::cout << "Sprout indeed!" << std::endl;
                        node_ptr tip = make_node(end);
                        connect(n, tip, vessel_type::sprout, cfg.r_sprout);
                    }
                }
            }
        }
    }

    void remove(edge_ptr e) {
        remove_node(e, e->src->edges);
        remove_node(e, e->dst->edges);
        edges_.erase(e);
        delete e;
    }

    void remove_node(edge_ptr e, std::vector<edge_ptr>& v) {
        using std::begin;
        using std::end;
        auto it = std::find(begin(v), end(v), e);
        v.erase(it);
    }

    bool inside_domain(point_type v) const {
        return 0 <= v.x && v.x <= 1 && 0 <= v.y && v.y <= 1 && 0 <= v.z && v.z <= 1;
    }

    point_type center(edge_ptr s) { return 0.5 * (s->src->pos + s->dst->pos); }

    point_type grad(value_type v) { return {v.dx, v.dy, v.dz}; }

    bool flip_coin(double p) { return rand(0, 1) < p; }

    double rand(double a, double b) {
        std::uniform_real_distribution<> dist{a, b};
        return dist(rng);
    }
};

class vasculature {
private:
    using value_array = ads::lin::tensor<double, 3>;
    value_array src;
    int sx, sy, sz;
    vessels vs;
    rasterizer raster;

public:
    vasculature(int sx, int sy, int sz, vessels vs)
    : src{{sx, sy, sz}}
    , sx{sx}
    , sy{sy}
    , sz{sz}
    , vs{std::move(vs)} {
        rasterize();
    }

    double source(double x, double y, double z) const {
        using std::max;
        int ix = coord(x, sx);
        int iy = coord(y, sy);
        int iz = coord(z, sz);

        return src(ix, iy, iz);
    }

    void to_file(const std::string& name) const {
        namespace out = ads::output;
        out::vtk output{ads::DEFAULT_FMT};

        std::vector<double> px;
        px.reserve(sx);
        for (int i = 0; i < sx; ++i) {
            px.push_back(i);
        }

        std::vector<double> py;
        py.reserve(sy);
        for (int i = 0; i < sy; ++i) {
            py.push_back(i);
        }

        std::vector<double> pz;
        pz.reserve(sz);
        for (int i = 0; i < sz; ++i) {
            pz.push_back(i);
        }

        auto rx = out::from_container(px);
        auto ry = out::from_container(py);
        auto rz = out::from_container(pz);
        auto grid = out::make_grid(rx, ry, rz);
        std::ofstream os{name};
        output.print(os, grid, src);
    }

    template <typename Tumor, typename TAF>
    void update(Tumor&& tumor, TAF&& taf, int iter, double dt) {
        vs.update(tumor, taf, iter, dt);
        recompute();
    }

private:
    int coord(double t, int s) const { return static_cast<int>(t * s); }

    void recompute() {
        clear();
        rasterize();
    }

    void clear() { zero(src); }

    void draw(const vessels::edge& e) {
        auto a = e.src->pos;
        auto b = e.dst->pos;
        raster.draw(a, b, e.radius, src);
    }

    void rasterize() {
        for (const auto* e : vs.edges()) {
            draw(*e);
        }
    }
};

}  // namespace tumor

#endif  // TUMOR_3D_VASCULATURE_HPP
