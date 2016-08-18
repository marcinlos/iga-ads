#ifndef PROBLEMS_TUMOR_VASCULATURE_VASCULATURE_HPP_
#define PROBLEMS_TUMOR_VASCULATURE_VASCULATURE_HPP_

#include <vector>
#include <queue>
#include <set>
#include <random>
#include <algorithm>
#include "ads/util/function_value.hpp"
#include "problems/tumor/vasculature/defs.hpp"
#include "problems/tumor/vasculature/rasterizer.hpp"
#include "problems/tumor/vasculature/plot.hpp"



namespace ads {
namespace tumor {
namespace vasc {


struct config {
    double init_stability = 0.5;
    double degeneration = 0.05;

//    double t_ec_sprout = 0.5;
//    double t_ec_migr = 100;

    double t_ec_sprout = 10;
    double t_ec_migr = 2;


    // my values
    double segment_length = 0.03;
    double t_ec_collapse = 10;
    double c_min = 0.2;
};


class vasculature {

public:
    struct segment;
    struct node;

    using segment_ptr = segment*;
    using node_ptr = node*;

    struct node {
        vector position;
        std::vector<segment_ptr> segments;

        node(vector p)
        : position(p)
        { }
    };

    struct segment {
        node_ptr begin, end;
        double stability;
    };

private:
    static constexpr std::size_t N = 500;
    using array = lin::tensor<double, Dim>;

    config cfg;
    array veins{{ N + 1, N + 1 }}, oxygen{{ N + 1, N + 1 }};
    std::vector<node_ptr> roots;

    std::set<node_ptr> nodes;
    std::set<segment_ptr> segments;

    std::mt19937 rng;

public:

    vasculature(std::vector<node_ptr> roots)
    : roots{ std::move(roots) }
    {
        std::queue<node_ptr> q;
        for (node_ptr node : this->roots) {
            q.push(node);
        }
        while (! q.empty()) {
            node_ptr n = q.front();
            q.pop();

            auto res = nodes.insert(n);
            if (res.second) {
                for (segment_ptr s : n->segments) {
                    segments.insert(s);
                    node_ptr other = s->begin != n ? s->begin : s->end;
                    q.push(other);
                }
            }
        }
    }

    void plot_veins(const std::string& file) const {
        plot(file, veins);
    }

    void plot_oxygen(const std::string& file) const {
        plot(file, oxygen);
    }

private:

    segment_ptr connect(node_ptr a, node_ptr b) {
        segment_ptr s = new segment{ a, b, cfg.init_stability };
        a->segments.push_back(s);
        b->segments.push_back(s);
        segments.insert(s);
        return s;
    }

    void remove(segment_ptr s) {
        remove_node_from_vector(s, s->begin->segments);
        remove_node_from_vector(s, s->end->segments);
        segments.erase(s);
        delete s;
    }

    void remove_node_from_vector(segment_ptr s, std::vector<segment_ptr>& v) {
        using std::begin;
        using std::end;
        auto it = std::find(begin(v), end(v), s);
        v.erase(it);
    }

    node_ptr make_node(vector p) {
        node_ptr n = new node{ p };
        nodes.insert(n);
        return n;
    }

public:
    void discretize() {
        zero(oxygen);
        zero(veins);
        for (segment_ptr s : segments) {
            draw_segment(s->begin->position, s->end->position, veins, 1);
        }
        blur(veins, oxygen, 3, 3);
    }

    const array& veins_grid() const {
        return veins;
    }

    const array& oxygen_grid() const {
        return oxygen;
    }

    double oxygen_level(double x, double y) const {
        int ix = static_cast<int>(x * N - 0.5);
        int iy = static_cast<int>(y * N - 0.5);
        return oxygen(ix, iy);
    }

    struct sprout {
        node_ptr tip;
        vector dir;
        double time;
    };

    std::vector<sprout> sprouts;
    std::set<node_ptr> sprout_nodes;

    using value_type = function_value_2d;

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


    template <typename Tumor, typename TAF>
    void update(Tumor&& tumor, TAF&& taf, double dt) {
        auto nodes_copy = nodes;
        for (node_ptr n : nodes_copy) {
            vector p = n->position;
            value_type c = taf(p.x, p.y);

            if (c.val > cfg.c_min) {
                if (flip_coin(dt / cfg.t_ec_sprout)) {
                    vector dir = normalized(grad(c));
                    vector end = p + cfg.segment_length * dir;
                    if (inside_domain(end)) {
                        node_ptr tip = make_node(end);
                        connect(n, tip);
                        sprouts.push_back( { tip, dir, 0 });
                    }
                }
            }
        }

        for (auto it = begin(sprouts); it != end(sprouts); ) {
            sprout& s = *it;
            node_ptr tip = s.tip;

            bool removed = false;
            if (flip_coin(dt / cfg.t_ec_migr)) {
                vector p = tip->position;
                value_type c = taf(p.x, p.y);
                vector dir = normalized(grad(c));
                vector end = p + cfg.segment_length * dir;
                if (inside_domain(end)) {
                    node_ptr new_tip = make_node(end);
                    connect(tip, new_tip);
                    s.tip = new_tip;

                    auto neighbor = find_neighbor(new_tip, tip, cfg.segment_length);
                    if (neighbor != nullptr) {
                        connect(neighbor, new_tip);
                        removed = true;
                        it = sprouts.erase(it);
                    }
                }
            }
            if (! removed)
                ++ it;
        }

        for (segment_ptr s : segments) {
            vector c = center(s);
            double b = tumor(c.x, c.y);
            if (b > 0.5) {
                s->stability -= b;
            }
            if (s->stability <= 0) {
                if (flip_coin(dt / cfg.t_ec_collapse)) {
//                    node_ptr begin = s->begin;
//                    node_ptr end = s->end;
                    remove(s);
                }
            }
        }
    }

    bool inside_domain(vector v) const {
        return 0 <= v.x && v.x <= 1 && 0 <= v.y && v.y <= 1;
    }

    vector center(segment_ptr s) {
        return 0.5 * (s->begin->position + s->end->position);
    }

    vector grad(value_type v) {
        return { v.dx, v.dy };
    }

    bool flip_coin(double p) {
        return rand(0, 1) < p;
    }

    double rand(double a, double b) {
        std::uniform_real_distribution<> dist{ a, b };
        return dist(rng);
    }

    void blur(array& src, array& dst, int r, double scale) const {
        double w = scale / ((2 * r + 1) * (2 * r + 1));
        for (int i = r; i <= N - r; ++ i) {
            for (int j = r; j <= N - r; ++ j) {
                double v = 0;
                for (int p = -r; p <= r; ++ p) {
                    for (int q = -r; q <= r; ++ q) {
                        v += src(i + p, j + q);
                    }
                }
                v *= w;
                dst(i, j) = v;
            }
        }
    }
};


}
}
}



#endif /* PROBLEMS_TUMOR_VASCULATURE_VASCULATURE_HPP_ */
