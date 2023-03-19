// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "vasculature.hpp"

namespace tumor::vasc {

vasculature::vasculature(std::vector<node_ptr> roots, const config& cfg)
: cfg{cfg}
, roots{std::move(roots)} {
    std::queue<node_ptr> q;
    for (node_ptr node : this->roots) {
        q.push(node);
    }
    while (!q.empty()) {
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

void vasculature::blur(array& src, array& dst, int r, double scale) const {
    double w = scale / ((2 * r + 1) * (2 * r + 1));
    for (unsigned i = r; i <= N - r; ++i) {
        for (unsigned j = r; j <= N - r; ++j) {
            double v = 0;
            for (int p = -r; p <= r; ++p) {
                for (int q = -r; q <= r; ++q) {
                    v += src(i + p, j + q);
                }
            }
            v *= w;
            dst(i, j) = v;
        }
    }
}

}  // namespace tumor::vasc
