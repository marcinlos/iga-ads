#include "ads/bspline/bspline.hpp"
#include "ads/util.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

namespace ads {
namespace bspline {

basis create_basis(double a, double b, int p, int elements) {
    int points = elements + 1;
    int knot_size = points + 2 * p; // clamped B-spline
    knot_vector knot(knot_size);

    for (int i = 0; i < p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }
    for (int i = 0; i < points; ++i) {
        knot[i + p] = lerp(i, elements, a, b);
    }

    return {std::move(knot), p};
}

basis create_basis(double a, double b, int p, int elements, int repeated_nodes) {
    int points = elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }

    for (int i = 1; i < points - 1; ++i) {
        for (int j = 0; j < r; ++ j) {
            knot[p + 1 + (i - 1) * r + j] = lerp(i, elements, a, b);
        }
    }
    return {std::move(knot), p};
}


basis create_basis_C0(double a, double b, int p, int elements) {
    int points = elements + 1;
    int knot_size = p * points + 2; // clamped B-spline with C^0 separators
    knot_vector knot(knot_size);

    knot[0] = a;
    knot[knot_size - 1] = b;

    for (int i = 0; i < points; ++i) {
        for (int j = 0; j < p; ++ j) {
            knot[1 + i * p + j] = lerp(i, elements, a, b);
        }
    }
    return {std::move(knot), p};
}

int find_span(double x, const basis& b) {
    int low = b.begin_idx();
    int high = b.end_idx();

    if (x >= b[high]) {
        return high - 1;
    } else if (x <= b[low]) {
        return low;
    }

    int idx = (low + high) / 2;
    while (x < b[idx] || x >= b[idx + 1]) {
        if (x < b[idx]) {
            high = idx;
        } else {
            low = idx;
        }
        idx = (low + high) / 2;
    }
    return idx;
}


void eval_basis(int i, double x, const basis& b, double* out, basis_eval_ctx& ctx) {
    auto left = ctx.left();
    auto right = ctx.right();

    out[0] = 1;
    for (int j = 1; j <= b.degree; ++j) {
        left[j] = x - b.knot[i + 1 - j];
        right[j] = b.knot[i + j] - x;
        double saved = 0;

        for (int r = 0; r < j; ++r) {
            double tmp = out[r] / (right[r + 1] + left[j - r]);
            out[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        out[j] = saved;
    }
}


void eval_basis_with_derivatives(int i, double x, const basis& b, double** out, int der, basis_eval_ctx& ctx) {
    auto& ndu = ctx.ndu;
    auto& a = ctx.a;
    auto left = ctx.left();
    auto right = ctx.right();

    ndu(0, 0) = 1;
    for (int j = 1; j <= b.degree; ++ j) {
        left[j] = x - b.knot[i + 1 - j];
        right[j] = b.knot[i + j] - x;
        double saved = 0;

        for (int r = 0; r < j; ++ r) {
            ndu(j, r) = right[r + 1] + left[j - r];
            double tmp = ndu(r, j - 1) / ndu(j, r);

            ndu(r, j) = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        ndu(j, j) = saved;
    }
    for (int j = 0; j <= b.degree; ++ j) {
        out[0][j] = ndu(j, b.degree);
    }
    for (int r = 0; r <= b.degree; ++ r) {
        int s1 = 0, s2 = 1;
        a(0, 0) = 1;
        for (int k = 1; k <= der; ++ k) {
            double d = 0;
            int rk = r - k, pk = b.degree - k;
            if (r >= k) {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? k - 1 : b.degree - r;
            for (int j = j1; j <= j2; ++ j) {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }
            if (r <= pk) {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }
            out[k][r] = d;
            std::swap(s1, s2);
        }
    }
    int r = b.degree;
    for (int k = 1; k <= der; ++ k) {
        for (int j = 0; j <= b.degree; ++ j) {
            out[k][j] *= r;
        }
        r *= (b.degree - k);
    }
}

std::vector<int> first_nonzero_dofs(const basis& b) {
    std::vector<int> dofs(b.elements());
    int p = b.degree;
    int e = 0;
    for (std::size_t i = p; i + 1 < b.knot_size() - p; ++ i) {
        if (b.knot[i] != b.knot[i + 1]) {
            dofs[e ++] = i - p;
        }
    }
    return dofs;
}

std::vector<std::pair<int, int>> elements_supporting_dofs(const basis& b) {
    std::vector<std::pair<int, int>> ranges(b.dofs());
    int p = b.degree;
    int e = 0;

    for (int i = 0; i < b.dofs(); ++ i) {
        auto ee = e - 1;
        for (int j = 0; j < p + 1; ++ j) {
            if (b.knot[i + j] != b.knot[i + j + 1]) {
                ++ ee;
            }
        }
        ranges[i].first = e;
        ranges[i].second = ee;
        if (b.knot[i] != b.knot[i + 1]) {
            ++ e;
        }
    }

    return ranges;

}

}
}
