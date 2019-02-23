#include "shishkin.hpp"



ads::bspline::basis create_basis(double a, double b, int p, int elements, int repeated_nodes, bool adapt, double d) {
    int points = elements + 1;
    int r = repeated_nodes + 1;
    int knot_size = 2 * (p + 1) + (points - 2) * r;
    ads::bspline::knot_vector knot(knot_size);

    for (int i = 0; i <= p; ++i) {
        knot[i] = a;
        knot[knot_size - i - 1] = b;
    }

    auto x0 = 0.5;
    auto y0 = 1 - d;

    for (int i = 1; i < points - 1; ++i) {
        auto t = ads::lerp(i, elements, 0.0, 1.0);

        auto s = adapt ? (t < x0 ? t / x0 * y0 : (t - x0) / (1 - x0) * (1 - y0) + y0) : t;
        for (int j = 0; j < r; ++ j) {
            knot[p + 1 + (i - 1) * r + j] = ads::lerp(s, a, b);
        }
    }

    return {std::move(knot), p};
}
