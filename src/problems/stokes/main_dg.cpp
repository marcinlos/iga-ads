#include <cstdlib>
#include <iostream>
#include <string>

#include "ads/bspline/bspline.hpp"
#include "problems/stokes/stokes.hpp"
// #include "problems/stokes/stokes_conforming.hpp"
// #include "problems/stokes/stokes_constrained.hpp"
#include "problems/stokes/stokes_dg.hpp"




using namespace ads;

using index_type = std::array<int, 2>;
using point_type = std::array<double, 2>;
using value_type = function_value_2d;
using value_pair = std::array<value_type, 2>;

double dot(value_type a, point_type b) {
    return a.dx * b[0] + a.dy * b[1];
}

value_type average(const value_pair& v) {
    return 0.5 * (v[0] + v[1]);
}

value_type jump(const value_pair& v) {
    return v[0] - v[1];
}

value_type eval_basis_at(point_type p, index_type span, index_type dof, const dimension& x, const dimension& y) {
    bspline::eval_ders_ctx cx{x.p, 1};
    bspline::eval_ders_ctx cy{y.p, 1};

    double** bvx = cx.basis_vals();
    double** bvy = cy.basis_vals();

    eval_basis_with_derivatives(span[0], p[0], x.B, bvx, 1, cx);
    eval_basis_with_derivatives(span[1], p[1], y.B, bvy, 1, cy);

    int offsetx = span[0] - x.p;
    int offsety = span[1] - y.p;

    int ix = dof[0] - offsetx;
    int iy = dof[1] - offsety;

    if (ix < 0 || ix > x.p || iy < 0 || iy > y.p) return {};

    auto value = bvx[0][ix] * bvy[0][iy];
    auto dx    = bvx[1][ix] * bvy[0][iy];
    auto dy    = bvx[0][ix] * bvy[1][iy];

    return { value, dx, dy };
}

value_pair eval_basis_at_edge(point_type p, boundary orientation, index_type dof, const dimension& x, const dimension& y) {
    int spanx = bspline::find_span(p[0], x.B);
    int spany = bspline::find_span(p[1], y.B);
    auto span = index_type{spanx, spany};

    auto val1 = eval_basis_at(p, span, dof, x, y);
    auto val0 = value_type{};

    if (orientation == boundary::vertical) {
        int spanx0 = bspline::find_span(p[0] - 1e-10, x.B);
        val0 = eval_basis_at(p, index_type{spanx0, spany}, dof, x, y);
    } else {
        int spany0 = bspline::find_span(p[1] - 1e-10, y.B);
        val0 = eval_basis_at(p, index_type{spanx, spany0}, dof, x, y);
    }
    return {val0, val1};
}

bool supported_in_1d(int dof, int e, const dimension& x) {
    auto xrange = x.basis.element_ranges[dof];
    return e >= xrange.first && e <= xrange.second;
}

bool touches_point(int dof, int sx, const dimension& x) {
    auto xrange = x.basis.element_ranges[dof];
    return sx >= xrange.first && sx <= xrange.second + 1;
}

template <typename Form>
double integrate_vface(int sx, int ey, index_type i, index_type j,
                       const dimension& Ux, const dimension& Uy,
                       const dimension& Vx, const dimension& Vy,
                       Form&& form) {

    if (! touches_point(i[0], sx, Ux)   || ! touches_point(j[0], sx, Vx))   return 0;
    if (! supported_in_1d(i[1], ey, Uy) || ! supported_in_1d(j[1], ey, Vy)) return 0;

    double x0 = Ux.basis.points[sx];
    auto n = point_type{1, 0};

    double val = 0;
    for (int q = 0; q < Uy.basis.quad_order; ++ q) {
        double w = Uy.basis.w[q];
        double J = Uy.basis.J[ey];

        auto x = point_type{x0, Uy.basis.x[ey][q]};
        auto uu = eval_basis_at_edge(x, boundary::vertical, i, Ux, Uy);
        auto vv = eval_basis_at_edge(x, boundary::vertical, j, Vx, Vy);

        auto fuv = form(uu, vv, x, n);
        val += fuv * w * J;
    }
    return val;
}

template <typename Form>
double integrate_hface(int sy, int ex, index_type i, index_type j,
                       const dimension& Ux, const dimension& Uy,
                       const dimension& Vx, const dimension& Vy,
                       Form&& form) {

    if (! supported_in_1d(i[0], ex, Ux) || ! supported_in_1d(j[0], ex, Vx)) return 0;
    if (! touches_point(i[1], sy, Uy)   || ! touches_point(j[1], sy, Vy))   return 0;

    double y0 = Uy.basis.points[sy];
    auto n = point_type{0, 1};

    double val = 0;
    for (int q = 0; q < Ux.basis.quad_order; ++ q) {
        double w = Ux.basis.w[q];
        double J = Ux.basis.J[ex];

        auto x = point_type{Ux.basis.x[ex][q], y0};
        auto uu = eval_basis_at_edge(x, boundary::horizontal, i, Ux, Uy);
        auto vv = eval_basis_at_edge(x, boundary::horizontal, j, Vx, Vy);

        auto fuv = form(uu, vv, x, n);
        val += fuv * w * J;
    }
    return val;
}


int main2(int argc, char* argv[]) {

    auto n = 4;
    auto p = 1;
    auto basis = bspline::create_basis(0, 1, p, n, p);
    auto dim = dimension{ basis, p + 1, 1, 1 };

    // for (auto v : basis.knot) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;

    // for (auto v : basis.points) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;

    // auto x = 0.5;
    // int span = bspline::find_span(x, basis);
    // int span2 = span - (p + 1);

    // auto ctx = bspline::eval_ders_ctx{p, 1};
    // auto ctx2 = bspline::eval_ders_ctx{p, 1};

    // auto data = ctx.basis_vals();
    // auto data2 = ctx2.basis_vals();

    // double u[] = {0, 1, 2, 0, 0, 0, 0, 0};

    // bspline::eval_basis_with_derivatives(span, x, basis, data, 1, ctx);
    // bspline::eval_basis_with_derivatives(span2, x, basis, data2, 1, ctx2);

    // int offset = span - p;
    // int offset2 = span2 - p;

    // std::cout << "span: " << span << std::endl;
    // std::cout << "offset:  " << offset << std::endl;
    // std::cout << "offset2: " << offset2 << std::endl;

    // double val = 0;
    // double val2 = 0;

    // for (int i = 0; i < basis.dofs_per_element(); ++ i) {
    //     int ii = i + offset;
    //     int ii2 = i + offset2;
    //     val += u[ii] * data[1][i];
    //     val2 += u[ii2] * data2[1][i];
    // }

    // std::cout << "val = " << val << ", val2 = " << val2 << std::endl;

    // for (int i = 0; i < basis.dofs(); ++ i) {
    //     std::cout << "dof " << i << " touches 0.5? " << touches_point(i, 2, dim) << std::endl;
    // }

    std::cout << "integral: " << integrate_vface(3, 1, {5, 3}, {4, 2}, dim, dim, dim, dim,
                                                 [](auto uu, auto vv, auto x, auto n) {
                                                     return jump(uu).val * dot(average(vv), n);
                                                 }) << std::endl;
    return 0;
}


int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: stokes <N> <p_trial> <C_trial> <p_test> <C_test>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int subdivision = 1;

    int p_trial = std::atoi(argv[2]);
    int C_trial = std::atoi(argv[3]);
    int p_test = std::atoi(argv[4]);
    int C_test = std::atoi(argv[5]);

    int quad = std::max(p_trial, p_test) + 1 + 1; // to integrate velocity

    timesteps_config steps{ 1, 0 };
    int ders = 2;
    int rep_trial = p_trial - 1 - C_trial;
    int rep_test = p_test - 1 - C_test;

    // Pressure spaces
    auto trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_x = dimension{ trial_basis_x, quad, ders, subdivision };

    auto trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto dtrial_y = dimension{ trial_basis_y, quad, ders, subdivision };

    auto test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_x = dimension{ test_basis_x, quad, ders, 1 };

    auto test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto dtest_y = dimension{ test_basis_y, quad, ders, 1 };

    // Velocity spaces
    auto U1_trial_basis_x = bspline::create_basis(0, 1, p_trial + 1, n, rep_trial);
    auto U1_dtrial_x = dimension{ U1_trial_basis_x, quad, ders, subdivision };

    auto U1_trial_basis_y = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U1_dtrial_y = dimension{ U1_trial_basis_y, quad, ders, subdivision };

    auto U1_test_basis_x = bspline::create_basis(0, 1, p_test + 1, subdivision*n, rep_test);
    auto U1_dtest_x = dimension{ U1_test_basis_x, quad, ders, 1 };

    auto U1_test_basis_y = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U1_dtest_y = dimension{ U1_test_basis_y, quad, ders, 1 };


    auto U2_trial_basis_x = bspline::create_basis(0, 1, p_trial, n, rep_trial);
    auto U2_dtrial_x = dimension{ U2_trial_basis_x, quad, ders, subdivision };

    auto U2_trial_basis_y = bspline::create_basis(0, 1, p_trial + 1, n, rep_trial);
    auto U2_dtrial_y = dimension{ U2_trial_basis_y, quad, ders, subdivision };

    auto U2_test_basis_x = bspline::create_basis(0, 1, p_test, subdivision*n, rep_test);
    auto U2_dtest_x = dimension{ U2_test_basis_x, quad, ders, 1 };

    auto U2_test_basis_y = bspline::create_basis(0, 1, p_test + 1, subdivision*n, rep_test);
    auto U2_dtest_y = dimension{ U2_test_basis_y, quad, ders, 1 };


    // Sanity check
    auto trial_dim = dtrial_x.B.dofs();
    auto test_dim = dtest_x.B.dofs();

    if (trial_dim > test_dim) {
        std::cerr << "Dimension of the trial space greater than that of test space ("
                  << trial_dim << " > " << test_dim << ")" << std::endl;
        std::exit(1);
    } else {
        std::cout << "dim(U) = " << trial_dim << ", dim(V) = " << test_dim << std::endl;
    }

    auto trial = space_set{
        U1_dtrial_x, U1_dtrial_y,
        U2_dtrial_x, U2_dtrial_y,
        dtrial_x, dtrial_y
    };

    auto test = space_set{
        U1_dtest_x, U1_dtest_y,
        U2_dtest_x, U2_dtest_y,
        dtest_x, dtest_y
    };

    auto sim = stokes_dg{trial, test, steps};

    sim.run();
}
