// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SIMULATION_BASIC_SIMULATION_2D_HPP
#define ADS_SIMULATION_BASIC_SIMULATION_2D_HPP

#include <array>
#include <cstddef>

#include <boost/range/counting_range.hpp>

#include "ads/lin/tensor.hpp"
#include "ads/simulation/boundary.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"

namespace ads {

class basic_simulation_2d {
public:
    virtual ~basic_simulation_2d() = default;

protected:
    using vector_type = lin::tensor<double, 2>;
    using vector_view = lin::tensor_view<double, 2>;
    using value_type = function_value_2d;

    using index_type = std::array<int, 2>;
    using index_1d_iter_type = boost::counting_iterator<int>;
    using index_iter_type = util::iter_product2<index_1d_iter_type, index_type>;
    using index_range = boost::iterator_range<index_iter_type>;

    using point_type = std::array<double, 2>;

    struct L2 {
        double operator()(value_type a) const { return a.val * a.val; }
    };

    struct H10 {
        double operator()(value_type a) const { return a.dx * a.dx + a.dy * a.dy; }
    };

    struct H1 {
        double operator()(value_type a) const { return a.val * a.val + a.dx * a.dx + a.dy * a.dy; }
    };

    value_type eval_basis(index_type e, index_type q, index_type a, const dimension& x,
                          const dimension& y) const {
        auto loc = dof_global_to_local(e, a, x, y);

        const auto& bx = x.basis;
        const auto& by = y.basis;

        double B1 = bx.b[e[0]][q[0]][0][loc[0]];
        double B2 = by.b[e[1]][q[1]][0][loc[1]];
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];

        double v = B1 * B2;
        double dxv = dB1 * B2;
        double dyv = B1 * dB2;

        return {v, dxv, dyv};
    }

    double laplacian(index_type e, index_type q, index_type a, const dimension& x,
                     const dimension& y) const {
        auto loc = dof_global_to_local(e, a, x, y);

        const auto& bx = x.basis;
        const auto& by = y.basis;

        double B1 = bx.b[e[0]][q[0]][0][loc[0]];
        double B2 = by.b[e[1]][q[1]][0][loc[1]];
        double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];
        double ddB2 = by.b[e[1]][q[1]][2][loc[1]];

        return B1 * ddB2 + ddB1 * B2;
    }

    template <typename Sol>
    value_type eval(const Sol& v, index_type e, index_type q, const dimension& x,
                    const dimension& y) const {
        value_type u{};
        for (auto b : dofs_on_element(e, x, y)) {
            double c = v(b[0], b[1]);
            value_type B = eval_basis(e, q, b, x, y);
            u += c * B;
        }
        return u;
    }

    index_range elements(const dimension& x, const dimension& y) const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices());
    }

    index_range quad_points(const dimension& x, const dimension& y) const {
        auto rx = boost::counting_range(0, x.basis.quad_order);
        auto ry = boost::counting_range(0, y.basis.quad_order);
        return util::product_range<index_type>(rx, ry);
    }

    index_range dofs_on_element(index_type e, const dimension& x, const dimension& y) const {
        auto rx = x.basis.dof_range(e[0]);
        auto ry = y.basis.dof_range(e[1]);
        return util::product_range<index_type>(rx, ry);
    }

    index_range elements_supporting_dof(index_type dof, const dimension& x,
                                        const dimension& y) const {
        auto rx = x.basis.element_range(dof[0]);
        auto ry = y.basis.element_range(dof[1]);
        return util::product_range<index_type>(rx, ry);
    }

    bool supported_in(index_type dof, index_type e, const dimension& x, const dimension& y) const {
        auto xrange = x.basis.element_ranges[dof[0]];
        auto yrange = y.basis.element_ranges[dof[1]];
        return e[0] >= xrange.first && e[0] <= xrange.second && e[1] >= yrange.first
            && e[1] <= yrange.second;
    }

    index_type dof_global_to_local(index_type e, index_type a, const dimension& x,
                                   const dimension& y) const {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        return {{a[0] - bx.first_dof(e[0]), a[1] - by.first_dof(e[1])}};
    }

    template <typename RHS>
    void update_global_rhs(RHS& global, const vector_type& local, index_type e, const dimension& x,
                           const dimension& y) const {
        for (auto a : dofs_on_element(e, x, y)) {
            auto loc = dof_global_to_local(e, a, x, y);
            global(a[0], a[1]) += local(loc[0], loc[1]);
        }
    }

    index_range dofs(const dimension& x, const dimension& y) const {
        auto rx = boost::counting_range(0, x.dofs());
        auto ry = boost::counting_range(0, y.dofs());
        return util::product_range<index_type>(rx, ry);
    }

    index_range internal_dofs(const dimension& x, const dimension& y) const {
        auto rx = boost::counting_range(1, x.dofs() - 1);
        auto ry = boost::counting_range(1, y.dofs() - 1);
        return util::product_range<index_type>(rx, ry);
    }

    double jacobian(index_type e, const dimension& x, const dimension& y) const {
        return x.basis.J[e[0]] * y.basis.J[e[1]];
    }

    double weight(index_type q, const dimension& x, const dimension& y) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]];
    }

    point_type point(index_type e, index_type q, const dimension& x, const dimension& y) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        return {px, py};
    }

    auto overlapping_dofs(int dof, int begin, int end, const dimension& x) const {
        using std::max;
        using std::min;

        auto minx = max(begin, dof - x.B.degree);
        auto maxx = min(end, dof + x.B.degree + 1);

        return boost::counting_range(minx, maxx);
    }

    index_range overlapping_dofs(index_type dof, const dimension& x, const dimension& y) const {
        auto rx = overlapping_dofs(dof[0], 0, x.dofs(), x);
        auto ry = overlapping_dofs(dof[1], 0, y.dofs(), y);
        return util::product_range<index_type>(rx, ry);
    }

    index_range overlapping_dofs(index_type dof, const dimension& Ux, const dimension& Uy,
                                 const dimension& Vx, const dimension& Vy) const {
        auto xrange = Ux.basis.element_ranges[dof[0]];
        auto yrange = Uy.basis.element_ranges[dof[1]];

        auto x0 = Vx.basis.first_dof(xrange.first);
        auto x1 = Vx.basis.last_dof(xrange.second) + 1;

        auto y0 = Vy.basis.first_dof(yrange.first);
        auto y1 = Vy.basis.last_dof(yrange.second) + 1;

        auto rx = boost::counting_range(x0, x1);
        auto ry = boost::counting_range(y0, y1);

        return util::product_range<index_type>(rx, ry);
    }

    index_range overlapping_internal_dofs(index_type dof, const dimension& x,
                                          const dimension& y) const {
        auto rx = overlapping_dofs(dof[0], 1, x.dofs() - 1, x);
        auto ry = overlapping_dofs(dof[1], 1, y.dofs() - 1, y);
        return util::product_range<index_type>(rx, ry);
    }

    int linear_index(index_type dof, const dimension& x, const dimension& y) const {
        auto order = reverse_ordering<2>({x.dofs(), y.dofs()});
        return order.linear_index(dof[0], dof[1]);
    }

    template <typename Fun>
    void for_boundary_dofs(const dimension& x, const dimension& y, Fun&& fun) const {
        for (auto jx = 0; jx < x.dofs(); ++jx) {
            fun({jx, 0});
            fun({jx, y.dofs() - 1});
        }
        for (auto jy = 1; jy < y.dofs() - 1; ++jy) {
            fun({0, jy});
            fun({x.dofs() - 1, jy});
        }
    }

    bool is_boundary(int dof, const dimension& x) const { return dof == 0 || dof == x.dofs() - 1; }

    bool is_boundary(index_type dof, const dimension& x, const dimension& y) const {
        return is_boundary(dof[0], x) || is_boundary(dof[1], y);
    }

    template <typename MT1, typename MT2>
    double kron(const MT1& A, const MT2& B, index_type i, index_type j) const {
        return A(i[0], j[0]) * B(i[1], j[1]);
    }

    template <typename RHS, typename Fun,
              typename = std::enable_if<std::is_arithmetic<std::result_of_t<Fun(double)>>{}>>
    void dirichlet_bc(RHS& u, boundary side, dimension& x, dimension& y, Fun&& fun) const {
        bool horizontal = side == boundary::top || side == boundary::bottom;
        auto& basis = horizontal ? x : y;
        const auto& other = horizontal ? y : x;

        lin::vector buf{{basis.dofs()}};
        compute_projection(buf, basis.basis, std::forward<Fun>(fun));
        lin::solve_with_factorized(basis.M, buf, basis.ctx);

        int idx = side == boundary::left || side == boundary::bottom ? 0 : other.dofs() - 1;
        for (int i = 0; i < basis.dofs(); ++i) {
            if (horizontal) {
                u(i, idx) = buf(i);
            } else {
                u(idx, i) = buf(i);
            }
        }
    }

    template <typename RHS>
    void dirichlet_bc(RHS& u, boundary side, dimension& x, dimension& y, double value) const {
        dirichlet_bc(u, side, x, y, [value](double) { return value; });
    }

    template <typename Norm, typename Fun>
    double norm(const dimension& Ux, const dimension& Uy, Norm&& norm, Fun&& fun) const {
        double val = 0;

        for (auto e : elements(Ux, Uy)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q, Ux, Uy);
                auto x = point(e, q, Ux, Uy);
                auto d = fun(x);
                val += norm(d) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Fun>
    double normL2(const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        return norm(Ux, Uy, L2{}, fun);
    }

    template <typename Fun>
    double normH1(const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        return norm(Ux, Uy, H1{}, fun);
    }

    template <typename Sol, typename Norm>
    double norm(const Sol& u, const dimension& Ux, const dimension& Uy, Norm&& norm) const {
        double val = 0;

        for (auto e : elements(Ux, Uy)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q, Ux, Uy);
                value_type uu = eval(u, e, q, Ux, Uy);
                val += norm(uu) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Sol>
    double normL2(const Sol& u, const dimension& Ux, const dimension& Uy) const {
        return norm(u, Ux, Uy, L2{});
    }

    template <typename Sol>
    double normH1(const Sol& u, const dimension& Ux, const dimension& Uy) const {
        return norm(u, Ux, Uy, H1{});
    }

    template <typename Sol, typename Fun, typename Norm>
    double error(const Sol& u, const dimension& Ux, const dimension& Uy, Norm&& norm,
                 Fun&& fun) const {
        double error = 0;

        for (auto e : elements(Ux, Uy)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q, Ux, Uy);
                auto x = point(e, q, Ux, Uy);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol, typename Fun>
    double errorL2(const Sol& u, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        return error(u, Ux, Uy, L2{}, fun);
    }

    template <typename Sol, typename Fun>
    double error_relative_L2(const Sol& u, const dimension& Ux, const dimension& Uy,
                             Fun&& fun) const {
        return error_relative(u, Ux, Uy, L2{}, fun);
    }

    template <typename Sol, typename Fun>
    double errorH1(const Sol& u, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        return error(u, Ux, Uy, H1{}, fun);
    }

    template <typename Sol, typename Fun, typename Norm>
    double error_relative(const Sol& u, const dimension& Ux, const dimension& Uy, Norm&& norm,
                          Fun&& fun) const {
        return error(u, Ux, Uy, norm, fun) / this->norm(Ux, Uy, norm, fun) * 100;
    }

    template <typename Sol, typename Fun>
    double error_relative_H1(const Sol& u, const dimension& Ux, const dimension& Uy,
                             Fun&& fun) const {
        return error_relative(u, Ux, Uy, H1{}, fun);
    }
};

}  // namespace ads

#endif  // ADS_SIMULATION_BASIC_SIMULATION_2D_HPP
