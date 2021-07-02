// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SIMULATION_BASIC_SIMULATION_3D_HPP
#define ADS_SIMULATION_BASIC_SIMULATION_3D_HPP

#include <array>
#include <cstddef>

#include <boost/range/counting_range.hpp>

#include "ads/lin/tensor.hpp"
#include "ads/simulation/boundary.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"


namespace ads {

class basic_simulation_3d {
public:
    virtual ~basic_simulation_3d() = default;

protected:
    using vector_type = lin::tensor<double, 3>;
    using vector_view = lin::tensor_view<double, 3>;
    using value_type = function_value_3d;

    using index_type = std::array<int, 3>;
    using index_1d_iter_type = boost::counting_iterator<int>;
    using index_iter_type = util::iter_product3<index_1d_iter_type, index_type>;
    using index_range = boost::iterator_range<index_iter_type>;

    using point_type = std::array<double, 3>;

    value_type eval_basis(index_type e, index_type q, index_type a,
            const dimension& x, const dimension& y, const dimension& z) const  {
        auto loc = dof_global_to_local(e, a, x, y, z);

        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        double B1  = bx.b[e[0]][q[0]][0][loc[0]];
        double B2  = by.b[e[1]][q[1]][0][loc[1]];
        double B3  = bz.b[e[2]][q[2]][0][loc[2]];
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];
        double dB3 = bz.b[e[2]][q[2]][1][loc[2]];

        double v = B1 * B2 * B3;
        double dxv = dB1 *  B2 *  B3;
        double dyv =  B1 * dB2 *  B3;
        double dzv =  B1 *  B2 * dB3;

        return { v, dxv, dyv, dzv };
    }

    template <typename Sol>
    value_type eval(const Sol& v, index_type e, index_type q,
            const dimension& x, const dimension& y, const dimension& z) const {
        value_type u{};
        for (auto b : dofs_on_element(e, x, y, z)) {
            double c = v(b[0], b[1], b[2]);
            value_type B = eval_basis(e, q, b, x, y, z);
            u += c * B;
        }
        return u;
    }

    index_range elements(const dimension& x, const dimension& y, const dimension& z) const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices(), z.element_indices());
    }

    index_range quad_points(const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = boost::counting_range(0, x.basis.quad_order);
        auto ry = boost::counting_range(0, y.basis.quad_order);
        auto rz = boost::counting_range(0, z.basis.quad_order);
        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range dofs_on_element(index_type e, const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = x.basis.dof_range(e[0]);
        auto ry = y.basis.dof_range(e[1]);
        auto rz = z.basis.dof_range(e[2]);
        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range elements_supporting_dof(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = x.basis.element_range(dof[0]);
        auto ry = y.basis.element_range(dof[1]);
        auto rz = z.basis.element_range(dof[2]);
        return util::product_range<index_type>(rx, ry, rz);
    }

    bool supported_in(index_type dof, index_type e, const dimension& x, const dimension& y, const dimension& z) const {
        auto xrange = x.basis.element_ranges[dof[0]];
        auto yrange = y.basis.element_ranges[dof[1]];
        auto zrange = z.basis.element_ranges[dof[2]];

        return e[0] >= xrange.first && e[0] <= xrange.second &&
               e[1] >= yrange.first && e[1] <= yrange.second &&
               e[2] >= zrange.first && e[2] <= zrange.second;
    }


    index_type dof_global_to_local(index_type e, index_type a,
            const dimension& x, const dimension& y, const dimension& z) const {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;
        return {{ a[0] - bx.first_dof(e[0]), a[1] - by.first_dof(e[1]), a[2] - bz.first_dof(e[2]) }};
    }

    template <typename RHS>
    void update_global_rhs(RHS& global, const vector_type& local, index_type e,
                           const dimension& x, const dimension& y, const dimension& z) const {
        for (auto a : dofs_on_element(e, x, y, z)) {
            auto loc = dof_global_to_local(e, a, x, y, z);
            global(a[0], a[1], a[2]) += local(loc[0], loc[1], loc[2]);
        }
    }

    index_range dofs(const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = boost::counting_range(0, x.dofs());
        auto ry = boost::counting_range(0, y.dofs());
        auto rz = boost::counting_range(0, z.dofs());
        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range internal_dofs(const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = boost::counting_range(1, x.dofs() - 1);
        auto ry = boost::counting_range(1, y.dofs() - 1);
        auto rz = boost::counting_range(1, z.dofs() - 1);
        return util::product_range<index_type>(rx, ry, rz);
    }

    double jacobian(index_type e, const dimension& x, const dimension& y, const dimension& z) const {
        return x.basis.J[e[0]] * y.basis.J[e[1]] * z.basis.J[e[2]];
    }

    double weight(index_type q, const dimension& x, const dimension& y, const dimension& z) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]] * z.basis.w[q[2]];
    }

    point_type point(index_type e, index_type q, const dimension& x, const dimension& y, const dimension& z) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        double pz = z.basis.x[e[2]][q[2]];
        return { px, py, pz };
    }

    auto overlapping_dofs(int dof, int begin, int end, const dimension& x) const {
        using std::min;
        using std::max;

        auto minx = max(begin, dof - x.B.degree);
        auto maxx = min(end, dof + x.B.degree + 1);

        return boost::counting_range(minx, maxx);
    }

    index_range overlapping_dofs(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = overlapping_dofs(dof[0], 0, x.dofs(), x);
        auto ry = overlapping_dofs(dof[1], 0, y.dofs(), y);
        auto rz = overlapping_dofs(dof[1], 0, z.dofs(), z);
        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range overlapping_dofs(index_type dof,
            const dimension& Ux, const dimension& Uy, const dimension& Uz,
            const dimension& Vx, const dimension& Vy, const dimension& Vz) const {
        auto xrange = Ux.basis.element_ranges[dof[0]];
        auto yrange = Uy.basis.element_ranges[dof[1]];
        auto zrange = Uz.basis.element_ranges[dof[2]];

        auto x0 = Vx.basis.first_dof(xrange.first);
        auto x1 = Vx.basis.last_dof(xrange.second) + 1;

        auto y0 = Vy.basis.first_dof(yrange.first);
        auto y1 = Vy.basis.last_dof(yrange.second) + 1;

        auto z0 = Vz.basis.first_dof(zrange.first);
        auto z1 = Vz.basis.last_dof(zrange.second) + 1;

        auto rx = boost::counting_range(x0, x1);
        auto ry = boost::counting_range(y0, y1);
        auto rz = boost::counting_range(z0, z1);

        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range overlapping_internal_dofs(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = overlapping_dofs(dof[0], 1, x.dofs() - 1, x);
        auto ry = overlapping_dofs(dof[1], 1, y.dofs() - 1, y);
        auto rz = overlapping_dofs(dof[2], 1, z.dofs() - 1, z);
        return util::product_range<index_type>(rx, ry, rz);
    }

    int linear_index(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        auto order = reverse_ordering<3>({x.dofs(), y.dofs(), z.dofs()});
        return order.linear_index(dof[0], dof[1], dof[2]);
    }

    bool is_boundary(int dof, const dimension& x) const {
        return dof == 0 || dof == x.dofs() - 1;
    }

    bool is_boundary(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        return is_boundary(dof[0], x) || is_boundary(dof[1], y) || is_boundary(dof[2], z);
    }

    template <typename Norm, typename Fun>
    double norm(const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm, Fun&& fun) const {
        double val = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                auto x = point(e, q, Ux, Uy, Uz);
                auto d = fun(x);
                val += norm(d) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Fun>
    double normL2(const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return norm(Ux, Uy, Uz, L2, fun);
    }

    template <typename Fun>
    double normH1(const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy + a.dz * a.dz; };
        return norm(Ux, Uy, Uz, H1, fun);
    }

    template <typename Sol, typename Norm>
    double norm(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm) const {
        double val = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                value_type uu = eval(u, e, q, Ux, Uy, Uz);
                val += norm(uu) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Sol>
    double normL2(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return norm(u, Ux, Uy, Uz, L2);
    }

    template <typename Sol>
    double normH1(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy + a.dz * a.dz; };
        return norm(u, Ux, Uy, Uz, H1);
    }

    template <typename Sol, typename Fun, typename Norm>
    double error(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                auto x = point(e, q, Ux, Uy, Uz);
                value_type uu = eval(u, e, q, Ux, Uy, Uz);

                auto d = uu - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol, typename Fun, typename Norm>
    double error_relative(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz,
            Norm&& norm, Fun&& fun) const {
        double error = 0;
        double ref_norm = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                auto x = point(e, q, Ux, Uy, Uz);
                value_type uu = eval(u, e, q, Ux, Uy, Uz);
                auto fx = fun(x);

                error += norm(uu - fx) * w * J;
                ref_norm += norm(fx) * w * J;
            }
        }
        return std::sqrt(error / ref_norm);
    }

    template <typename Sol, typename Fun>
    double errorL2(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return error(u, Ux, Uy, Uz, L2, fun);
    }

    template <typename Sol, typename Fun>
    double error_relative_L2(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return error_relative(u, Ux, Uy, Uz, L2, fun);
    }

    template <typename Sol, typename Fun>
    double errorH1(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy + a.dz * a.dz; };
        return error(u, Ux, Uy, Uz, H1, fun);
    }

    template <typename Sol, typename Fun>
    double error_relative_H1(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy + a.dz * a.dz; };
        return error_relative(u, Ux, Uy, Uz, H1, fun);
    }

    template <typename Sol>
    double norm_rot(const Sol& X, const Sol& Y, const Sol& Z,
            const dimension& Ux, const dimension& Uy, const dimension& Uz) const {
        double val = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                auto x = eval(X, e, q, Ux, Uy, Uz);
                auto y = eval(Y, e, q, Ux, Uy, Uz);
                auto z = eval(Z, e, q, Ux, Uy, Uz);

                auto rx = z.dy - y.dz;
                auto ry = x.dz - z.dx;
                auto rz = y.dx - x.dy;

                val += (rx * rx + ry * ry + rz * rz) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Sol, typename FunX, typename FunY, typename FunZ>
    double error_rot(const Sol& X, const Sol& Y, const Sol& Z,
            const dimension& U1x, const dimension& U1y, const dimension& U1z,
            const dimension& U2x, const dimension& U2y, const dimension& U2z,
            const dimension& U3x, const dimension& U3y, const dimension& U3z,
            FunX&& fx, FunY&& fy, FunZ&& fz) const {
        double val = 0;

        for (auto e : elements(U1x, U1y, U1z)) {
            double J = jacobian(e, U1x, U1y, U1z);
            for (auto q : quad_points(U1x, U1y, U1z)) {
                auto p = point(e, q, U1x, U2y, U3z);
                double w = weight(q, U1x, U2y, U3z);
                auto x = eval(X, e, q, U1x, U1y, U1z) - fx(p);
                auto y = eval(Y, e, q, U2x, U2y, U2z) - fy(p);
                auto z = eval(Z, e, q, U3x, U3y, U3z) - fz(p);

                auto rx = z.dy - y.dz;
                auto ry = x.dz - z.dx;
                auto rz = y.dx - x.dy;

                val += (rx * rx + ry * ry + rz * rz) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Sol>
    double norm_div(const Sol& X, const Sol& Y, const Sol& Z,
            const dimension& Ux, const dimension& Uy, const dimension& Uz) const {
        double val = 0;

        for (auto e : elements(Ux, Uy, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weight(q, Ux, Uy, Uz);
                auto x = eval(X, e, q, Ux, Uy, Uz);
                auto y = eval(Y, e, q, Ux, Uy, Uz);
                auto z = eval(Z, e, q, Ux, Uy, Uz);

                auto v = x.dx + y.dy + z.dz;
                val += v * v * w * J;
            }
        }
        return std::sqrt(val);
    }

};

}

#endif // ADS_SIMULATION_BASIC_SIMULATION_3D_HPP
