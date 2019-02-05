#ifndef ADS_SIMULATION_BASIC_SIMULATION_3D_HPP_
#define ADS_SIMULATION_BASIC_SIMULATION_3D_HPP_

#include <cstddef>
#include <array>
#include <boost/range/counting_range.hpp>

#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"



namespace ads {

class basic_simulation_3d {
protected:
    using vector_type = lin::tensor<double, 3>;
    using vector_view = lin::tensor_view<double, 3>;
    using value_type = function_value_3d;

    using index_type = std::array<int, 3>;
    using index_1d_iter_type = boost::counting_iterator<int>;
    using index_iter_type = util::iter_product3<index_1d_iter_type, index_type>;
    using index_range = boost::iterator_range<index_iter_type>;

    using point_type = std::array<double, 3>;


    struct L2 {
        double operator ()(value_type a) const {
            return a.val * a.val;
        }
    };

    struct H10 {
        double operator ()(value_type a) const {
            return a.dx * a.dx + a.dy * a.dy + a.dz * a.dz;
        }
    };

    struct H1 {
        double operator ()(value_type a) const {
            return a.val * a.val + a.dx * a.dx + a.dy * a.dy + a.dz * a.dz;
        }
    };

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
        double dB3 = by.b[e[2]][q[2]][1][loc[2]];


        double v = B1 * B2 * B3;
        double dxv = dB1 *  B2 * B3;
        double dyv =  B1 * dB2 * B3;
        double dzv =  B1 * B2 * dB3;


        return { v, dxv, dyv, dzv };
    }

    double laplacian(index_type e, index_type q, index_type a,
                     const dimension& x, const dimension& y, const dimension& z) const {
        auto loc = dof_global_to_local(e, a, x, y, z);

        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;


        double B1  = bx.b[e[0]][q[0]][0][loc[0]];
        double B2  = by.b[e[1]][q[1]][0][loc[1]];
        double B3  = bz.b[e[2]][q[2]][0][loc[2]];

        double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];
        double ddB2 = by.b[e[1]][q[1]][2][loc[1]];
        double ddB3 = by.b[e[2]][q[2]][2][loc[2]];


        return ddB1 * B2 * B3 + B1 * ddB2 * B3 + B1 * B2 * ddB3;
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

        return e[0] >= xrange.first && e[0] <= xrange.second
            && e[1] >= yrange.first && e[1] <= yrange.second
            && e[2] >= zrange.first && e[1] <= zrange.second;
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
            global(a[0], a[1]) += local(loc[0], loc[1], loc[2]);
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

    double weigth(index_type q, const dimension& x, const dimension& y, const dimension& z) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]] * z.basis.w[q[2]];
    }

    point_type point(index_type e, index_type q, const dimension& x, const dimension& y, const dimension& z) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        double pz = z.basis.x[e[2]][q[2]];
        return { px, py, pz };
    }

private:
    auto overlapping_dofs(int dof, int begin, int end, const dimension& x) const {
        using std::min;
        using std::max;

        auto minx = max(begin, dof - x.B.degree);
        auto maxx = min(end, dof + x.B.degree + 1);

        return boost::counting_range(minx, maxx);
    }

protected:
    index_range overlapping_dofs(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        auto rx = overlapping_dofs(dof[0], 0, x.dofs(), x);
        auto ry = overlapping_dofs(dof[1], 0, y.dofs(), y);
        auto rz = overlapping_dofs(dof[2], 0, z.dofs(), z);
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

    template <typename Fun>
    void for_boundary_dofs(const dimension& x, const dimension& y, const dimension& z, Fun&& fun) const {
        for (auto jx = 0; jx < x.dofs(); ++ jx) {
            for (auto jy = 0; jy < y.dofs(); ++ jy) {
                fun({jx, jy, 0});
                fun({jx, jy, z.dofs() - 1});
            }
        }
        for (auto jx = 0; jx < x.dofs(); ++ jx) {
            for (auto jz = 1; jz < z.dofs() - 1; ++ jz) {
                fun({jx, 0, jz});
                fun({jx, y.dofs() - 1, jz});
            }
        }
        for (auto jy = 1; jy < y.dofs() - 1; ++ jy) {
            for (auto jz = 1; jz < z.dofs() - 1; ++ jz) {
                fun({0, jy, jz});
                fun({x.dofs() - 1, jy, jz});
            }
        }
    }

    bool is_boundary(int dof, const dimension& x) const {
        return dof == 0 || dof == x.dofs() - 1;
    }

    bool is_boundary(index_type dof, const dimension& x, const dimension& y, const dimension& z) const {
        return is_boundary(dof[0], x) || is_boundary(dof[1], y) || is_boundary(dof[2], z);
    }

    template <typename Norm, typename Sol>
    double norm(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm) const {
        double val = 0;

        for (auto e : elements(Ux, Ux, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weigth(q, Ux, Uy, Uz);
                value_type uu = eval(u, e, q, Ux, Uy, Uz);
                val += norm(uu) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Norm, typename Fun>
    double norm(const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm, Fun&& fun) const {
        double val = 0;

        for (auto e : elements(Ux, Ux, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weigth(q, Ux, Uy, Uz);
                auto x = point(e, q, Ux, Uy, Uz);
                auto d = fun(x);
                val += norm(d) * w * J;
            }
        }
        return std::sqrt(val);
    }

    template <typename Sol, typename Fun, typename Norm>
    double error(const Sol& u, const dimension& Ux, const dimension& Uy, const dimension& Uz, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(Ux, Ux, Uz)) {
            double J = jacobian(e, Ux, Uy, Uz);
            for (auto q : quad_points(Ux, Uy, Uz)) {
                double w = weigth(q, Ux, Uy, Uz);
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
        return error(u, Ux, Uy, Uz, norm, fun) / this->norm(Ux, Uy, Uz, norm, fun) * 100;
    }

};

}


#endif // ADS_SIMULATION_BASIC_SIMULATION_3D_HPP_
