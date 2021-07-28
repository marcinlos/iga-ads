// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_SIMULATION_SIMULATION_3D_HPP
#define ADS_SIMULATION_SIMULATION_3D_HPP

#include <array>
#include <cstddef>

#include <boost/range/counting_range.hpp>

#include "ads/lin/tensor.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/simulation/simulation_base.hpp"
#include "ads/solver.hpp"
#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"
#include "basic_simulation_3d.hpp"

namespace ads {

class simulation_3d : public basic_simulation_3d, public simulation_base {
public:
    using basic_simulation_3d::dof_global_to_local;
    using basic_simulation_3d::dofs;
    using basic_simulation_3d::dofs_on_element;
    using basic_simulation_3d::elements;
    using basic_simulation_3d::elements_supporting_dof;
    using basic_simulation_3d::eval;
    using basic_simulation_3d::eval_basis;
    using basic_simulation_3d::jacobian;
    using basic_simulation_3d::point;
    using basic_simulation_3d::quad_points;
    using basic_simulation_3d::update_global_rhs;
    using basic_simulation_3d::weight;

    dimension x, y, z;
    vector_type buffer;

    void solve(vector_type& rhs) { ads_solve(rhs, buffer, x.data(), y.data(), z.data()); }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, z.basis, f);
    }

    double grad_dot(value_type a, value_type b) const {
        return a.dx * b.dx + a.dy * b.dy + a.dz * b.dz;
    }

    std::array<int, 3> shape() const { return {x.dofs(), y.dofs(), z.dofs()}; }

    std::array<int, 3> local_shape() const {
        return {x.basis.dofs_per_element(), y.basis.dofs_per_element(), z.basis.dofs_per_element()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
        z.factorize_matrix();
    }

    index_range elements() const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices(),
                                               z.element_indices());
    }

    index_range quad_points() const {
        auto rx = boost::counting_range(0, x.basis.quad_order);
        auto ry = boost::counting_range(0, y.basis.quad_order);
        auto rz = boost::counting_range(0, z.basis.quad_order);
        return util::product_range<index_type>(rx, ry, rz);
    }

    index_range dofs_on_element(index_type e) const {
        auto rx = x.basis.dof_range(e[0]);
        auto ry = y.basis.dof_range(e[1]);
        auto rz = z.basis.dof_range(e[2]);
        return util::product_range<index_type>(rx, ry, rz);
    }

    double jacobian(index_type e) const {
        return x.basis.J[e[0]] * y.basis.J[e[1]] * z.basis.J[e[2]];
    }

    double weight(index_type q) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]] * z.basis.w[q[2]];
    }

    point_type point(index_type e, index_type q) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        double pz = z.basis.x[e[2]][q[2]];
        return {px, py, pz};
    }

    value_type eval_basis(index_type e, index_type q, index_type a) const {
        auto loc = dof_global_to_local(e, a);

        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        double B1 = bx.b[e[0]][q[0]][0][loc[0]];
        double B2 = by.b[e[1]][q[1]][0][loc[1]];
        double B3 = bz.b[e[2]][q[2]][0][loc[2]];
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];
        double dB3 = bz.b[e[2]][q[2]][1][loc[2]];

        double v = B1 * B2 * B3;
        double dxv = dB1 * B2 * B3;
        double dyv = B1 * dB2 * B3;
        double dzv = B1 * B2 * dB3;

        return {v, dxv, dyv, dzv};
    }

    value_type eval_fun(const vector_type& v, index_type e, index_type q) const {
        value_type u{};
        for (auto b : dofs_on_element(e)) {
            double c = v(b[0], b[1], b[2]);
            value_type B = eval_basis(e, q, b);
            u += c * B;
        }
        return u;
    }

    index_type dof_global_to_local(index_type e, index_type a) const {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        return {{a[0] - bx.first_dof(e[0]), a[1] - by.first_dof(e[1]), a[2] - bz.first_dof(e[2])}};
    }

    vector_type element_rhs() const { return vector_type{local_shape()}; }

    void update_global_rhs(vector_type& global, const vector_type& local, index_type e) const {
        for (auto a : dofs_on_element(e)) {
            auto loc = dof_global_to_local(e, a);
            global(a[0], a[1], a[2]) += local(loc[0], loc[1], loc[2]);
        }
    }

    explicit simulation_3d(const config_3d& config);

    simulation_3d(dimension x, dimension y, dimension z, const timesteps_config& steps);
};

}  // namespace ads

#endif  // ADS_SIMULATION_SIMULATION_3D_HPP
