#ifndef ADS_SIMULATION_SIMULATION_2D_HPP_
#define ADS_SIMULATION_SIMULATION_2D_HPP_

#include <array>
#include <cstddef>
#include <array>
#include <boost/range/counting_range.hpp>

#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"
#include "ads/simulation/basic_simulation_2d.hpp"


namespace ads {


class simulation_2d : public basic_simulation_2d, public simulation_base {
protected:
    using basic_simulation_2d::eval;
    using basic_simulation_2d::eval_basis;
    using basic_simulation_2d::elements;
    using basic_simulation_2d::quad_points;
    using basic_simulation_2d::dofs_on_element;
    using basic_simulation_2d::elements_supporting_dof;
    using basic_simulation_2d::dof_global_to_local;
    using basic_simulation_2d::update_global_rhs;
    using basic_simulation_2d::dofs;
    using basic_simulation_2d::jacobian;
    using basic_simulation_2d::weigth;
    using basic_simulation_2d::point;


    dimension x, y;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, f);
    }

    double grad_dot(value_type a, value_type b) const {
        return a.dx * b.dx + a.dy * b.dy;
    }

    std::array<std::size_t, 2> shape() const {
        return {x.dofs(), y.dofs()};
    }

    std::array<std::size_t, 2> local_shape() const {
        return {x.basis.dofs_per_element(), y.basis.dofs_per_element()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
    }

    index_range elements() const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices());
    }

    index_range dofs() const {
        auto rx = boost::counting_range(0, x.basis.dofs);
        auto ry = boost::counting_range(0, y.basis.dofs);
        return util::product_range<index_type>(rx, ry);
    }

    index_range quad_points() const {
        auto rx = boost::counting_range(0, x.basis.quad_order);
        auto ry = boost::counting_range(0, y.basis.quad_order);
        return util::product_range<index_type>(rx, ry);
    }

    index_range dofs_on_element(index_type e) const {
        auto rx = x.basis.dof_range(e[0]);
        auto ry = y.basis.dof_range(e[1]);
        return util::product_range<index_type>(rx, ry);
    }

    index_range elements_supporting_dof(index_type dof) const {
        auto rx = x.basis.element_range(dof[0]);
        auto ry = y.basis.element_range(dof[1]);
        return util::product_range<index_type>(rx, ry);
    }

    double jacobian(index_type e) const {
        return x.basis.J[e[0]] * y.basis.J[e[1]];
    }

    double weigth(index_type q) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]];
    }

    point_type point(index_type e, index_type q) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        return { px, py };
    }

    value_type eval_basis(index_type e, index_type q, index_type a) const  {
        auto loc = dof_global_to_local(e, a);

        const auto& bx = x.basis;
        const auto& by = y.basis;

        double B1  = bx.b[e[0]][q[0]][0][loc[0]];
        double B2  = by.b[e[1]][q[1]][0][loc[1]];
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];

        double v = B1 * B2;
        double dxv = dB1 *  B2;
        double dyv =  B1 * dB2;

        return { v, dxv, dyv };
    }

    value_type eval_fun(const vector_type& v, index_type e, index_type q) const {
        value_type u{};
        for (auto b : dofs_on_element(e)) {
            double c = v(b[0], b[1]);
            value_type B = eval_basis(e, q, b);
            u += c * B;
        }
        return u;
    }

    index_type dof_global_to_local(index_type e, index_type a) const {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        return {{ a[0] - bx.first_dof(e[0]), a[1] - by.first_dof(e[1]) }};
    }

    vector_type element_rhs() const {
        return vector_type{local_shape()};
    }

    void update_global_rhs(vector_type& global, const vector_type& local, index_type e) const {
        for (auto a : dofs_on_element(e)) {
            auto loc = dof_global_to_local(e, a);
            global(a[0], a[1]) += local(loc[0], loc[1]);
        }
    }


public:
    explicit simulation_2d(const config_2d& config);

    simulation_2d(dimension x, dimension y, const timesteps_config& steps);
};

}


#endif /* ADS_SIMULATION_SIMULATION_2D_HPP_ */
