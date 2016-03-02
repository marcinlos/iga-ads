#ifndef ADS_SIMULATION_SIMULATION_3D_HPP_
#define ADS_SIMULATION_SIMULATION_3D_HPP_

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


namespace ads {


class simulation_3d : public simulation_base {
protected:
    using vector_type = lin::tensor<double, 3>;
    using value_type = function_value_3d;

    using index_type = std::array<int, 3>;
    using index_1d_iter_type = boost::counting_iterator<int>;
    using index_iter_type = util::iter_product3<index_1d_iter_type, index_type>;
    using index_range = boost::iterator_range<index_iter_type>;

    using point_type = std::array<double, 3>;

    dimension x, y, z;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data(), z.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, z.basis, f);
    }

    double grad_dot(value_type a, value_type b) const {
        return a.dx * b.dx + a.dy * b.dy + a.dz * b.dz;
    }

    std::array<std::size_t, 3> shape() const {
        return {x.dofs(), y.dofs(), z.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
        z.factorize_matrix();
    }

    index_range elements() const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices(), z.element_indices());
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

    double weigth(index_type q) const {
        return x.basis.w[q[0]] * y.basis.w[q[1]] * z.basis.w[q[2]];
    }

    point_type point(index_type e, index_type q) const {
        double px = x.basis.x[e[0]][q[0]];
        double py = y.basis.x[e[1]][q[1]];
        double pz = z.basis.x[e[2]][q[2]];
        return { px, py, pz };
    }

    value_type eval_basis(index_type e, index_type q, index_type a) {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        int first1 = bx.first_dof(e[0]);
        int first2 = by.first_dof(e[1]);
        int first3 = bz.first_dof(e[2]);

        double B1  = bx.b[e[0]][q[0]][0][a[0] - first1];
        double B2  = by.b[e[1]][q[1]][0][a[1] - first2];
        double B3  = bz.b[e[2]][q[2]][0][a[2] - first3];
        double dB1 = bx.b[e[0]][q[0]][1][a[0] - first1];
        double dB2 = by.b[e[1]][q[1]][1][a[1] - first2];
        double dB3 = bz.b[e[2]][q[2]][1][a[2] - first3];

        double v = B1 * B2 * B3;
        double dxv = dB1 *  B2 *  B3;
        double dyv =  B1 * dB2 *  B3;
        double dzv =  B1 *  B2 * dB3;

        return { v, dxv, dyv, dzv };
    }

    value_type eval_fun(const vector_type& v, index_type e, index_type q) {
        value_type u{};
        for (auto b : dofs_on_element(e)) {
            double c = v(b[0], b[1], b[2]);
            value_type B = eval_basis(e, q, b);
            u += c * B;
        }
        return u;
    }

public:
    simulation_3d(const config_3d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , z{config.z, config.derivatives}
    , buffer{shape()}
    { }
};


}


#endif /* ADS_SIMULATION_SIMULATION_3D_HPP_ */
