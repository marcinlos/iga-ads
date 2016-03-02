#ifndef ADS_SIMULATION_SIMULATION_1D_HPP_
#define ADS_SIMULATION_SIMULATION_1D_HPP_

#include <array>
#include <cstddef>
#include "ads/util/function_value.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"


namespace ads {


class simulation_1d : public simulation_base {
protected:
    using vector_type = lin::vector;
    using value_type = function_value_1d;

    using index_type = int;
    using index_iter_type = boost::counting_iterator<int>;
    using index_range = boost::iterator_range<index_iter_type>;

    using point_type = double;

    dimension x;

    void solve(vector_type& rhs) {
        ads_solve(rhs, x.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, f);
    }

    double grad_dot(value_type a, value_type b) const {
        return a.dx * b.dx;
    }

    std::array<std::size_t, 1> shape() const {
        return {x.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
    }

    index_range elements() const {
        return x.element_indices();
    }

    index_range quad_points() const {
        return boost::counting_range(0, x.basis.quad_order);
    }

    index_range dofs_on_element(index_type e) const {
        return x.basis.dof_range(e);
    }

    double jacobian(index_type e) const {
        return x.basis.J[e];
    }

    double weigth(index_type q) const {
        return x.basis.w[q];
    }

    point_type point(index_type e, index_type q) const {
        return x.basis.x[e][q];
    }

    value_type eval_basis(index_type e, index_type q, index_type a) {
        const auto& bx = x.basis;
        int first = bx.first_dof(e);

        double v  = bx.b[e][q][0][a - first];
        double dv = bx.b[e][q][1][a - first];

        return { v, dv };
    }

    value_type eval_fun(const vector_type& v, index_type e, index_type q) {
        int first = x.basis.first_dof(e);
        int last  = x.basis.last_dof(e);

        value_type u{};
        for (int b1 = first; b1 <= last; ++ b1) {
            double c = v(b1);
            value_type B = eval_basis(e, q, b1);
            u += c * B;
        }
        return u;
    }

public:
    simulation_1d(const config_1d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    { }
};


}


#endif /* ADS_SIMULATION_SIMULATION_1D_HPP_ */
