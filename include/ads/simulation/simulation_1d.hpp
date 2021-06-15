#ifndef ADS_SIMULATION_SIMULATION_1D_HPP_
#define ADS_SIMULATION_SIMULATION_1D_HPP_

#include <array>
#include <cstddef>

#include "ads/lin/tensor.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/simulation/simulation_base.hpp"
#include "ads/solver.hpp"
#include "ads/util/function_value.hpp"


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

    std::array<std::size_t, 1> local_shape() const {
        return {x.basis.dofs_per_element()};
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

    value_type eval_basis(index_type e, index_type q, index_type a) const {
        auto loc = dof_global_to_local(e, a);

        const auto& bx = x.basis;

        double v  = bx.b[e][q][0][loc];
        double dv = bx.b[e][q][1][loc];

        return { v, dv };
    }

    value_type eval_fun(const vector_type& v, index_type e, index_type q) const  {
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

    index_type dof_global_to_local(index_type e, index_type a) const {
        return  a - x.basis.first_dof(e) ;
    }

    vector_type element_rhs() const {
        return vector_type{local_shape()};
    }

    void update_global_rhs(vector_type& global, const vector_type& local, index_type e) const {
        for (auto a : dofs_on_element(e)) {
            auto loc = dof_global_to_local(e, a);
            global(a) += local(loc);
        }
    }

public:
    explicit simulation_1d(const config_1d& config);

    simulation_1d(dimension x, const timesteps_config& steps);

};

}

#endif /* ADS_SIMULATION_SIMULATION_1D_HPP_ */
