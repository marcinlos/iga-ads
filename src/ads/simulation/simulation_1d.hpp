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

    dimension x;

    void solve(vector_type& rhs) {
        ads_solve(rhs, x.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, f);
    }

    std::array<std::size_t, 1> shape() const {
        return {x.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
    }

    value_type eval_basis(int e, int q, int a) {
        const auto& bx = x.basis;
        int first = bx.first_dof(e);

        double v  = bx.b[e][q][0][a - first];
        double dv = bx.b[e][q][1][a - first];

        return { v, dv };
    }

    value_type eval_fun(const vector_type& v, int e, int q) {
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
