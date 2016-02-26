#ifndef ADS_SIMULATION_SIMULATION_2D_HPP_
#define ADS_SIMULATION_SIMULATION_2D_HPP_

#include <array>
#include <cstddef>
#include "ads/util/function_value.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"


namespace ads {


class simulation_2d : public simulation_base {
protected:
    using vector_type = lin::tensor<double, 2>;
    using value_type = function_value_2d;

    dimension x, y;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, f);
    }

    std::array<std::size_t, 2> shape() const {
        return {x.dofs(), y.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
    }

    value_type eval_basis(int e1, int e2, int q1, int q2, int a1, int a2) {
        const auto& bx = x.basis;
        const auto& by = y.basis;

        int first1 = bx.first_dof(e1);
        int first2 = by.first_dof(e2);

        double B1  = bx.b[e1][q1][0][a1 - first1];
        double B2  = by.b[e2][q2][0][a2 - first2];
        double dB1 = bx.b[e1][q1][1][a1 - first1];
        double dB2 = by.b[e2][q2][1][a2 - first2];

        double v = B1 * B2;
        double dxv = dB1 *  B2;
        double dyv =  B1 * dB2;

        return { v, dxv, dyv };
    }

    value_type eval_fun(const vector_type& v, int e1, int e2, int q1, int q2) {
        int first1 = x.basis.first_dof(e1);
        int last1  = x.basis.last_dof(e1);
        int first2 = y.basis.first_dof(e2);
        int last2  = y.basis.last_dof(e2);

        value_type u{};
        for (int b1 = first1; b1 <= last1; ++ b1) {
        for (int b2 = first2; b2 <= last2; ++ b2) {
            double c = v(b1, b2);
            value_type B = eval_basis(e1, e2, q1, q2, b1, b2);
            u += c * B;
        }
        }
        return u;
    }


public:
    simulation_2d(const config_2d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , buffer{shape()}
    { }
};


#endif /* ADS_SIMULATION_SIMULATION_2D_HPP_ */
