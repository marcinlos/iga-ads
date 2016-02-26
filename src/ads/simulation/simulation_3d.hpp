#ifndef ADS_SIMULATION_SIMULATION_3D_HPP_
#define ADS_SIMULATION_SIMULATION_3D_HPP_

#include <array>
#include <cstddef>
#include "ads/util/function_value.hpp"
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

    dimension x, y, z;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data(), z.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, z.basis, f);
    }

    std::array<std::size_t, 3> shape() const {
        return {x.dofs(), y.dofs(), z.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
        z.factorize_matrix();
    }

    value_type eval_basis(int e1, int e2, int e3, int q1, int q2, int q3, int a1, int a2, int a3) {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        int first1 = bx.first_dof(e1);
        int first2 = by.first_dof(e2);
        int first3 = bz.first_dof(e3);

        double B1  = bx.b[e1][q1][0][a1 - first1];
        double B2  = by.b[e2][q2][0][a2 - first2];
        double B3  = bz.b[e3][q3][0][a3 - first3];
        double dB1 = bx.b[e1][q1][1][a1 - first1];
        double dB2 = by.b[e2][q2][1][a2 - first2];
        double dB3 = bz.b[e3][q3][1][a3 - first3];

        double v = B1 * B2 * B3;
        double dxv = dB1 *  B2 *  B3;
        double dyv =  B1 * dB2 *  B3;
        double dzv =  B1 *  B2 * dB3;

        return { v, dxv, dyv, dzv };
    }

    value_type eval_fun(const vector_type& v, int e1, int e2, int e3, int q1, int q2, int q3) {
        int first1 = x.basis.first_dof(e1);
        int last1  = x.basis.last_dof(e1);
        int first2 = y.basis.first_dof(e2);
        int last2  = y.basis.last_dof(e2);
        int first3 = z.basis.first_dof(e3);
        int last3  = z.basis.last_dof(e3);

        value_type u{};
        for (int b1 = first1; b1 <= last1; ++ b1) {
        for (int b2 = first2; b2 <= last2; ++ b2) {
        for (int b3 = first3; b3 <= last3; ++ b3) {
            double c = v(b1, b2, b3);
            value_type B = eval_basis(e1, e2, e3, q1, q2, q3, b1, b2, b3);
            u += c * B;
        }
        }
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
