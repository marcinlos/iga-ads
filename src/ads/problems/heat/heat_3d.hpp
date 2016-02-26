#ifndef ADS_PROBLEMS_HEAT_HEAT_3D_HPP_
#define ADS_PROBLEMS_HEAT_HEAT_3D_HPP_

#include "ads/output_manager.hpp"
#include "ads/problems/heat.hpp"

#include "ads/simulation.hpp"
#include <complex>

namespace ads {
namespace problems {


class heat_3d : public simulation_3d {
private:
    using Base = simulation_1d;
    vector_type u, u_prev;

public:
    heat_3d(const config_3d& config)
    : simulation_3d{config}
    , u{shape()}
    , u_prev{shape()}
    { }

    double init_state(double x, double y, double z) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double dz = z - 0.5;
        double r2 = std::min(8 * (dx * dx + dy * dy + dz * dz), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    };

private:
    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y, double z) { return init_state(x, y, z); };
        projection(u, init);
        solve(u);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs();
        solve(u);
    }

    void compute_rhs() {
        auto& rhs = u;

        const auto& bx = x.basis;
        const auto& by = y.basis;
        const auto& bz = z.basis;

        zero(rhs);
        for (element_id e1 = 0; e1 < bx.elements; ++ e1) {
        for (element_id e2 = 0; e2 < by.elements; ++ e2) {
        for (element_id e3 = 0; e3 < bz.elements; ++ e3) {
            double J = bx.J[e1] * by.J[e2] * bz.J[e3];
            int first1 = bx.first_dof(e1);
            int last1 = bx.last_dof(e1);
            int first2 = by.first_dof(e2);
            int last2 = by.last_dof(e2);
            int first3 = bz.first_dof(e3);
            int last3 = bz.last_dof(e3);

            for (int q1 = 0; q1 < bx.quad_order; ++ q1) {
            for (int q2 = 0; q2 < by.quad_order; ++ q2) {
            for (int q3 = 0; q3 < bz.quad_order; ++ q3) {
                double w = bx.w[q1] * by.w[q2] * bz.w[q3];

                for (int a1 = first1; a1 <= last1; ++ a1) {
                for (int a2 = first2; a2 <= last2; ++ a2) {
                for (int a3 = first3; a3 <= last3; ++ a3) {
                    value_type v = eval_basis(e1, e2, e3, q1, q2, q3, a1, a2, a3);
                    value_type u = eval_fun(u_prev, e1, e2, e3, q1, q2, q3);

                    double gradient_prod = u.dx * v.dx + u.dy * v.dy + u.dz * v.dz;
                    double val = u.val * v.val - steps.dt * gradient_prod;
                    rhs(a1, a2, a3) += val * w * J;
                }
                }
                }
            }
            }
            }
        }
        }
        }
    }

};

}
}




#endif /* ADS_PROBLEMS_HEAT_HEAT_3D_HPP_ */
