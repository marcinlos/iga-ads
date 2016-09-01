#ifndef ADS_PROBLEMS_HEAT_HEAT_3D_HPP_
#define ADS_PROBLEMS_HEAT_HEAT_3D_HPP_

#include "ads/simulation.hpp"

namespace ads {
namespace problems {


class heat_3d : public simulation_3d {
private:
    using Base = simulation_3d;
    vector_type u, u_prev;

public:
    heat_3d(const config_3d& config)
    : Base{config}
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

        zero(rhs);
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);

                    double gradient_prod = u.dx * v.dx + u.dy * v.dy + u.dz * v.dz;
                    double val = u.val * v.val - steps.dt * gradient_prod;
                    rhs(a[0], a[1], a[2]) += val * w * J;
                }
            }
        }
    }
};

}
}


#endif /* ADS_PROBLEMS_HEAT_HEAT_3D_HPP_ */
