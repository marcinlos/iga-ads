#ifndef ADS_PROBLEMS_HEAT_HEAT_2D_HPP_
#define ADS_PROBLEMS_HEAT_HEAT_2D_HPP_


#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"


namespace ads {
namespace problems {


class heat_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    output_manager<2> output;
public:
    heat_2d(const config_2d& config)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , output{ x.B, y.B, 200 }
    { }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(8 * (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1) * 0;
    };

private:
    void solve(vector_type& v) {
        lin::vector buf{{ y.dofs() }};
        compute_projection(buf, y.basis, [](double y) {
            return std::sin(y * M_PI);
        });
        for (int i = 0; i < y.dofs(); ++ i) {
            v(0, i) = buf(i);
        }
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);

        output.to_file(u, "init.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs();
        solve(u);
    }

    void after_step(int iter, double /*t*/) override {
        if (iter % 100 == 0) {
            output.to_file(u, "out_%d.data", iter);
        }
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

                    double gradient_prod = grad_dot(u, v);
                    double val = u.val * v.val - steps.dt * gradient_prod;
                    rhs(a[0], a[1]) += val * w * J;
                }
            }
        }
    }
};

}
}




#endif /* ADS_PROBLEMS_HEAT_HEAT_2D_HPP_ */
