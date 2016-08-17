#ifndef ADS_PROBLEMS_HEAT_HEAT_1D_HPP_
#define ADS_PROBLEMS_HEAT_HEAT_1D_HPP_

#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"


namespace ads {
namespace problems {


class heat_1d : public simulation_1d {
private:
    using Base = simulation_1d;

    vector_type u, u_prev;
    output_manager<1> output;

public:
    heat_1d(const config_1d& config)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , output{x.B, 1000}
    { }

private:
    void solve(vector_type& v) {
        v(0) = 0;
        v(x.dofs() - 1) = 0;
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        Base::prepare_matrices();
    }

    double init_state(double x) {
        double dx = x - 0.5;
        double r2 = std::min(8 * (dx * dx), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x) { return init_state(x); };
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
        const auto& bx = x.basis;
        auto& rhs = u;

        zero(rhs);
        for (element_id e = 0; e < bx.elements; ++ e) {
            double J = bx.J[e];
            int first = bx.first_dof(e);
            int last = bx.last_dof(e);

            for (int q = 0; q < bx.quad_order; ++ q) {
                double w = bx.w[q];

                for (int a = first; a <= last; ++ a) {

                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);

                    double gradient_prod = u.dx * v.dx;
                    double val = u.val * v.val - steps.dt * gradient_prod;
                    rhs(a) += val * w * J;
                }
            }
        }
    }

};



}
}



#endif /* ADS_PROBLEMS_HEAT_HEAT_1D_HPP_ */
