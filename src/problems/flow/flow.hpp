#ifndef ADS_PROBLEMS_FLOW_FLOW_HPP_
#define ADS_PROBLEMS_FLOW_FLOW_HPP_

#include <cmath>
#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"
#include "problems/flow/pumps.hpp"
#include "problems/flow/environment.hpp"


namespace ads {
namespace problems {


class flow : public simulation_3d {
private:
    using Base = simulation_3d;
    vector_type u, u_prev;

    galois_executor executor{4};

    environment env {1};
    lin::tensor<double, 6> kq;

public:
    flow(const config_3d& config)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , kq{{
        x.basis.elements, y.basis.elements, z.basis.elements,
        x.basis.quad_order + 1, y.basis.quad_order + 1, z.basis.quad_order + 1
    }}
    { }

    double init_state(double x, double y, double z) {
        double r = 0.1;
        double R = 0.25;
        return ads::bump(r, R, x, y, z);
    };

private:
    void before() override {
        fill_permeability_map();
        prepare_matrices();

        auto init = [this](double x, double y, double z) { return init_state(x, y, z); };
        projection(u, init);
        solve(u);
    }

    void fill_permeability_map() {
        for (auto e : elements()) {
            for (auto q : quad_points()) {
                auto x = point(e, q);
                kq(e[0], e[1], e[2], q[0], q[1], q[2]) = env.permeability(x[0], x[1], x[2]);
            }
        }
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
    }

    void compute_rhs(double t) {
        auto& rhs = u;

        zero(rhs);
        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);

                double mi = 10;
                double k = permeability(e, q);
                value_type u = eval_fun(u_prev, e, q);
                double h = forcing(x, t);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double val = - k * std::exp(mi * u.val) * grad_dot(u, v) + h * v.val;
                    U(aa[0], aa[1], aa[2]) += (u.val * v.val + steps.dt * val) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, U, e); });
        });
    }

    double energy(const vector_type& u) const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                value_type a = eval_fun(u, e, q);
                E += a.val * a.val * w * J;
            }
        }
        return E;
    }

    void after_step(int iter, double /*t*/) override {
        if (iter % 10 == 0) {
            std::cout << "Energy: " << energy(u) << std::endl;
        }
    }

    double permeability(index_type e, index_type q) const {
        return kq(e[0], e[1], e[2], q[0], q[1], q[2]);
    }

    double forcing(point_type x, double /*t*/) const {
        using std::sin;
        double pi2 = 2 * M_PI;
        return sin(pi2 * x[0]) * sin(pi2 * x[1]) * sin(pi2 * x[2]);
    }
};

}
}

#endif /* ADS_PROBLEMS_FLOW_FLOW_HPP_ */
