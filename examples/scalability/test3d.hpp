// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_PROBLEMS_SCALABILITY_TEST_3D_HPP_
#define ADS_PROBLEMS_SCALABILITY_TEST_3D_HPP_

#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"

namespace ads {
namespace problems {


class scalability_3d : public simulation_3d {
private:
    using Base = simulation_3d;
    vector_type u, u_prev;

    galois_executor executor;
    galois::StatTimer integration_timer{"integration"};

public:
    scalability_3d(const config_3d& config, int threads)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , executor{threads}
    { }

private:
    void solve(vector_type& v) {
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        for (int i = 0; i < x.dofs(); ++ i) {
            for (int j = 0; j < y.dofs(); ++ j) {
                for (int k = 0; k < z.dofs(); ++ k) {
                    u(i, j, k) = 1;
                }
            }
        }
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

    double forcing(double x, double y, double z) const {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double dz = z - 0.5;
        double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        return std::exp(- r) + 1 + std::cos(M_PI * x) * std::cos(M_PI * y) * std::cos(M_PI * z);
    }

    void compute_rhs() {
        integration_timer.start();
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval_fun(u_prev, e, q);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double gradient_prod = grad_dot(u, v);
                    double val = u.val * v.val - steps.dt * (gradient_prod - forcing(x[0], x[1], x[2]));
                    U(aa[0], aa[1], aa[2]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
        integration_timer.stop();
    }

    virtual void after() override {
        auto total = static_cast<double>(integration_timer.get());
        auto avg = total / steps.step_count;
        std::cout << "{ 'integration' : " << avg  << "}" << std::endl;
    }
};

}
}




#endif /* ADS_PROBLEMS_SCALABILITY_TEST_3D_HPP_ */
