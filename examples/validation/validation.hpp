// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef VALIDATION_VALIDATION_HPP
#define VALIDATION_VALIDATION_HPP

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"


namespace ads::problems {

class validation : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    output_manager<2> output;
    galois_executor executor{8};

public:
    validation(const config_2d& config)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , output{ x.B, y.B, 200 }
    { }

    double init_state(double x, double y) {
        return fi(x, y) * sc(0);
    };

private:
    void solve(vector_type& v) {
        for (int i = 0; i < y.dofs(); ++ i) {
            v(0, i) = 0;
            v(x.dofs() - 1, i) = 0;
        }
        for (int i = 0; i < x.dofs(); ++ i) {
            v(i, 0) = 0;
            v(i, y.dofs() - 1) = 0;
        }
        Base::solve(v);
    }

    static constexpr double k = 2 * M_PI * M_PI;

    value_type solution(double x, double y, double t) const {
        return value_type{
            sc(t) * fi(x, y),
            sc(t) * M_PI * std::cos(x * M_PI) * std::sin(y * M_PI),
            sc(t) * M_PI * std::sin(x * M_PI) * std::cos(y * M_PI)
        };
    }

    double fi(double x, double y) const {
        return std::sin(x * M_PI) * std::sin(y * M_PI);
    }

    double sc(double t) const {
        return std::exp(-k*t);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
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

    void after() override {
        double T = steps.dt * steps.step_count;
        std::cout << errorL2(T) << "  " << errorH1(T) << std::endl;
    }

    void after_step(int /*iter*/, double /*t*/) override {
        // if (iter % 1000 == 0) {
        //     output.to_file(u, "out_%d.data", iter);
        //     validate(t);
        // }
    }

    void compute_rhs() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type u = eval_fun(u_prev, e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double gradient_prod = grad_dot(u, v);
                    double val = u.val * v.val - steps.dt * gradient_prod ;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    double errorL2(double t) const {
        auto sol = [&](point_type x) { return solution(x[0], x[1], t); };
        return Base::errorL2(u, x, y, sol) / normL2(x, y, sol) * 100;
    }

    double errorH1(double t) const {
        auto sol = [&](point_type x) { return solution(x[0], x[1], t); };
        return Base::errorH1(u, x, y, sol) / normH1(x, y, sol) * 100;
    }
};

}

#endif // VALIDATION_VALIDATION_HPP
