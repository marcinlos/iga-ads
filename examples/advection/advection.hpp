// SPDX-FileCopyrightText: 2026 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADVECTION_ADVECTION_HPP
#define ADVECTION_ADVECTION_HPP

#include <fmt/base.h>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"

namespace ads {

class advection : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    output_manager<2> output;
    galois_executor executor{8};
    double error_norm_sq = 0.0;

public:
    explicit advection(const config_2d& config)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , output{x.B, y.B, 200} { }

    double init_state(double x, double y) { return fi(x, y) * sc(0); }

private:
    void solve(vector_type& v) {
        for (int i = 0; i < y.dofs(); ++i) {
            v(0, i) = 0;
            v(x.dofs() - 1, i) = 0;
        }
        for (int i = 0; i < x.dofs(); ++i) {
            v(i, 0) = 0;
            v(i, y.dofs() - 1) = 0;
        }
        Base::solve(v);
    }

    static constexpr double k = 1;

    value_type solution(double x, double y, double t) const {
        return value_type{
            sc(t) * fi(x, y),
            sc(t) * M_PI * std::cos(x * M_PI) * std::sin(y * M_PI),
            sc(t) * M_PI * std::sin(x * M_PI) * std::cos(y * M_PI),
        };
    }

    double fi(double x, double y) const { return std::sin(x * M_PI) * std::sin(y * M_PI); }

    double sc(double t) const { return std::exp(-k * t); }

    double f(double t, double x, double y) const {
        using std::sin, std::cos, std::exp;
        auto const space = (2 * M_PI * M_PI - 1) * sin(M_PI * x) * sin(M_PI * y)
            + M_PI * cos(M_PI * x) * sin(M_PI * y) + M_PI * sin(M_PI * x) * cos(M_PI * y);
        return exp(-t) * space;
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

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
    }

    void after() override {
        double T = steps.dt * steps.step_count;
        // std::cout << errorL2(T) << "  " << errorH1(T) << std::endl;
        auto error_V = std::sqrt(error_norm_sq);
        fmt::print("Error: {}\n", error_V);
    }

    void after_step(int iter, double t) override {
        if (iter % 1000 == 0) {
            output.to_file(u, "out_%d.data", iter);
            // validate(t);
        }

        auto err = errorH1(t);
        error_norm_sq += err * err * steps.dt;
    }

    void compute_rhs(double t) {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                double fx = f(t, x[0], x[1]);

                value_type u = eval_fun(u_prev, e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double gradient_prod = grad_dot(u, v);
                    double rhs = - gradient_prod - (u.dx + u.dy) * v.val + fx * v.val;
                    double val = u.val * v.val + steps.dt * rhs;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    double errorL2_rel(double t) const {
        auto sol = [&](point_type x) { return solution(x[0], x[1], t); };
        return Base::errorL2(u, x, y, sol) / normL2(x, y, sol) * 100;
    }

    double errorH1(double t) const {
        auto sol = [&](point_type x) { return solution(x[0], x[1], t); };
        // return Base::errorH1(u, x, y, sol);
        return error_par(u, x, y, sol);
    }

    template <typename Sol, typename Fun>
    double error_par(const Sol& u, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        double error = 0;

        executor.for_each(elements(), [&](index_type e) {
            double local_error = 0;
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q, Ux, Uy);
                auto x = point(e, q, Ux, Uy);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - fun(x);
                local_error += (d.dx * d.dx + d.dy * d.dy) * w * J;
            }

            executor.synchronized([&]() { error += local_error; });
        });
        return std::sqrt(error);
    }
};

}  // namespace ads

#endif  //  ADVECTION_ADVECTION_HPP
