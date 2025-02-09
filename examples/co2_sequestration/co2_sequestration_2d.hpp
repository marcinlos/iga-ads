// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CO2_SEQUESTRATION_2D_HPP
#define CO2_SEQUESTRATION_2D_HPP

#include <galois/Timer.h>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"

namespace ads::problems {

class co2_sequestration_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type p, s, p_prev, s_prev;

    output_manager<2> output;
    galois_executor executor{4};
    galois::StatTimer integration_timer{"integration"};
    ads::mumps::solver solver;
    // parameters
    const double mu_w; // brine viscosity
    const double mu_g; // gas viscosity
    const double K; // permeability tensor
    const double phi; // porosity
    const double rho_w; // brine density
    const double rho_g; // gas density
    const double g; // gravitational acceleration
    const double qg_x; // gas injection location - x
    const double qg_y; // gas injection location - y
    const double qg_rate; // gas injection rate
    bool verbose;

public:
    explicit co2_sequestration_2d(const config_2d& config, 
    const double mu_w, const double mu_g, const double K, const double phi, 
    const double rho_w, const double rho_g, const double g, 
    const double qg_x, const double qg_y, const double qg_rate, bool verbose)
    : Base{config}
    , p{shape()}
    , s{shape()}
    , p_prev{shape()}
    , s_prev{shape()}
    , mu_w{mu_w}
    , mu_g{mu_g}
    , K{K}
    , phi{phi}
    , rho_w{rho_w}
    , rho_g{rho_g}
    , g{g}
    , qg_x{qg_x}
    , qg_y{qg_y} 
    , qg_rate{qg_rate}
    , verbose{verbose}
    , output{x.B, y.B, 2 * config.x.elements, 2 * config.y.elements} { }

    double init_state(double x, double y) {
        /*double dx = x - 50;
        double dy = y - 8;
        double r2 = std::min(0.25 * (dx * dx + dy * dy), 1.0);
        return 0.5 * ((r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1));
        */
        return 0;
    };

    double source_g(double x, double y, double t) {
        double dx = x - qg_x;
        double dy = y - qg_y;
        double r2 = std::min(0.5 * (dx * dx + dy * dy), 1.0);
        return qg_rate * ((r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1));
    }

    double source_w(double x, double y, double t) {
        double val = 0;
        return val;
    }

private:
    void solve(vector_type& v) {
        // lin::vector buf{{y.dofs()}};
        // compute_projection(buf, y.basis, [](double y) { return std::sin(y * M_PI); });
        // for (int i = 0; i < y.dofs(); ++i) {
        //    v(0, i) = buf(i);
        // }
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        //y.fix_left();
        //y.fix_right();

        // !!all the bcs are Dirichlet for now - this is not correct for the final version!!
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(p, init);
        solve(p);

        projection(s, init);
        solve(s);

        output.to_file(p, "p.init.data");
        output.to_file(s, "s.init.data");
        if (verbose) {
            std::cout << "Initial projection computed" << std::endl;}
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(p, p_prev);
        swap(s, s_prev);
    }

    void step(int /*iter*/, double t) override {

        //solve for p
        compute_rhs_p(t);
        //compute_rhs_simple(t);
        //p(0, 0) = 0;

        dirichlet_bc(p, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(p, boundary::right, x, y, [](double t) { return 0; });
        
        ads::mumps::problem problem_p(p.data(), p.size());
        assemble_problem(problem_p);
        solver.solve(problem_p);

        // once p is solved, we can solve for s and move to the next iteration afterwards

        compute_rhs(t);
        dirichlet_bc(s, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(s, boundary::right, x, y, [](double t) { return 0; });
        solve(s);
    }

    void after_step(int iter, double /*t*/) override {
        if (iter % 10 == 0) {
            output.to_file(p, "p.out_%d.data", iter);
            output.to_file(s, "s.out_%d.data", iter);
            if (verbose) {
                std::cout << "Iteration " << iter << " passed" << std::endl;}
        }
    }

    void compute_rhs_simple(double t) {
        integration_timer.start();
        auto& rhs = p;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    //double const n = 8.0 / 100.0;
                    //double const m = 11.0 / 16.0;
                    //auto const pi = M_PI;

                    //auto const lambda = pi * pi * (n * n + m * m);
                    auto fval = source_g(x[0], x[1], t);

                    double val = fval * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
        integration_timer.stop();
    }

    void compute_rhs(double t) {
        integration_timer.start();
        auto& rhs = s;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type s = eval_fun(s_prev, e, q);
                    value_type p_here = eval_fun(p, e, q);

                    double s_val = std::clamp(s.val, 0.0, 1.0);
                    double term_1 = s_val * grad_dot(p_here, v) * K / mu_g;
                    double term_2 = s_val * v.dy * K * rho_g * g / mu_g;
                    double term_3 = v.val * source_g(x[0], x[1], t);
                    double term_extra = -1 * s.val * v.val * (x[1] >= 63);

                    double val = (term_2 + term_3 - term_1) * steps.dt / phi + s_val * v.val;
                    // brute forcing the upper ceiling on saturation FOR NOW
                    val = val + term_extra;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
        integration_timer.stop();
    }

    void compute_rhs_p(double t) {
        auto& rhs = p;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type s = eval_fun(s_prev, e, q);

                    double term_1 = K * s.dy * v.val * g * rho_w / mu_w; 
                    double term_2 = K * s.dy * v.val * g * rho_g / mu_g; 
                    double term_3 = (source_w(x[0], x[1], t) + source_g(x[0], x[1], t)) * v.val;

                    double val = term_1 - term_2 - term_3;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    void after() override {
        std::cout << "integration: " << static_cast<double>(integration_timer.get()) << std::endl;
    }

    void assemble_problem(ads::mumps::problem& problem) {
        executor.for_each(dofs(x, y), [&](auto a) {
            std::vector<std::tuple<int, int, double>> vals_buf; 
            for (auto b : overlapping_dofs(a, x, y)) {
                if (is_fixed(a, x, y))
                    continue;

                double val = 0;
                for (auto e : elements_supporting_dof(a, x, y)) {
                    if (!supported_in(b, e, x, y))
                        continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(x, y)) {
                        double w = weight(q, x, y);
                        value_type ww = eval_basis(e, q, a, x, y);
                        value_type uu = eval_basis(e, q, b, x, y);

                        double s = std::clamp(eval_fun(s_prev, e, q).val, 0.0, 1.0);
                        double diff_1 = (1 - s) / mu_w;
                        double diff_2 = s / mu_g;
                        double bwu = -1 * K * (diff_1 + diff_2) * grad_dot(uu, ww);

                        //double bwu = grad_dot(uu, ww);
                        val += bwu * w * J;
                    }
                }

                if (val != 0) {
                    int i = linear_index(a, x, y) + 1;
                    int j = linear_index(b, x, y) + 1;
                    vals_buf.push_back({i, j, val});
                    // executor.synchronized([&]() { problem.add(i, j, val); });
                }
            }

            // 1's for Dirichlet BC
            for_boundary_dofs(x, y, [&](index_type dof) {
                if (is_fixed(dof, x, y)) {
                    int i = linear_index(dof, x, y) + 1;
                    vals_buf.push_back({i, i, 1});
                    //executor.synchronized([&]() { problem.add(i, i, 1); });
                }
            });

            executor.synchronized([&]() {
                for (auto& [i, j, val] : vals_buf) {
                    problem.add(i, j, val);
                }
            });
        });
    }

    bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
        //return false; //dof[0] == 0 && dof[1] == 0; //|| dof[1] == y.dofs() - 1;
        return dof[0] == 0 || dof[0] == x.dofs() - 1; 
    }

    //bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
    //    return dof[0] == 0 || dof[0] == x.dofs() - 1 || dof[1] == 0 || dof[1] == y.dofs() - 1;
    //}

};

}  // namespace ads::problems

#endif  // HEAT_HEAT_2D_HPP
