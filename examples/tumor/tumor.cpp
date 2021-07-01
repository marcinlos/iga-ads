// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "tumor.hpp"


using ads::config_2d;

namespace tumor {

    tumor_2d::tumor_2d(const config_2d& config, const params& params, int save_every, vasc::vasculature vasculature)
    : Base{config}
    , now{ shape() }
    , prev{ shape() }
    , p{ params }
    , vasculature{ std::move(vasculature) }
    , output{ x.B, y.B, 200 }
    , save_every{ save_every }
    , xctx{ x.B.degree }
    , yctx{ x.B.degree }
    , xdctx{ x.B.degree, 1 }
    , ydctx{ x.B.degree, 1 }
    { }

    void tumor_2d::save_to_file(int iter) {
        output.to_file(now.b, "tumor_%d.data", iter);
        output.to_file(now.c, "taf_%d.data", iter);
        output.to_file(now.M, "ecm_%d.data", iter);
        output.to_file(now.A, "degraded_ecm_%d.data", iter);
        plot_vasculature(iter);
    }

    void tumor_2d::before() {
        prepare_matrices();

        auto tumor = [this](double x, double y) { return init_tumor(x, y); };
        auto m = [this](double x, double y) { return init_M(x, y); };
        projection(now.b, tumor);
        projection(now.A, constant(0));
        projection(now.M, m);

        solve_all();

        save_to_file(0);

        vasculature.discretize();
        plot_vasculature(0);
    }

    void tumor_2d::solve_all() {
        apply_boundary_conditions(now.b);
        solve(now.b);

        apply_boundary_conditions(now.c);
        solve(now.c);

        apply_boundary_conditions(now.M);
        solve(now.M);

        apply_boundary_conditions(now.A);
        solve(now.A);
    }

    void tumor_2d::compute_rhs() {
        executor.for_each(elements(), [&](index_type e) {
            state<Dim> loc{ local_shape() };

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    value_type b = ensure_positive(eval_fun(prev.b, e, q));
                    value_type c = ensure_positive(eval_fun(prev.c, e, q));
                    value_type M = ensure_positive(eval_fun(prev.M, e, q));
                    value_type A = ensure_positive(eval_fun(prev.A, e, q));

                    auto nx = normalize(x);
                    double o = vasculature.oxygen_level(nx[0], nx[1]) * p.o_max;

                    // tumor density
                    double b_src = 0;
                    double b_sink = 0;

                    // tumor source/sink
                    if (o >= p.o_prol_TC) {
                        double tbA = p.tau_b * A.val;
                        double s = 1 + tbA / (1 + tbA) * p.P_b;
                        b_src = b.val / p.t_prol_TC * s * (1 - b.val / p.c_b_max);
                    }
                    if (o < p.o_death_TC) {
                        b_sink = -b.val / p.t_death_TC;
                    }

                    double D_b = p.skin.diffusion(x[0], x[1], x[1]);
                    double grad_Av = grad_dot(A, v);
                    double grad_Pv = b.val >= p.c_b_norm ? grad_dot(b, v) / (p.c_b_max - p.c_b_norm) : 0;

                    double divJv = D_b * b.val * (grad_Pv + p.r_b * grad_Av);

                    double bv = - divJv + (b_src + b_sink) * v.val;
                    val(loc.b, aa) += (b.val * v.val + bv * steps.dt) * w * J;

                    // ECM evolution
                    double Mv = - p.beta_m * M.val * b.val * v.val;
                    val(loc.M, aa) += (M.val * v.val + Mv * steps.dt) * w * J;

                    double Av = (p.gamma_a * M.val * b.val - p.gamma_oA * A.val) * v.val - p.chi_aA * grad_Av;
                    val(loc.A, aa) += (A.val * v.val + Av * steps.dt) * w * J;

                    // TAF
                    double c_src = o < p.o_death_TC ? b.val * (1 - c.val): 0;
                    double cv = - p.diff_c * grad_dot(c, v) + (c_src - p.cons_c * o * c.val) * v.val;
                    val(loc.c, aa) += (c.val * v.val + cv * steps.dt) * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(now.b, loc.b, e);
                update_global_rhs(now.M, loc.M, e);
                update_global_rhs(now.A, loc.A, e);
                update_global_rhs(now.c, loc.c, e);
            });
        });
    }

    void tumor_2d::update_vasculature(int iter) {
        int next_iter = iter + 1;
        if (next_iter % vasc_update_every == 0) {
            auto taf = [&,this](double px, double py) {
                double xx = x.a + px * (x.b - x.a);
                double yy = y.a + py * (y.b - y.a);
                return ads::bspline::eval_ders(xx, yy, now.c, this->x.B, this->y.B, xdctx, ydctx);
            };

            auto tumor = [&,this](double px, double py) {
                double xx = x.a + px * (x.b - x.a);
                double yy = y.a + py * (y.b - y.a);
                return ads::bspline::eval(xx, yy, now.b, this->x.B, this->y.B, xctx, yctx);
            };

            vasculature.update(tumor, taf, steps.dt);
            vasculature.discretize();
        }
    }

}
