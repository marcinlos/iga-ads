#ifndef ADS_PROBLEMS_TUMOR_TUMOR_HPP_
#define ADS_PROBLEMS_TUMOR_TUMOR_HPP_


#include "problems/tumor/skin.hpp"
#include "problems/tumor/state.hpp"
#include "problems/tumor/vasculature.hpp"

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"

#include "ads/executor/galois.hpp"


#include <cmath>
#include <boost/format.hpp>


namespace ads {
namespace tumor {



struct params {
    double c_b_norm = 1; // normal concentration of tumor cells
    double c_b_max = 2; // maximal concentration of tumor cells
    double tau_b = 0.5; // ?

    double o_prol_TC = 0.1;   // oxygen proliferation threshold
    double o_death_TC = 0.01; // oxygen survival threshold
    double t_prol_TC = 10;
    double t_death_TC = 100;
    double P_b = 0.001; // stimulated mitosis rate
    double r_b = 1e-4;  // chemoattractant sensitivity, 1-3-5 x 10^-4

    double beta_m = 1;      // ???? some ECM decay coefficient
    double gamma_a = 0.5;   // production rate of attractants
    double chi_aA = 0.01;   // diffusion of degraded ECM
    double gamma_oA = 0.01; // decay of degraded ECM

    double rho_n = 0.28;  // haptotactic cell migration
    double D_n = 0.0003;  // diffusion of endothelial cells
    double chi_n = 0.38;  // chemotactic cell migration
    double delta_n = 0.6; // chemotactic constant

    double beta_f = 0.05; // production rate of fibronectin
    double gamma_f = 0.1; // degradation rate of fibronectin

    double alpha_m = 0.000001; // MDA generation rate
    double epsilon_m = 0.01;   // MDA diffusion coefficient
    double upsilon_m = 3;      // MDA degradation rage

    double diff_c = 0.01; // TAF diffusion rate
    double cons_c = 0.3;  // TAF consumption

    skin_model skin;

    double init_M = 0.015;
};

class tumor : public simulation_2d {
private:
    using Base = simulation_2d;

    state now, prev;
    params p;

    vasc::vasculature vasculature;

    output_manager<2> output;

    int save_every = 100;

    int vasc_update_every = 10;

    galois_executor executor{8};



    bspline::eval_ctx xctx;
    bspline::eval_ctx yctx;
    bspline::eval_ders_ctx xdctx;
    bspline::eval_ders_ctx ydctx;

public:
    tumor(const config_2d& config, const params& params, vasc::vasculature vasculature)
    : Base{config}
    , now{ shape() }
    , prev{ shape() }
    , p{ params }
    , vasculature{ std::move(vasculature) }
    , output{ x.B, y.B, 100 }
    , xctx{ x.B.degree }
    , yctx{ x.B.degree }
    , xdctx{ x.B.degree, 1 }
    , ydctx{ x.B.degree, 1 }
    { }


private:

    struct constant {
        double value;

        constant(double v)
        : value{ v }
        { }

        double operator ()(double, double) const {
            return value;
        }
    };

    double init_tumor(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(12 * (dx * dx + dy * dy), 1.0);
        return 0.8 * (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

    double init_fibronectin(double x, double y) {
        auto lay = p.skin.layer_at(x, y, x);
        return lay == skin_model::layer::dermis || lay == skin_model::layer::hypodermis ? 0.8 : 0;
    }

    double init_M(double x, double y) {
        auto lay = p.skin.layer_at(x, y, x);
        return lay == skin_model::layer::dermis ? 1 : p.init_M;
    }

    void dirichlet() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();
    }

    void apply_boundary_conditions(vector_type& v) {
        for (int i = 0; i < y.dofs(); ++ i) {
            v(0, i) = 0;
            v(x.dofs() - 1, i) = 0;
        }
        for (int i = 0; i < x.dofs(); ++ i) {
            v(i, 0) = 0;
            v(i, y.dofs() - 1) = 0;
        }
    }

    void prepare_matrices() {
        dirichlet();
        Base::prepare_matrices();
    }

    void plot_vasculature(int iter) {
        using boost::format;
        vasculature.plot_veins(str(format("vasculature_%d.data") % iter));
        vasculature.plot_oxygen(str(format("oxygen_%d.data") % iter));
    }

    void before() override {
        prepare_matrices();

        auto tumor = [this](double x, double y) { return init_tumor(x, y); };
        auto fibro = [this](double x, double y) { return init_fibronectin(x, y); };
        auto m = [this](double x, double y) { return init_M(x, y); };
        projection(now.b, tumor);
        projection(now.A, constant(0));
        projection(now.M, m);
        projection(now.f, fibro);

        solve_all();

        save_to_file(0);

        vasculature.discretize();
        plot_vasculature(0);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);

        now.clear();
    }

    void step(int iter, double /*t*/) override {
        compute_rhs();
        solve_all();
        update_vasculature(iter);
    }

    void save_to_file(int iter) {
        output.to_file(now.b, "tumor_%d.data", iter);
        output.to_file(now.n, "endothelial_%d.data", iter);
        output.to_file(now.f, "fibronectin_%d.data", iter);
        output.to_file(now.m, "mde_%d.data", iter);
        output.to_file(now.c, "taf_%d.data", iter);
        output.to_file(now.M, "ecm_%d.data", iter);
        output.to_file(now.A, "degraded_ecm_%d.data", iter);
        plot_vasculature(iter);
    }

    void after_step(int iter, double /*t*/) override {
        int next_iter = iter + 1;
        if (next_iter % save_every == 0) {
            save_to_file(next_iter);
        }
    }


    void solve_all() {
        apply_boundary_conditions(now.b);
        solve(now.b);

        apply_boundary_conditions(now.c);
        solve(now.c);

        apply_boundary_conditions(now.n);
        solve(now.n);

        apply_boundary_conditions(now.f);
        solve(now.f);

        apply_boundary_conditions(now.m);
        solve(now.m);

        apply_boundary_conditions(now.M);
        solve(now.M);

        apply_boundary_conditions(now.A);
        solve(now.A);
    }

    value_type ensure_positive(value_type v) const {
        if (v.val < 0) {
            v.val = 0;
        }
        return v;
    }

    void compute_rhs() {
        executor.for_each(elements(), [&](index_type e) {
            state loc{ local_shape() };

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    value_type b = ensure_positive(eval_fun(prev.b, e, q));
                    value_type c = ensure_positive(eval_fun(prev.c, e, q));
                    value_type n = ensure_positive(eval_fun(prev.n, e, q));
                    value_type f = ensure_positive(eval_fun(prev.f, e, q));
                    value_type m = ensure_positive(eval_fun(prev.m, e, q));

                    value_type M = ensure_positive(eval_fun(prev.M, e, q));
                    value_type A = ensure_positive(eval_fun(prev.A, e, q));

                    double o = vasculature.oxygen_level(x[0], x[1]);

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

//                    double D_b = p.skin.diffusion(x[0], x[1], x[2]);
                    double D_b = p.skin.diffusion(x[0], x[1], x[0]) * 0.01; // ????
                    double grad_Av = grad_dot(A, v);
                    double grad_Pv = b.val >= p.c_b_norm ? grad_dot(b, v) / (p.c_b_max - p.c_b_norm) : 0;

                    double divJv = D_b * b.val * (grad_Pv + p.r_b * grad_Av);

                    double bv = - divJv + (b_src + b_sink) * v.val;
                    val(loc.b, aa) += (b.val * v.val + bv * steps.dt) * w * J;

                    // endothelial cells
                    double X = p.chi_n / (1 + p.delta_n * c.val);
//                    double nv = - p.D_n * grad_dot(n, v) + grad_dot(n.val * X * c, v) + p.rho_n * grad_dot(n.val * f, v);
                    double nv = grad_dot(-p.D_n * n + n.val * X * c + n.val * f, v);
                    val(loc.n, aa) += (n.val * v.val + nv * steps.dt) * w * J;

                    // fibronectin
                    double fv = p.beta_f * n.val - p.gamma_f * m.val * f.val;
                    val(loc.f, aa) += (f.val * v.val + fv * steps.dt) * w * J;

                    // MDE
                    double mv = p.alpha_m * n.val - p.epsilon_m * grad_dot(m, v) - p.upsilon_m * m.val;
                    val(loc.m, aa) += (m.val * v.val + mv * steps.dt) * w * J;

                    // ECM evolution
                    double Mv = - p.beta_m * M.val * b.val * v.val;
                    val(loc.M, aa) += (M.val * v.val + Mv * steps.dt) * w * J;

                    double Av = (p.gamma_a * M.val * b.val - p.gamma_oA * A.val) * v.val - p.chi_aA * grad_Av;
                    val(loc.A, aa) += (A.val * v.val + Av * steps.dt) * w * J;

                    // TAF
                    double c_src = o < p.o_death_TC ? b.val * (1 - c.val): 0;
                    double cv = - p.diff_c * grad_dot(c, v) + c_src - p.cons_c * o * c.val;
                    val(loc.c, aa) += (c.val * v.val + cv * steps.dt) * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(now.b, loc.b, e);
                update_global_rhs(now.n, loc.n, e);
                update_global_rhs(now.f, loc.f, e);
                update_global_rhs(now.m, loc.m, e);
                update_global_rhs(now.M, loc.M, e);
                update_global_rhs(now.A, loc.A, e);
                update_global_rhs(now.c, loc.c, e);
            });
        });
    }

    void update_vasculature(int iter) {
        int next_iter = iter + 1;
        if (next_iter % vasc_update_every == 0) {
            auto taf = [&,this](double x, double y) {
                return bspline::eval_ders(x, y, now.c, this->x.B, this->y.B, xdctx, ydctx);
            };

            auto tumor = [&,this](double x, double y) {
                return bspline::eval(x, y, now.b, this->x.B, this->y.B, xctx, yctx);
            };

            vasculature.update(tumor, taf, steps.dt * vasc_update_every);

            vasculature.discretize();
        }
    }

    double& val(vector_type& v, index_type idx) {
        return v(idx[0], idx[1]);
    }
};


}
}



#endif /* ADS_PROBLEMS_TUMOR_TUMOR_HPP_ */
