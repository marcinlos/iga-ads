#ifndef ADS_PROBLEMS_TUMOR_3D_TUMOR_3D_HPP
#define ADS_PROBLEMS_TUMOR_3D_TUMOR_3D_HPP

#include <algorithm>

#include "problems/tumor/skin.hpp"
#include "problems/tumor/state.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/params.hpp"

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"

#include "ads/executor/galois.hpp"


namespace tumor {

    class tumor_3d : public ads::simulation_3d {
    private:
        static constexpr std::size_t Dim = 3;
        using Base = ads::simulation_3d;

        state<Dim> now, prev;
        params p;

        ads::output_manager<3> output;

        ads::galois_executor executor{4};

    public:
        tumor_3d(const ads::config_3d& config, const params& params)
        : Base{config}
        , now{ shape() }
        , prev{ shape() }
        , p{ params }
        , output{ x.B, y.B, z.B, 50 }
        { }

    private:

        auto constant(double c) const {
            return [c](double, double, double) { return c; };
        }

        void dirichlet() {
            x.fix_left();
            x.fix_right();
            y.fix_left();
            y.fix_right();
            z.fix_left();
            z.fix_right();
        }

        void zero_bc(vector_type& v) {
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < y.dofs(); ++ j) {
                    v(i, j, 0) = 0;
                    v(i, j, z.dofs() - 1) = 0;
                }
            }
            for (int i = 0; i < y.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(0, i, j) = 0;
                    v(x.dofs() - 1, i, j) = 0;
                }
            }
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(i, 0, j) = 0;
                    v(i, y.dofs() - 1, j) = 0;
                }
            }
        }

        template <typename Function>
        void apply_bc(vector_type& v, Function f) const {
            using side_type = ads::lin::tensor<double, 2>;
            side_type xy{{ x.dofs(), y.dofs() }};

            ads::compute_projection(xy, x.basis, y.basis, [&](double x, double y) { return f(x, y, z.a); });
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < y.dofs(); ++ j) {
                    v(i, j, 0) = xy(i, j);
                }
            }

            zero(xy);
            ads::compute_projection(xy, x.basis, y.basis, [&](double x, double y) { return f(x, y, z.b); });
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < y.dofs(); ++ j) {
                    v(i, j, z.dofs() - 1) = xy(i, j);
                }
            }

            side_type yz{{ y.dofs(), z.dofs() }};

            ads::compute_projection(yz, y.basis, z.basis, [&](double y, double z) { return f(x.a, y, z); });
            for (int i = 0; i < y.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(0, i, j) = yz(i, j);
                }
            }

            zero(yz);
            ads::compute_projection(yz, y.basis, z.basis, [&](double y, double z) { return f(x.b, y, z); });
            for (int i = 0; i < y.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(x.dofs() - 1, i, j) = yz(i, j);
                }
            }

            side_type xz{{ x.dofs(), z.dofs() }};

            ads::compute_projection(xz, x.basis, z.basis, [&](double x, double z) { return f(x, y.a, z); });
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(i, 0, j) = yz(i, j);
                }
            }

            zero(xz);
            ads::compute_projection(xz, x.basis, z.basis, [&](double x, double z) { return f(x, y.b, z); });
            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < z.dofs(); ++ j) {
                    v(i, y.dofs() - 1, j) = yz(i, j);
                }
            }

            using edge_type = ads::lin::tensor<double, 1>;
            edge_type ex{{ x.dofs() }}, ey{{ y.dofs() }}, ez{{ z.dofs() }};

            ads::compute_projection(ex, x.basis, [&](double t) { return f(t, y.a, z.a); });
            for (int i = 0; i < x.dofs(); ++ i) {
                v(i, 0, 0) = ex(i);
            }
            zero(ex);

            ads::compute_projection(ex, x.basis, [&](double t) { return f(t, y.b, z.a); });
            for (int i = 0; i < x.dofs(); ++ i) {
                v(i, y.dofs() - 1, 0) = ex(i);
            }
            zero(ex);

            ads::compute_projection(ex, x.basis, [&](double t) { return f(t, y.a, z.b); });
            for (int i = 0; i < x.dofs(); ++ i) {
                v(i, 0, z.dofs() - 1) = ex(i);
            }
            zero(ex);

            ads::compute_projection(ex, x.basis, [&](double t) { return f(t, y.b, z.b); });
            for (int i = 0; i < x.dofs(); ++ i) {
                v(i, y.dofs() - 1, z.dofs() - 1) = ex(i);
            }
            zero(ex);

            ads::compute_projection(ey, y.basis, [&](double t) { return f(x.a, t, z.a); });
            for (int i = 0; i < y.dofs(); ++ i) {
                v(0, i, 0) = ey(i);
            }
            zero(ey);

            ads::compute_projection(ey, y.basis, [&](double t) { return f(x.b, t, z.a); });
            for (int i = 0; i < y.dofs(); ++ i) {
                v(x.dofs() - 1, i, 0) = ey(i);
            }
            zero(ey);

            ads::compute_projection(ey, y.basis, [&](double t) { return f(x.a, t, z.b); });
            for (int i = 0; i < y.dofs(); ++ i) {
                v(0, i, z.dofs() - 1) = ey(i);
            }
            zero(ey);

            ads::compute_projection(ey, y.basis, [&](double t) { return f(x.b, t, z.b); });
            for (int i = 0; i < y.dofs(); ++ i) {
                v(x.dofs() - 1, i, z.dofs() - 1) = ey(i);
            }
            zero(ey);

            ads::compute_projection(ez, z.basis, [&](double t) { return f(x.a, y.a, t); });
            for (int i = 0; i < z.dofs(); ++ i) {
                v(0, 0, i) = ez(i);
            }
            zero(ez);

            ads::compute_projection(ez, z.basis, [&](double t) { return f(x.b, y.a, t); });
            for (int i = 0; i < z.dofs(); ++ i) {
                v(x.dofs() - 1, 0, i) = ez(i);
            }
            zero(ez);

            ads::compute_projection(ez, z.basis, [&](double t) { return f(x.a, y.b, t); });
            for (int i = 0; i < z.dofs(); ++ i) {
                v(0, y.dofs() - 1, i) = ez(i);
            }
            zero(ez);

            ads::compute_projection(ez, z.basis, [&](double t) { return f(x.b, y.b, t); });
            for (int i = 0; i < z.dofs(); ++ i) {
                v(x.dofs() - 1, y.dofs() - 1, i) = ez(i);
            }
            zero(ez);

            v(0, 0, 0)                                  = f(x.a, y.a, z.a);
            v(0, 0, z.dofs() - 1)                       = f(x.a, y.a, z.b);
            v(0, y.dofs() - 1, 0)                       = f(x.a, y.b, z.a);
            v(0, y.dofs() - 1, z.dofs() - 1)            = f(x.a, y.b, z.b);
            v(x.dofs() - 1, 0, 0)                       = f(x.b, y.a, z.a);
            v(x.dofs() - 1, 0, z.dofs() - 1)            = f(x.b, y.a, z.b);
            v(x.dofs() - 1, y.dofs() - 1, 0)            = f(x.b, y.b, z.a);
            v(x.dofs() - 1, y.dofs() - 1, z.dofs() - 1) = f(x.b, y.b, z.b);
        }

        void prepare_matrices() {
            dirichlet();
            Base::prepare_matrices();
        }

        void before() override {
            prepare_matrices();

            projection(now.b, constant(0));
            projection(now.c, constant(0));
            projection(now.o, constant(0));

            projection(now.A, constant(0));
            projection(now.M, [this](double x, double y, double z) { return p.skin.init_M(x, y, z); });

            solve_all();

            save_to_file(0);
        }

        void step(int /*iter*/, double /*t*/) override {
            compute_rhs();
            solve_all();
        }


        void solve_all() {
            zero_bc(now.b);
            solve(now.b);

            zero_bc(now.c);
            solve(now.c);

            zero_bc(now.o);
            solve(now.o);

            apply_bc(now.M, [this](double x, double y, double z) { return p.skin.init_M(x, y, z); });
            solve(now.M);

            zero_bc(now.A);
            solve(now.A);
        }

        void compute_rhs() {
            now.clear();
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution(e);
                executor.synchronized([&] {
                    apply_local_contribution(local, e);
                });
            });
        }

        state<Dim> local_contribution(index_type e) const {
            auto local = state<Dim>{ local_shape() };

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    value_type b = ensure_positive(eval_fun(prev.b, e, q));
                    value_type c = ensure_positive(eval_fun(prev.c, e, q));
                    value_type o = ensure_positive(eval_fun(prev.o, e, q));

                    value_type M = ensure_positive(eval_fun(prev.M, e, q));
                    value_type A = ensure_positive(eval_fun(prev.A, e, q));

                    // tumor density
                    double b_src = 0;
                    double b_sink = 0;

                    // tumor source/sink
                    if (o.val >= p.o_prol_TC) {
                        double tbA = p.tau_b * A.val;
                        double s = 1 + tbA / (1 + tbA) * p.P_b;
                        b_src = b.val / p.t_prol_TC * s * (1 - b.val / p.c_b_max);
                    }
                    if (o.val < p.o_death_TC) {
                        b_sink = -b.val / p.t_death_TC;
                    }

                    double D_b = p.skin.diffusion(x[0], x[1], x[2]);
                    double grad_Av = grad_dot(A, v);
                    double grad_Pv = b.val >= p.c_b_norm ? grad_dot(b, v) / (p.c_b_max - p.c_b_norm) : 0;

                    double divJv = D_b * b.val * (grad_Pv + p.r_b * grad_Av);

                    double bv = - divJv + (b_src + b_sink) * v.val;
                    ref(local.b, aa) += (b.val * v.val + bv * steps.dt) * w * J;

                    // ECM evolution
                    double Mv = - p.beta_m * M.val * b.val * v.val;
                    ref(local.M, aa) += (M.val * v.val + Mv * steps.dt) * w * J;

                    double Av = (p.gamma_a * M.val * b.val - p.gamma_oA * A.val) * v.val - p.chi_aA * grad_Av;
                    ref(local.A, aa) += (A.val * v.val + Av * steps.dt) * w * J;

                    // TAF
                    double c_src = o.val < p.o_death_TC ? b.val * (1 - c.val): 0;
                    double cv = - p.diff_c * grad_dot(c, v) + c_src - p.cons_c * o.val * c.val;
                    ref(local.c, aa) += (c.val * v.val + cv * steps.dt) * w * J;

                    // oxygen
                    using std::max;
                    double o2src = 1.0;
                    double ov = - p.alpha_0 * grad_dot(o, v) + (- p.gamma_T * b.val * o.val + p.alpha_1 * max(0.0, 60 - o.val) * o2src) * v.val;
                    ref(local.o, aa) += (o.val * v.val + ov * steps.dt) * w * J;
                }
            }
            return local;
        }


        void apply_local_contribution(const state<Dim>& loc, index_type e) {
            update_global_rhs(now.b, loc.b, e);
            update_global_rhs(now.c, loc.c, e);
            update_global_rhs(now.o, loc.o, e);
            update_global_rhs(now.A, loc.A, e);
            update_global_rhs(now.M, loc.M, e);
        }

        double& ref(vector_type& v, index_type idx) const {
            return v(idx[0], idx[1], idx[2]);
        }

        value_type ensure_positive(value_type v) const {
            v.val = std::max(v.val, 0.0);
            return v;
        }

        void save_to_file(int iter) {
            output.to_file(now.b, "tumor_%d.vti", iter);
            output.to_file(now.c, "taf_%d.vti", iter);
            output.to_file(now.o, "oxygen_%d.vti", iter);
            output.to_file(now.M, "ecm_%d.vti", iter);
            output.to_file(now.A, "degraded_ecm_%d.vti", iter);
        }

    };

}


#endif /* ADS_PROBLEMS_TUMOR_3D_TUMOR_3D_HPP */
