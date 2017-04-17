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

#include "problems/tumor/3d/vasculature.hpp"


namespace tumor {

    class tumor_3d : public ads::simulation_3d {
    private:
        static constexpr std::size_t Dim = 3;
        using Base = ads::simulation_3d;

        state<Dim> now, prev;
        state<Dim> k1, k2, k3, k4;
        params p;
        vasculature vasc;

        ads::bspline::eval_ctx xctx;
        ads::bspline::eval_ctx yctx;
        ads::bspline::eval_ctx zctx;

        ads::bspline::eval_ders_ctx xdctx;
        ads::bspline::eval_ders_ctx ydctx;
        ads::bspline::eval_ders_ctx zdctx;

        ads::output_manager<3> output;

        ads::galois_executor executor;

    public:
        tumor_3d(const ads::config_3d& config, const params& params, vasculature vasc, int threads)
        : Base{config}
        , now{ shape() }
        , prev{ shape() }
        , k1{ shape() }, k2{ shape() }, k3{ shape() }, k4{ shape() }
        , p{ params }
        , vasc{ std::move(vasc) }
        , xctx{ x.B.degree }
        , yctx{ y.B.degree }
        , zctx{ z.B.degree }
        , xdctx{ x.B.degree, 1 }
        , ydctx{ y.B.degree, 1 }
        , zdctx{ z.B.degree, 1}
        , output{ x.B, y.B, z.B, 50 }
        , executor{ threads }
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

            projection(now.b, [](double x, double y, double z) {
                double dx = (x - 2500) / 400;
                double dy = (y - 2500) / 400;
                double dz = (z - 2440) / 400;
                double r2 = std::min(1.0, dx*dx + dy*dy + dz*dz);
                return 1.2 * (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
            });

            projection(now.c, constant(0));
            projection(now.o, constant(0));

            projection(now.A, constant(0));
            projection(now.M, [this](double x, double y, double z) { return p.skin.init_M(x, y, z); });

            solve_all(now);

            save_to_file(0);
        }

        void before_step(int /*iter*/, double /*t*/) override {
            using std::swap;
            swap(now, prev);

            now.clear();
        }

        void step(int iter, double /*t*/) override {
            double h = steps.dt;
            rk_step(prev, k1, h/2);
            rk_step(k1, k2, h/2);
            rk_step(k2, k3, h);
            rk_step(k3, k4, -h/2);

            for (int i = 0; i < x.dofs(); ++ i) {
                for (int j = 0; j < y.dofs(); ++ j) {
                    for (int k = 0; k < z.dofs(); ++ k) {
                        now.b(i, j, k) = 1.0 / 3 * (k1.b(i, j, k) + 2 * k2.b(i, j, k) + k3.b(i, j, k) - k4.b(i, j, k));
                        now.c(i, j, k) = 1.0 / 3 * (k1.c(i, j, k) + 2 * k2.c(i, j, k) + k3.c(i, j, k) - k4.c(i, j, k));
                        now.o(i, j, k) = 1.0 / 3 * (k1.o(i, j, k) + 2 * k2.o(i, j, k) + k3.o(i, j, k) - k4.o(i, j, k));
                        now.M(i, j, k) = 1.0 / 3 * (k1.M(i, j, k) + 2 * k2.M(i, j, k) + k3.M(i, j, k) - k4.M(i, j, k));
                        now.A(i, j, k) = 1.0 / 3 * (k1.A(i, j, k) + 2 * k2.A(i, j, k) + k3.A(i, j, k) - k4.A(i, j, k));
                    }
                }
            }

            update_vasculature(iter);
        }

        void rk_step(const state<Dim>& prev, state<Dim>& next, double h) {
            next.clear();
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution(prev, e, h);
                executor.synchronized([&] { apply_local_contribution(next, local, e); });
            });
            solve_all(next);
        }

        void after_step(int iter, double /*t*/) override {
            std::cout << "Iter " << iter << " done" << std::endl;
            if ((iter + 1) % 100 == 0) {
                save_to_file(iter + 1);
            }
        }


        void solve_all(state<Dim>& s) {
            zero_bc(s.b);
            solve(s.b);

            zero_bc(s.c);
            solve(s.c);

            zero_bc(s.o);
            solve(s.o);

            apply_bc(s.M, [this](double x, double y, double z) { return p.skin.init_M(x, y, z); });
            solve(s.M);

            zero_bc(s.A);
            solve(s.A);
        }

        state<Dim> local_contribution(const state<Dim>& s, index_type e, double h) const {
            auto local = state<Dim>{ local_shape() };

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    value_type b = ensure_positive(eval_fun(s.b, e, q));
                    value_type c = ensure_positive(eval_fun(s.c, e, q));
                    value_type o = ensure_positive(eval_fun(s.o, e, q));

                    value_type M = ensure_positive(eval_fun(s.M, e, q));
                    value_type A = ensure_positive(eval_fun(s.A, e, q));

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
                    ref(local.b, aa) += (b.val * v.val + bv * h) * w * J;

                    // ECM evolution
                    double Mv = - p.beta_m * M.val * b.val * v.val;
                    ref(local.M, aa) += (M.val * v.val + Mv * h) * w * J;

                    double Av = (p.gamma_a * M.val * b.val - p.gamma_oA * A.val) * v.val - p.chi_aA * grad_Av;
                    ref(local.A, aa) += (A.val * v.val + Av * h) * w * J;

                    // TAF
                    double c_src = o.val < p.o_death_TC ? b.val * (1 - c.val) : 0;
                    double cv = - p.diff_c * grad_dot(c, v) + (c_src - p.cons_c * c.val * o.val) * v.val;
                    ref(local.c, aa) += (c.val * v.val + cv * h) * w * J;

                    // oxygen
                    using std::max;
                    double o2src = oxygen(x);
                    double o_rhs = - p.gamma_T * b.val * o.val + p.alpha_1 * max(0.0, p.o_max - o.val) * o2src;
                    double ov = - p.alpha_0 * grad_dot(o, v) + o_rhs * v.val;
                    ref(local.o, aa) += (o.val * v.val + ov * h) * w * J;
                }
            }
            return local;
        }

        double oxygen(point_type p) const {
            double xx = (p[0] - x.a) / (x.b - x.a);
            double yy = (p[1] - y.a) / (y.b - y.a);
            double zz = (p[2] - z.a) / (z.b - z.a);
            return vasc.source(xx, yy, zz);
        }

        void apply_local_contribution(state<Dim>& rhs, const state<Dim>& loc, index_type e) const {
            update_global_rhs(rhs.b, loc.b, e);
            update_global_rhs(rhs.c, loc.c, e);
            update_global_rhs(rhs.o, loc.o, e);
            update_global_rhs(rhs.A, loc.A, e);
            update_global_rhs(rhs.M, loc.M, e);
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

            using boost::format;
            vasc.to_file(str(format("vasculature_%d.vti") % iter));
        }

        void update_vasculature(int iter) {
            auto taf = [&,this](double px, double py, double pz) {
                double xx = x.a + px * (x.b - x.a);
                double yy = y.a + py * (y.b - y.a);
                double zz = z.a + pz * (z.b - z.a);
                return ads::bspline::eval_ders(xx, yy, zz, now.c, this->x.B, this->y.B, this->z.B, xdctx, ydctx, zdctx);
            };

            auto tumor = [&,this](double px, double py, double pz) {
                double xx = x.a + px * (x.b - x.a);
                double yy = y.a + py * (y.b - y.a);
                double zz = z.a + pz * (z.b - z.a);
                return ads::bspline::eval(xx, yy, zz, now.b, this->x.B, this->y.B, this->z.B, xctx, yctx, zctx);
            };
            vasc.update(tumor, taf, iter, steps.dt);
        }
    };

}


#endif /* ADS_PROBLEMS_TUMOR_3D_TUMOR_3D_HPP */
