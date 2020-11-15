#ifndef PROBLEMS_ELASTICITY_POURIA_HPP_
#define PROBLEMS_ELASTICITY_POURIA_HPP_

#include <cmath>
#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"


namespace problems {

    class elasticity_pouria : public ads::simulation_2d {
        using Base = simulation_2d;
        using tensor = double[2][2];

        struct state {
            vector_type ux, uy;
            vector_type vx, vy;
            state(std::array<std::size_t, 2> shape)
            : ux{ shape }, uy{ shape }
            , vx{ shape }, vy{ shape }
            { }
        };

        state now, prev, pprev;
        vector_type energy;

        ads::output_manager<2> output;

        ads::galois_executor executor{8};

        ads::lin::band_matrix Dx, Dy;
        ads::lin::band_matrix Kx, Ky;

        static constexpr double lambda = 1;
        static constexpr double mi = 1;

        int save_every = 1;

        template <typename Fun>
        void for_all(state& s, Fun fun) {
            fun(s.ux);
            fun(s.uy);
            fun(s.vx);
            fun(s.vy);
        }

    public:
        elasticity_pouria(const ads::config_2d& config, int save_every)
        : Base{ config }
        , now{ shape() }, prev{ shape() }, pprev{ shape() }
        , energy{ shape() }
        , output{ x.B, y.B, 50 }
        , Dx{x.p, x.p, x.B.dofs()}
        , Dy{y.p, y.p, y.B.dofs()}
        , Kx{x.p, x.p, x.B.dofs()}
        , Ky{y.p, y.p, y.B.dofs()}
        , save_every{save_every}
        {
            double hh = steps.dt * steps.dt * 1.0 / 3 * mi;
            matrix(Kx, x.basis, hh);
            matrix(Ky, y.basis, hh);

            double h = steps.dt * steps.dt * 1.0 / 3 * (lambda + 2 * mi);
            matrix(Dx, x.basis, h);
            matrix(Dy, y.basis, h);
        }

    private:

        void matrix(ads::lin::band_matrix& K, const ads::basis_data& d, double h) {
            for (ads::element_id e = 0; e < d.elements; ++ e) {
                for (int q = 0; q < d.quad_order; ++ q) {
                    int first = d.first_dof(e);
                    int last = d.last_dof(e);
                    for (int a = 0; a + first <= last; ++ a) {
                        for (int b = 0; b + first <= last; ++ b) {
                            int ia = a + first;
                            int ib = b + first;
                            auto va = d.b[e][q][0][a];
                            auto vb = d.b[e][q][0][b];
                            auto da = d.b[e][q][1][a];
                            auto db = d.b[e][q][1][b];
                            K(ia, ib) += (va * vb + h * da * db) * d.w[q] * d.J[e];
                        }
                    }
                }
            }
        }


        void before() override {
            prepare_matrices();
            ads::lin::factorize(Kx, x.ctx);
            ads::lin::factorize(Ky, y.ctx);
            ads::lin::factorize(Dx, x.ctx);
            ads::lin::factorize(Dy, y.ctx);
        }

        void compute_rhs(double t) {
            for_all(now, [](vector_type& a) { zero(a); });
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution(e, t);
                executor.synchronized([&] {
                    apply_local_contribution(local, e);
                });
            });
        }

        state local_contribution(index_type e, double t) const {
            auto local = state{ local_shape() };

            double J = jacobian(e);
            for (auto q : quad_points()) {
                auto x = point(e, q);
                double w = weigth(q);
                value_type vx = eval_fun(prev.vx, e, q);
                value_type vy = eval_fun(prev.vy, e, q);

                value_type ux = eval_fun(prev.ux, e, q);
                value_type uy = eval_fun(prev.uy, e, q);

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,        },
                };
                tensor s{};
                stress_tensor(s, eps);
                auto F = force(x, t);

                for (auto a : dofs_on_element(e)) {
                    value_type b = eval_basis(e, q, a);

                    double rho = 1;
                    double axb = (-s[0][0] * b.dx - s[0][1] * b.dy + F[0] * b.val) / rho;
                    double ayb = (-s[1][0] * b.dx - s[1][1] * b.dy + F[1] * b.val) / rho;

                    double dt = steps.dt;
                    double t2 = dt * dt / 2;
                    auto aa = dof_global_to_local(e, a);
                    ref(local.ux, aa) += ((ux.val + dt * vx.val) * b.val + t2 * axb) * w * J;
                    ref(local.uy, aa) += ((uy.val + dt * vy.val) * b.val + t2 * ayb) * w * J;

                    ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                    ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                }
            }
            return local;
        }

        state local_contribution_x(index_type e, double t) const {
            auto local = state{ local_shape() };
            double dt = steps.dt / 3.0;

            double J = jacobian(e);
            for (auto q : quad_points()) {
                auto x = point(e, q);
                double w = weigth(q);
                value_type vx = eval_fun(prev.vx, e, q);
                value_type vy = eval_fun(prev.vy, e, q);

                value_type ux = eval_fun(prev.ux, e, q) + dt * vx;
                value_type uy = eval_fun(prev.uy, e, q) + dt * vy;

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,        },
                };
                tensor s{};
                stress_tensor(s, eps);
                auto F = force(x, t);

                for (auto a : dofs_on_element(e)) {
                    value_type b = eval_basis(e, q, a);

                    double rho = 1;
                    double axb = (-s[0][0] * b.dx - s[0][1] * b.dy + F[0] * b.val) / rho;
                    double ayb = (-s[1][0] * b.dx - s[1][1] * b.dy + F[1] * b.val) / rho;

                    axb += (lambda + 2 * mi) * dt * vx.dx * b.dx;
                    ayb +=                mi * dt * vy.dx * b.dx;

                    auto aa = dof_global_to_local(e, a);
                    ref(local.ux, aa) += ux.val * b.val * w * J;
                    ref(local.uy, aa) += uy.val * b.val * w * J;

                    ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                    ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                }
            }
            return local;
        }
        state local_contribution_y(index_type e, double t) const {
            auto local = state{ local_shape() };
            double dt = steps.dt / 3.0;

            double J = jacobian(e);
            for (auto q : quad_points()) {
                auto x = point(e, q);
                double w = weigth(q);
                value_type vx = eval_fun(prev.vx, e, q);
                value_type vy = eval_fun(prev.vy, e, q);

                value_type ux = eval_fun(prev.ux, e, q) + dt * vx;
                value_type uy = eval_fun(prev.uy, e, q) + dt * vy;

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,        },
                };
                tensor s{};
                stress_tensor(s, eps);
                auto F = force(x, t);

                for (auto a : dofs_on_element(e)) {
                    value_type b = eval_basis(e, q, a);

                    double rho = 1;
                    double axb = (-s[0][0] * b.dx - s[0][1] * b.dy + F[0] * b.val) / rho;
                    double ayb = (-s[1][0] * b.dx - s[1][1] * b.dy + F[1] * b.val) / rho;

                    axb +=                mi * dt * vx.dx * b.dx;
                    ayb += (lambda + 2 * mi) * dt * vy.dx * b.dx;

                    auto aa = dof_global_to_local(e, a);
                    ref(local.ux, aa) += ux.val * b.val * w * J;
                    ref(local.uy, aa) += uy.val * b.val * w * J;

                    ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                    ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                }
            }
            return local;
        }

        void apply_local_contribution(const state& loc, index_type e) {
            update_global_rhs(now.ux, loc.ux, e);
            update_global_rhs(now.uy, loc.uy, e);
            update_global_rhs(now.vx, loc.vx, e);
            update_global_rhs(now.vy, loc.vy, e);
        }

        double kinetic_energy() const {
            double E = 0;
            executor.for_each(elements(), [&](index_type e) {
                double Eloc = 0;
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    double w = weigth(q);
                    value_type vx = eval_fun(now.vx, e, q);
                    value_type vy = eval_fun(now.vy, e, q);
                    Eloc += w * J * 0.5 * (vx.val * vx.val + vy.val * vy.val);
                }
                executor.synchronized([&] { E += Eloc;  });
            });
            return E;
        }

        double potential_energy() const {
            double E = 0;
            executor.for_each(elements(), [&](index_type e) {
                double Eloc = 0;
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    double w = weigth(q);
                    value_type ux = eval_fun(now.ux, e, q);
                    value_type uy = eval_fun(now.uy, e, q);

                    tensor eps = {
                        {         ux.dx,         0.5 * (ux.dy + uy.dx) },
                        { 0.5 * (ux.dy + uy.dx),         uy.dy,        },
                    };
                    tensor s{};
                    stress_tensor(s, eps);
                    double U = 0;
                    for (int i = 0; i < 2; ++ i) {
                        for (int j = 0; j < 2; ++ j) {
                            U += 0.5 * s[i][j] * eps[i][j];
                        }
                    }
                    Eloc += w * J * U;
                }
                executor.synchronized([&] { E += Eloc;  });
            });
            return E;
        }

        void compute_potential_energy() {
            zero(energy);
            executor.for_each(elements(), [&](index_type e) {
                auto Eloc = element_rhs();
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    double w = weigth(q);
                    value_type ux = eval_fun(now.ux, e, q);
                    value_type uy = eval_fun(now.uy, e, q);

                    tensor eps = {
                        {         ux.dx,         0.5 * (ux.dy + uy.dx) },
                        { 0.5 * (ux.dy + uy.dx),         uy.dy,        },
                    };
                    tensor s{};
                    stress_tensor(s, eps);
                    double U = 0;
                    for (int i = 0; i < 2; ++ i) {
                        for (int j = 0; j < 2; ++ j) {
                            U += 0.5 * s[i][j] * eps[i][j];
                        }
                    }
                    for (auto a : dofs_on_element(e)) {
                        auto aa = dof_global_to_local(e, a);
                        value_type v = eval_basis(e, q, a);
                        Eloc(aa[0], aa[1]) += w * J * v.val * U * 1e8;
                    }
                }
                executor.synchronized([&]() {
                    update_global_rhs(energy, Eloc, e);
                });
            });
            solve(energy);
        }

        double total() const {
            double E = 0;
            for (auto e : elements()) {
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    double w = weigth(q);
                    value_type vx = eval_fun(now.vx, e, q);
                    E += w * J * vx.val;
                }
            }
            return E;
        }

        void stress_tensor(tensor& s, const tensor& eps) const {
            double tr = eps[0][0] + eps[1][1];

            for (int i = 0; i < 2; ++ i) {
                for (int j = 0; j < 2; ++ j) {
                    s[i][j] = 2 * mi * eps[i][j];
                }
                s[i][i] += lambda * tr;
            }
        }

        point_type force(point_type x, double t) const {
            using std::pow;
            constexpr double t0 = 1;
            double tt = t / t0;
            double f = tt < 1 ? pow(tt * (1 - tt), 2) : 0;
            double r = pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
            double a = - 10 * f * std::exp(- 10 * r);
            return {a, a};
        }

        void before_step(int /*iter*/, double /*t*/) override {
            // using std::swap;
            // swap(now, prev);

            pprev = prev;
            prev  = now;
        }

        void step2(int /*iter*/, double t) {
            compute_rhs(t);
            for_all(now, [this](vector_type& a) { solve(a); });
        }

        void step(int /*iter*/, double t) {
            using std::swap;

            for_all(now, [](vector_type& a) { zero(a); });
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution_x(e, t);
                executor.synchronized([&] {
                    apply_local_contribution(local, e);
                });
            });

            ads_solve(now.ux, buffer, x.data(), y.data());
            ads_solve(now.uy, buffer, x.data(), y.data());

            ads_solve(now.vx, buffer, ads::dim_data{Dx, x.ctx}, y.data());
            ads_solve(now.vy, buffer, ads::dim_data{Kx, x.ctx}, y.data());

            swap(now, prev);

            for_all(now, [](vector_type& a) { zero(a); });
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution_y(e, t);
                executor.synchronized([&] {
                    apply_local_contribution(local, e);
                });
            });

            ads_solve(now.ux, buffer, x.data(), y.data());
            ads_solve(now.uy, buffer, x.data(), y.data());

            ads_solve(now.vx, buffer, x.data(), ads::dim_data{Ky, y.ctx});
            ads_solve(now.vy, buffer, x.data(), ads::dim_data{Dy, y.ctx});
        }

        void after_step(int iter, double t) override {
            iter += 1;
            if (iter % save_every == 0) {
                std::cout << "** Iteration " << iter << ", t = " << t + steps.dt << std::endl;

                double Ek = kinetic_energy();
                double Ep = potential_energy();
                compute_potential_energy();

                // std::cout << "Kinetic energy: " << Ek << std::endl;
                // std::cout << "Potential energy: " << Ep << std::endl;
                // std::cout << "Total energy: " << Ek + Ep << std::endl;

                // std::cout << "Total disp:   : " << total() << std::endl;
                // std::cout << std::endl;

                output.to_file(now.ux, "ux_%d.data", iter);
                output.to_file(now.uy, "uy_%d.data", iter);
                output.to_file(energy, "energy_%d.data", iter);

                // output.to_file("out_%d.data", iter*10,
                //                output.evaluate(now.ux),
                //                output.evaluate(now.uy),
                //                output.evaluate(now.uz),
                //                output.evaluate(energy));

                std::cout << iter << " " << t << " " << Ek << " " << Ep << " " << Ek + Ep << std::endl;
                // std::cout << "Step " << (iter + 1) << ", t = " << t << std::endl;

                // int num = (iter + 1) / save_every;
                // auto name = str(boost::format("out_%d.data") % num);
                // std::ofstream sol(name);
                // for (int i = 0; i < x.dofs(); ++ i) {
                //     for (int j = 0; j < y.dofs(); ++ j) {
                //         for (int k = 0; k < z.dofs(); ++ k) {
                //             sol << i << " " << j << " " << " " << k << " "
                //                 << now.ux(i, j, k) << " "
                //                 << now.uy(i, j, k) << " "
                //                 << now.uz(i, j, k) << std::endl;
                //         }
                //     }
                // }
            }
        }


        double& ref(vector_type& v, index_type idx) const {
            return v(idx[0], idx[1]);
        }
    };
}



#endif /* PROBLEMS_ELASTICITY_POURIA_HPP_ */
