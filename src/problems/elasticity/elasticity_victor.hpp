#ifndef PROBLEMS_ELASTICITY_ELASTICITY_VICTOR_HPP_
#define PROBLEMS_ELASTICITY_ELASTICITY_VICTOR_HPP_

#include <cmath>
#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"


namespace problems {

class elasticity_victor : public ads::simulation_3d {
    using Base = simulation_3d;
    using tensor = double[3][3];

    struct state {
        vector_type ux, uy, uz;
        vector_type vx, vy, vz;
        vector_type ax, ay, az;

        state(std::array<std::size_t, 3> shape)
            : ux{ shape }, uy{ shape }, uz{ shape }
            , vx{ shape }, vy{ shape }, vz{ shape }
            , ax{ shape }, ay{ shape }, az{ shape }
            { }
    };

    state now, prev;
    vector_type energy;

    ads::output_manager<3> output;

    ads::galois_executor executor{8};

    ads::lin::band_matrix Kx, Ky, Kz;

    static constexpr double rho = 1;
    static constexpr double lambda = -1;
    static constexpr double mi = 1;

    int save_every = 1;

    template <typename Fun>
    void for_all(state& s, Fun fun) {
        fun(s.ux);
        fun(s.uy);
        fun(s.uz);
        fun(s.vx);
        fun(s.vy);
        fun(s.vz);
        fun(s.ax);
        fun(s.ay);
        fun(s.az);
    }

public:
    elasticity_victor(const ads::config_3d& config, int save_every)
        : Base{ config }
        , now{ shape() }, prev{ shape() }
        , energy{ shape() }
        , output{ x.B, y.B, z.B, 50 }
        , Kx{x.p, x.p, x.B.dofs()}
        , Ky{y.p, y.p, y.B.dofs()}
        , Kz{z.p, z.p, z.B.dofs()}
        , save_every{save_every}
        {
            auto dt = steps.dt / 3;
            double hh = 0.5 * dt * dt;
            matrix(Kx, x.basis, hh / rho);
            matrix(Ky, y.basis, hh / rho);
            matrix(Kz, z.basis, hh / rho);
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

    void add(vector_type& left, const vector_type& right, double a) const {
        for (auto ix = 0; ix < x.dofs(); ++ ix) {
            for (auto iy = 0; iy < y.dofs(); ++ iy) {
                for (auto iz = 0; iz < z.dofs(); ++ iz) {
                    left(ix, iy, iz) += a * right(ix, iy, iz);
                }
            }
        }
    }

    void before() override {
        prepare_matrices();
        ads::lin::factorize(Kx, x.ctx);
        ads::lin::factorize(Ky, y.ctx);
        ads::lin::factorize(Kz, z.ctx);
    }

    state local_contribution_x(index_type e, double t) const {
        auto local = state{ local_shape() };
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weigth(q);
            value_type ux = eval_fun(prev.ux, e, q);
            value_type uy = eval_fun(prev.uy, e, q);
            value_type uz = eval_fun(prev.uz, e, q);

            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ax = eval_fun(prev.ax, e, q);
            value_type ay = eval_fun(prev.ay, e, q);
            value_type az = eval_fun(prev.az, e, q);

            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double tt = 0.5 * dt * dt;
                auto dax = - tt * (ax.dy * b.dy + ax.dz * b.dz) - mi * grad_dot(ux, b) - dt * grad_dot(vx, b) +  F[0] * b.val;
                auto day = - tt * (ay.dy * b.dy + ay.dz * b.dz) - mi * grad_dot(uy, b) - dt * grad_dot(vy, b) +  F[1] * b.val;
                auto daz = - tt * (az.dy * b.dy + az.dz * b.dz) - mi * grad_dot(uz, b) - dt * grad_dot(vz, b) +  F[2] * b.val;

                // Version 1
                // dax -= (mi + lambda) * (ux.dx * b.dx + uy.dx * b.dy + uz.dx * b.dz);
                // day -= (mi + lambda) * (ux.dy * b.dx + uy.dy * b.dy + uz.dy * b.dz);
                // daz -= (mi + lambda) * (ux.dz * b.dx + uy.dz * b.dy + uz.dz * b.dz);

                // Version 2
                // auto tr = ux.dx + uy.dy + uz.dz;
                // dax -= (mi + lambda) * tr * b.dx;
                // day -= (mi + lambda) * tr * b.dy;
                // daz -= (mi + lambda) * tr * b.dz;

                auto aa = dof_global_to_local(e, a);
                ref(local.ax, aa) += dax / rho * w * J;
                ref(local.ay, aa) += day / rho * w * J;
                ref(local.az, aa) += daz / rho * w * J;
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
            value_type ux = eval_fun(prev.ux, e, q);
            value_type uy = eval_fun(prev.uy, e, q);
            value_type uz = eval_fun(prev.uz, e, q);

            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ax = eval_fun(prev.ax, e, q);
            value_type ay = eval_fun(prev.ay, e, q);
            value_type az = eval_fun(prev.az, e, q);

            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double tt = 0.5 * dt * dt;
                auto dax = - tt * (ax.dx * b.dx + ax.dz * b.dz) - mi * grad_dot(ux, b) - dt * grad_dot(vx, b) +  F[0] * b.val;
                auto day = - tt * (ay.dx * b.dx + ay.dz * b.dz) - mi * grad_dot(uy, b) - dt * grad_dot(vy, b) +  F[1] * b.val;
                auto daz = - tt * (az.dx * b.dx + az.dz * b.dz) - mi * grad_dot(uz, b) - dt * grad_dot(vz, b) +  F[2] * b.val;

                // Version 1
                // dax -= (mi + lambda) * (ux.dx * b.dx + uy.dx * b.dy + uz.dx * b.dz);
                // day -= (mi + lambda) * (ux.dy * b.dx + uy.dy * b.dy + uz.dy * b.dz);
                // daz -= (mi + lambda) * (ux.dz * b.dx + uy.dz * b.dy + uz.dz * b.dz);

                // Version 2
                // auto tr = ux.dx + uy.dy + uz.dz;
                // dax -= (mi + lambda) * tr * b.dx;
                // day -= (mi + lambda) * tr * b.dy;
                // daz -= (mi + lambda) * tr * b.dz;

                auto aa = dof_global_to_local(e, a);
                ref(local.ax, aa) += dax / rho * w * J;
                ref(local.ay, aa) += day / rho * w * J;
                ref(local.az, aa) += daz / rho * w * J;
            }
        }
        return local;
    }

    state local_contribution_z(index_type e, double t) const {
        auto local = state{ local_shape() };
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weigth(q);
            value_type ux = eval_fun(prev.ux, e, q);
            value_type uy = eval_fun(prev.uy, e, q);
            value_type uz = eval_fun(prev.uz, e, q);

            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ax = eval_fun(prev.ax, e, q);
            value_type ay = eval_fun(prev.ay, e, q);
            value_type az = eval_fun(prev.az, e, q);

            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double tt = 0.5 * dt * dt;
                auto dax = - tt * (ax.dx * b.dx + ax.dy * b.dy) - mi * grad_dot(ux, b) - dt * grad_dot(vx, b) +  F[0] * b.val;
                auto day = - tt * (ay.dx * b.dx + ay.dy * b.dy) - mi * grad_dot(uy, b) - dt * grad_dot(vy, b) +  F[1] * b.val;
                auto daz = - tt * (az.dx * b.dx + az.dy * b.dy) - mi * grad_dot(uz, b) - dt * grad_dot(vz, b) +  F[2] * b.val;

                // Version 1
                // dax -= (mi + lambda) * (ux.dx * b.dx + uy.dx * b.dy + uz.dx * b.dz);
                // day -= (mi + lambda) * (ux.dy * b.dx + uy.dy * b.dy + uz.dy * b.dz);
                // daz -= (mi + lambda) * (ux.dz * b.dx + uy.dz * b.dy + uz.dz * b.dz);

                // Version 2
                // auto tr = ux.dx + uy.dy + uz.dz;
                // dax -= (mi + lambda) * tr * b.dx;
                // day -= (mi + lambda) * tr * b.dy;
                // daz -= (mi + lambda) * tr * b.dz;

                auto aa = dof_global_to_local(e, a);
                ref(local.ax, aa) += dax / rho * w * J;
                ref(local.ay, aa) += day / rho * w * J;
                ref(local.az, aa) += daz / rho * w * J;
            }
        }
        return local;
    }

    state local_contribution_truncated(index_type e, double t) const {
        auto local = state{ local_shape() };
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weigth(q);
            value_type ux = eval_fun(prev.ux, e, q);
            value_type uy = eval_fun(prev.uy, e, q);
            value_type uz = eval_fun(prev.uz, e, q);

            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ax = eval_fun(prev.ax, e, q);
            value_type ay = eval_fun(prev.ay, e, q);
            value_type az = eval_fun(prev.az, e, q);

            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double tt = 0.5 * dt * dt;
                auto dax = - mi * grad_dot(ux, b) - dt * grad_dot(vx, b) +  F[0] * b.val;
                auto day = - mi * grad_dot(uy, b) - dt * grad_dot(vy, b) +  F[1] * b.val;
                auto daz = - mi * grad_dot(uz, b) - dt * grad_dot(vz, b) +  F[2] * b.val;

                auto aa = dof_global_to_local(e, a);
                ref(local.ax, aa) += dax / rho * w * J;
                ref(local.ay, aa) += day / rho * w * J;
                ref(local.az, aa) += daz / rho * w * J;
            }
        }
        return local;
    }

    void apply_local_contribution(const state& loc, index_type e) {
        update_global_rhs(now.ax, loc.ax, e);
        update_global_rhs(now.ay, loc.ay, e);
        update_global_rhs(now.az, loc.az, e);
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
                value_type vz = eval_fun(now.vz, e, q);
                Eloc += w * J * 0.5 * (vx.val * vx.val + vy.val * vy.val + vz.val * vz.val);
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
                value_type uz = eval_fun(now.uz, e, q);

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,         0.5 * (uy.dz + uz.dy) },
                    { 0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy),         uz.dz         }
                };
                tensor s{};
                stress_tensor(s, eps);
                double U = 0;
                for (int i = 0; i < 3; ++ i) {
                    for (int j = 0; j < 3; ++ j) {
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
                value_type uz = eval_fun(now.uz, e, q);

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,         0.5 * (uy.dz + uz.dy) },
                    { 0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy),         uz.dz         }
                };
                tensor s{};
                stress_tensor(s, eps);
                double U = 0;
                for (int i = 0; i < 3; ++ i) {
                    for (int j = 0; j < 3; ++ j) {
                        U += 0.5 * s[i][j] * eps[i][j];
                    }
                }
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    Eloc(aa[0], aa[1], aa[2]) += w * J * v.val * U;
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
        double tr = eps[0][0] + eps[1][1] + eps[2][2];

        for (int i = 0; i < 3; ++ i) {
            for (int j = 0; j < 3; ++ j) {
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
        double r = pow(x[0] - 1, 2) + pow(x[1] - 1, 2) + pow(x[2] - 1, 2);
        double a = - 1 * f * std::exp(- 10 * r);
        return {a, a, a};
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

    void update_u_and_v() {
        using std::swap;
        swap(now.vx, prev.vx);
        swap(now.vy, prev.vy);
        swap(now.vz, prev.vz);

        swap(now.ux, prev.ux);
        swap(now.uy, prev.uy);
        swap(now.uz, prev.uz);

        auto dt = steps.dt / 3;
        auto tt = 0.5 * dt * dt;

        // add(now.vx, now.ax, dt);
        // add(now.vy, now.ay, dt);
        // add(now.vz, now.az, dt);

        // add(now.ux, now.vx, dt);
        // add(now.uy, now.vy, dt);
        // add(now.uz, now.vz, dt);

        // add(now.ux, now.ax, -tt);
        // add(now.uy, now.ay, -tt);
        // add(now.uz, now.az, -tt);

        double beta = 0.25;
        double gamma = 0.5;

        add(now.ux, prev.vx, dt);
        add(now.uy, prev.vy, dt);
        add(now.uz, prev.vz, dt);
        add(now.ux, prev.ax, tt * (1 - 2 * beta));
        add(now.uy, prev.ay, tt * (1 - 2 * beta));
        add(now.uz, prev.az, tt * (1 - 2 * beta));
        add(now.ux, now.ax, tt * 2 * beta);
        add(now.uy, now.ay, tt * 2 * beta);
        add(now.uz, now.az, tt * 2 * beta);

        add(now.vx, prev.ax, dt * (1 - gamma));
        add(now.vy, prev.ay, dt * (1 - gamma));
        add(now.vz, prev.az, dt * (1 - gamma));

        add(now.vx, now.ax, dt * gamma);
        add(now.vy, now.ay, dt * gamma);
        add(now.vz, now.az, dt * gamma);

    }

    void step_split(int /*iter*/, double t) {
        using std::swap;

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_x(e, t);
            executor.synchronized([&] {
                apply_local_contribution(local, e);
            });
        });

        ads_solve(now.ax, buffer, ads::dim_data{Kx, x.ctx}, y.data(), z.data());
        ads_solve(now.ay, buffer, ads::dim_data{Kx, x.ctx}, y.data(), z.data());
        ads_solve(now.az, buffer, ads::dim_data{Kx, x.ctx}, y.data(), z.data());
        update_u_and_v();

        swap(now, prev);

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_y(e, t);
            executor.synchronized([&] {
                apply_local_contribution(local, e);
            });
        });

        ads_solve(now.ax, buffer, x.data(), ads::dim_data{Ky, y.ctx}, z.data());
        ads_solve(now.ay, buffer, x.data(), ads::dim_data{Ky, y.ctx}, z.data());
        ads_solve(now.az, buffer, x.data(), ads::dim_data{Ky, y.ctx}, z.data());
        update_u_and_v();

        swap(now, prev);

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_z(e, t);
            executor.synchronized([&] {
                apply_local_contribution(local, e);
            });
        });

        ads_solve(now.ax, buffer, x.data(), y.data(), ads::dim_data{Kz, z.ctx});
        ads_solve(now.ay, buffer, x.data(), y.data(), ads::dim_data{Kz, z.ctx});
        ads_solve(now.az, buffer, x.data(), y.data(), ads::dim_data{Kz, z.ctx});
        update_u_and_v();
    }

    void step(int /*iter*/, double t) override {
        using std::swap;

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_truncated(e, t);
            executor.synchronized([&] {
                apply_local_contribution(local, e);
            });
        });

        auto data_x = ads::dim_data{Kx, x.ctx};
        auto data_y = ads::dim_data{Ky, y.ctx};
        auto data_z = ads::dim_data{Kz, z.ctx};

        ads_solve(now.ax, buffer, data_x, data_y, data_z);
        ads_solve(now.ay, buffer, data_x, data_y, data_z);
        ads_solve(now.az, buffer, data_x, data_y, data_z);

        update_u_and_v();
    }

    void after_step(int iter, double t) override {
        iter += 1;
        if (iter % save_every == 0) {
            std::cout << "** Iteration " << iter << ", t = " << t + steps.dt << std::endl;

            double Ek = kinetic_energy();
            double Ep = potential_energy();
            compute_potential_energy();

            std::cout << "Kinetic energy: " << Ek << std::endl;
            std::cout << "Potential energy: " << Ep << std::endl;
            std::cout << "Total energy: " << Ek + Ep << std::endl;

            std::cout << "Total disp:   : " << total() << std::endl;
            std::cout << std::endl;

            // output.to_file("out_%d.vti", iter,
            //                output.evaluate(now.ux),
            //                output.evaluate(now.uy),
            //                output.evaluate(now.uz),
            //                output.evaluate(energy));
            // std::cout << iter << " " << t << " " << Ek << " " << Ep << " " << Ek + Ep << std::endl;

            // std::cout << "Step " << iter << ", t = " << t << std::endl;
            // int num = iter / save_every;
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
        return v(idx[0], idx[1], idx[2]);
    }
};
}




#endif /* PROBLEMS_ELASTICITY_ELASTICITY_VICTOR_HPP_ */
