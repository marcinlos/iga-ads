// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ELASTICITY_IMPLICIT_HPP
#define ELASTICITY_IMPLICIT_HPP

#include <cmath>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"

namespace problems {

class implicit_elasticity : public ads::simulation_3d {
    using Base = simulation_3d;
    using tensor = double[3][3];

    struct state {
        vector_type ux, uy, uz;
        vector_type vx, vy, vz;
        explicit state(std::array<int, 3> shape)
        : ux{shape}
        , uy{shape}
        , uz{shape}
        , vx{shape}
        , vy{shape}
        , vz{shape} { }
    };

    state now, prev;
    vector_type energy;

    ads::output_manager<3> output;

    ads::galois_executor executor{8};

    ads::lin::band_matrix Dx, Dy, Dz;
    ads::lin::band_matrix Kx, Ky, Kz;

    static constexpr double lambda = 1;
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
    }

public:
    implicit_elasticity(const ads::config_3d& config, int save_every)
    : Base{config}
    , now{shape()}
    , prev{shape()}
    , energy{shape()}
    , output{x.B, y.B, z.B, 50}
    , Dx{x.p, x.p, x.B.dofs()}
    , Dy{y.p, y.p, y.B.dofs()}
    , Dz{z.p, z.p, z.B.dofs()}
    , Kx{x.p, x.p, x.B.dofs()}
    , Ky{y.p, y.p, y.B.dofs()}
    , Kz{z.p, z.p, z.B.dofs()}
    , save_every{save_every} {
        double hh = steps.dt * steps.dt * 1.0 / 3 * mi;
        matrix(Kx, x.basis, hh);
        matrix(Ky, y.basis, hh);
        matrix(Kz, z.basis, hh);

        double h = steps.dt * steps.dt * 1.0 / 3 * (lambda + 2 * mi);
        matrix(Dx, x.basis, h);
        matrix(Dy, y.basis, h);
        matrix(Dz, z.basis, h);
    }

private:
    void matrix(ads::lin::band_matrix& K, const ads::basis_data& d, double h) {
        for (ads::element_id e = 0; e < d.elements; ++e) {
            for (int q = 0; q < d.quad_order; ++q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++a) {
                    for (int b = 0; b + first <= last; ++b) {
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
        ads::lin::factorize(Kz, z.ctx);
        ads::lin::factorize(Dx, x.ctx);
        ads::lin::factorize(Dy, y.ctx);
        ads::lin::factorize(Dz, z.ctx);
    }

    void compute_rhs(double t) {
        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution(e, t);
            executor.synchronized([&] { apply_local_contribution(local, e); });
        });
    }

    state local_contribution(index_type e, double t) const {
        auto local = state{local_shape()};

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weight(q);
            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ux = eval_fun(prev.ux, e, q);
            value_type uy = eval_fun(prev.uy, e, q);
            value_type uz = eval_fun(prev.uz, e, q);

            tensor eps = {
                {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
            };
            tensor s{};
            stress_tensor(s, eps);
            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double rho = 1;
                double axb =
                    (-s[0][0] * b.dx - s[0][1] * b.dy - s[0][2] * b.dz + F[0] * b.val) / rho;
                double ayb =
                    (-s[1][0] * b.dx - s[1][1] * b.dy - s[1][2] * b.dz + F[1] * b.val) / rho;
                double azb =
                    (-s[2][0] * b.dx - s[2][1] * b.dy - s[2][2] * b.dz + F[2] * b.val) / rho;

                double dt = steps.dt;
                double t2 = dt * dt / 2;
                auto aa = dof_global_to_local(e, a);
                ref(local.ux, aa) += ((ux.val + dt * vx.val) * b.val + t2 * axb) * w * J;
                ref(local.uy, aa) += ((uy.val + dt * vy.val) * b.val + t2 * ayb) * w * J;
                ref(local.uz, aa) += ((uz.val + dt * vz.val) * b.val + t2 * azb) * w * J;

                ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                ref(local.vz, aa) += (vz.val * b.val + dt * azb) * w * J;
            }
        }
        return local;
    }

    state local_contribution_x(index_type e, double t) const {
        auto local = state{local_shape()};
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weight(q);
            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ux = eval_fun(prev.ux, e, q) + dt * vx;
            value_type uy = eval_fun(prev.uy, e, q) + dt * vy;
            value_type uz = eval_fun(prev.uz, e, q) + dt * vz;

            tensor eps = {
                {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
            };
            tensor s{};
            stress_tensor(s, eps);
            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double rho = 1;
                double axb =
                    (-s[0][0] * b.dx - s[0][1] * b.dy - s[0][2] * b.dz + F[0] * b.val) / rho;
                double ayb =
                    (-s[1][0] * b.dx - s[1][1] * b.dy - s[1][2] * b.dz + F[1] * b.val) / rho;
                double azb =
                    (-s[2][0] * b.dx - s[2][1] * b.dy - s[2][2] * b.dz + F[2] * b.val) / rho;

                axb += (lambda + 2 * mi) * dt * vx.dx * b.dx;
                ayb += mi * dt * vy.dx * b.dx;
                azb += mi * dt * vz.dx * b.dx;

                auto aa = dof_global_to_local(e, a);
                ref(local.ux, aa) += ux.val * b.val * w * J;
                ref(local.uy, aa) += uy.val * b.val * w * J;
                ref(local.uz, aa) += uz.val * b.val * w * J;

                ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                ref(local.vz, aa) += (vz.val * b.val + dt * azb) * w * J;
            }
        }
        return local;
    }
    state local_contribution_y(index_type e, double t) const {
        auto local = state{local_shape()};
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weight(q);
            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ux = eval_fun(prev.ux, e, q) + dt * vx;
            value_type uy = eval_fun(prev.uy, e, q) + dt * vy;
            value_type uz = eval_fun(prev.uz, e, q) + dt * vz;

            tensor eps = {
                {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
            };
            tensor s{};
            stress_tensor(s, eps);
            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double rho = 1;
                double axb =
                    (-s[0][0] * b.dx - s[0][1] * b.dy - s[0][2] * b.dz + F[0] * b.val) / rho;
                double ayb =
                    (-s[1][0] * b.dx - s[1][1] * b.dy - s[1][2] * b.dz + F[1] * b.val) / rho;
                double azb =
                    (-s[2][0] * b.dx - s[2][1] * b.dy - s[2][2] * b.dz + F[2] * b.val) / rho;

                axb += mi * dt * vx.dx * b.dx;
                ayb += (lambda + 2 * mi) * dt * vy.dx * b.dx;
                azb += mi * dt * vz.dx * b.dx;

                auto aa = dof_global_to_local(e, a);
                ref(local.ux, aa) += ux.val * b.val * w * J;
                ref(local.uy, aa) += uy.val * b.val * w * J;
                ref(local.uz, aa) += uz.val * b.val * w * J;

                ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                ref(local.vz, aa) += (vz.val * b.val + dt * azb) * w * J;
            }
        }
        return local;
    }

    state local_contribution_z(index_type e, double t) const {
        auto local = state{local_shape()};
        double dt = steps.dt / 3.0;

        double J = jacobian(e);
        for (auto q : quad_points()) {
            auto x = point(e, q);
            double w = weight(q);
            value_type vx = eval_fun(prev.vx, e, q);
            value_type vy = eval_fun(prev.vy, e, q);
            value_type vz = eval_fun(prev.vz, e, q);

            value_type ux = eval_fun(prev.ux, e, q) + dt * vx;
            value_type uy = eval_fun(prev.uy, e, q) + dt * vy;
            value_type uz = eval_fun(prev.uz, e, q) + dt * vz;

            tensor eps = {
                {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
            };
            tensor s{};
            stress_tensor(s, eps);
            auto F = force(x, t);

            for (auto a : dofs_on_element(e)) {
                value_type b = eval_basis(e, q, a);

                double rho = 1;
                double axb =
                    (-s[0][0] * b.dx - s[0][1] * b.dy - s[0][2] * b.dz + F[0] * b.val) / rho;
                double ayb =
                    (-s[1][0] * b.dx - s[1][1] * b.dy - s[1][2] * b.dz + F[1] * b.val) / rho;
                double azb =
                    (-s[2][0] * b.dx - s[2][1] * b.dy - s[2][2] * b.dz + F[2] * b.val) / rho;

                axb += mi * dt * vx.dx * b.dx;
                ayb += mi * dt * vy.dx * b.dx;
                azb += (lambda + 2 * mi) * dt * vz.dx * b.dx;

                auto aa = dof_global_to_local(e, a);
                ref(local.ux, aa) += ux.val * b.val * w * J;
                ref(local.uy, aa) += uy.val * b.val * w * J;
                ref(local.uz, aa) += uz.val * b.val * w * J;

                ref(local.vx, aa) += (vx.val * b.val + dt * axb) * w * J;
                ref(local.vy, aa) += (vy.val * b.val + dt * ayb) * w * J;
                ref(local.vz, aa) += (vz.val * b.val + dt * azb) * w * J;
            }
        }
        return local;
    }

    void apply_local_contribution(const state& loc, index_type e) {
        update_global_rhs(now.ux, loc.ux, e);
        update_global_rhs(now.uy, loc.uy, e);
        update_global_rhs(now.uz, loc.uz, e);
        update_global_rhs(now.vx, loc.vx, e);
        update_global_rhs(now.vy, loc.vy, e);
        update_global_rhs(now.vz, loc.vz, e);
    }

    double kinetic_energy() const {
        double E = 0;
        executor.for_each(elements(), [&](index_type e) {
            double Eloc = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type vx = eval_fun(now.vx, e, q);
                value_type vy = eval_fun(now.vy, e, q);
                value_type vz = eval_fun(now.vz, e, q);
                Eloc += w * J * 0.5 * (vx.val * vx.val + vy.val * vy.val + vz.val * vz.val);
            }
            executor.synchronized([&] { E += Eloc; });
        });
        return E;
    }

    double potential_energy() const {
        double E = 0;
        executor.for_each(elements(), [&](index_type e) {
            double Eloc = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type ux = eval_fun(now.ux, e, q);
                value_type uy = eval_fun(now.uy, e, q);
                value_type uz = eval_fun(now.uz, e, q);

                tensor eps = {
                    {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                    {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                    {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
                };
                tensor s{};
                stress_tensor(s, eps);
                double U = 0;
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        U += 0.5 * s[i][j] * eps[i][j];
                    }
                }
                Eloc += w * J * U;
            }
            executor.synchronized([&] { E += Eloc; });
        });
        return E;
    }

    void compute_potential_energy() {
        zero(energy);
        executor.for_each(elements(), [&](index_type e) {
            auto Eloc = element_rhs();
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type ux = eval_fun(now.ux, e, q);
                value_type uy = eval_fun(now.uy, e, q);
                value_type uz = eval_fun(now.uz, e, q);

                tensor eps = {
                    {ux.dx, 0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx)},
                    {0.5 * (ux.dy + uy.dx), uy.dy, 0.5 * (uy.dz + uz.dy)},
                    {0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy), uz.dz},
                };
                tensor s{};
                stress_tensor(s, eps);
                double U = 0;
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        U += 0.5 * s[i][j] * eps[i][j];
                    }
                }
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    Eloc(aa[0], aa[1], aa[2]) += w * J * v.val * U;
                }
            }
            executor.synchronized([&]() { update_global_rhs(energy, Eloc, e); });
        });
        solve(energy);
    }

    double total() const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type vx = eval_fun(now.vx, e, q);
                E += w * J * vx.val;
            }
        }
        return E;
    }

    void stress_tensor(tensor& s, const tensor& eps) const {
        double tr = eps[0][0] + eps[1][1] + eps[2][2];

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
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
        double a = -10 * f * std::exp(-10 * r);
        return {a, a, a};
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

    void step2(int /*iter*/, double t) {
        compute_rhs(t);
        for_all(now, [this](vector_type& a) { solve(a); });
    }

    void step(int /*iter*/, double t) override {
        using std::swap;

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_x(e, t);
            executor.synchronized([&] { apply_local_contribution(local, e); });
        });

        ads_solve(now.ux, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uy, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uz, buffer, x.data(), y.data(), z.data());

        ads_solve(now.vx, buffer, ads::dim_data{Dx, x.ctx}, y.data(), z.data());
        ads_solve(now.vy, buffer, ads::dim_data{Kx, x.ctx}, y.data(), z.data());
        ads_solve(now.vz, buffer, ads::dim_data{Kx, x.ctx}, y.data(), z.data());

        swap(now, prev);

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_y(e, t);
            executor.synchronized([&] { apply_local_contribution(local, e); });
        });

        ads_solve(now.ux, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uy, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uz, buffer, x.data(), y.data(), z.data());

        ads_solve(now.vx, buffer, x.data(), ads::dim_data{Ky, y.ctx}, z.data());
        ads_solve(now.vy, buffer, x.data(), ads::dim_data{Dy, y.ctx}, z.data());
        ads_solve(now.vz, buffer, x.data(), ads::dim_data{Ky, y.ctx}, z.data());

        swap(now, prev);

        for_all(now, [](vector_type& a) { zero(a); });
        executor.for_each(elements(), [&](index_type e) {
            auto local = local_contribution_z(e, t);
            executor.synchronized([&] { apply_local_contribution(local, e); });
        });

        ads_solve(now.ux, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uy, buffer, x.data(), y.data(), z.data());
        ads_solve(now.uz, buffer, x.data(), y.data(), z.data());

        ads_solve(now.vx, buffer, x.data(), y.data(), ads::dim_data{Kz, z.ctx});
        ads_solve(now.vy, buffer, x.data(), y.data(), ads::dim_data{Kz, z.ctx});
        ads_solve(now.vz, buffer, x.data(), y.data(), ads::dim_data{Dz, z.ctx});
    }

    void after_step(int iter, double t) override {
        iter += 1;
        if (iter % save_every == 0) {
            std::cout << "** Iteration " << iter << ", t = " << t + steps.dt << std::endl;

            double Ek = kinetic_energy();
            double Ep = potential_energy();
            // compute_potential_energy();

            std::cout << "Kinetic energy: " << Ek << std::endl;
            std::cout << "Potential energy: " << Ep << std::endl;
            std::cout << "Total energy: " << Ek + Ep << std::endl;

            std::cout << "Total disp:   : " << total() << std::endl;
            std::cout << std::endl;

            // output.to_file("out_%d.vti", iter*10,
            // output.evaluate(now.ux),
            // output.evaluate(now.uy),
            // output.evaluate(now.uz),
            // output.evaluate(energy));
            // std::cout << iter << " " << t << " " << Ek << " " << Ep << " " << Ek + Ep <<
            // std::endl; std::cout << "Step " << (iter + 1) << ", t = " << t << std::endl;

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

    double& ref(vector_type& v, index_type idx) const { return v(idx[0], idx[1], idx[2]); }
};
}  // namespace problems

#endif  // ELASTICITY_IMPLICIT_HPP
