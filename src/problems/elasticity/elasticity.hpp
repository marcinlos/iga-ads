#ifndef PROBLEMS_ELASTICITY_ELASTICITY_HPP_
#define PROBLEMS_ELASTICITY_ELASTICITY_HPP_

#include <cmath>
#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"


namespace problems {

    class linear_elasticity : public ads::simulation_3d {
        using Base = simulation_3d;
        using tensor = double[3][3];

        struct state {
            vector_type ux, uy, uz;
            vector_type vx, vy, vz;
            state(std::array<std::size_t, 3> shape)
            : ux{ shape }, uy{ shape }, uz{ shape }
            , vx{ shape }, vy{ shape }, vz{ shape }
            { }
        };

        state now, prev;
        vector_type energy;

        ads::output_manager<3> output;

        ads::galois_executor executor{8};

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
        linear_elasticity(const ads::config_3d& config)
        : Base{ config }
        , now{ shape() }, prev{ shape() }
        , energy{ shape() }
        , output{ x.B, y.B, z.B, 50 }
        { }

    private:

        void before() override {
            prepare_matrices();
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
                value_type ux = eval_fun(prev.ux, e, q);
                value_type uy = eval_fun(prev.uy, e, q);
                value_type uz = eval_fun(prev.uz, e, q);
                value_type vx = eval_fun(prev.vx, e, q);
                value_type vy = eval_fun(prev.vy, e, q);
                value_type vz = eval_fun(prev.vz, e, q);

                tensor eps = {
                    {         ux.dx,         0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx) },
                    { 0.5 * (ux.dy + uy.dx),         uy.dy,         0.5 * (uy.dz + uz.dy) },
                    { 0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy),         uz.dz         }
                };
                tensor s{};
                stress_tensor(s, eps);
                auto F = force(x, t);

                for (auto a : dofs_on_element(e)) {
                    value_type b = eval_basis(e, q, a);

                    double rho = 1;
                    double axb = (-s[0][0] * b.dx - s[0][1] * b.dy - s[0][2] * b.dz + F[0] * b.val) / rho;
                    double ayb = (-s[1][0] * b.dx - s[1][1] * b.dy - s[1][2] * b.dz + F[1] * b.val) / rho;
                    double azb = (-s[2][0] * b.dx - s[2][1] * b.dy - s[2][2] * b.dz + F[2] * b.val) / rho;

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
            for (auto e : elements()) {
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    double w = weigth(q);
                    value_type vx = eval_fun(now.vx, e, q);
                    value_type vy = eval_fun(now.vy, e, q);
                    value_type vz = eval_fun(now.vz, e, q);
                    E += w * J * 0.5 * (vx.val * vx.val + vy.val * vy.val + vz.val * vz.val);
                }
            }
            return E;
        }

        double potential_energy() const {
            double E = 0;
            for (auto e : elements()) {
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
                    E += w * J * U;
                }
            }
            return E;
        }

        void compute_potential_energy() {
            zero(energy);
            for (auto e : elements()) {
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
                        value_type v = eval_basis(e, q, a);
                        energy(a[0], a[1], a[2]) += w * J * v.val * U;
                    }
                }
            }
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
            double lambda = 1;
            double mi = 1;
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
            constexpr double t0 = 0.02;
            double tt = t / t0;
            double f = tt < 1 ? pow(tt * (1 - tt), 2) : 0;
            double r = pow(x[0] - 1, 2) + pow(x[1] - 1, 2) + pow(x[2] - 1, 2);
            double a = - 10 * f * std::exp(- 10 * r);
            return {a, a, a};
        }

        void before_step(int /*iter*/, double /*t*/) override {
            using std::swap;
            swap(now, prev);
        }

        void step(int /*iter*/, double t) override {
            compute_rhs(t);
            for_all(now, [this](vector_type& a) { solve(a); });
        }

        void after_step(int iter, double t) override {
            if (iter % 10 == 0) {
                std::cout << "** Iteration " << iter << ", t = " << t << std::endl;

                double Ek = kinetic_energy();
                double Ep = potential_energy();
                compute_potential_energy();

                std::cout << "Kinetic energy: " << Ek << std::endl;
                std::cout << "Potential energy: " << Ep << std::endl;
                std::cout << "Total energy: " << Ek + Ep << std::endl;

                std::cout << "Total disp:   : " << total() << std::endl;
                std::cout << std::endl;

                output.to_file("out_%d.vti", iter*10,
                               output.evaluate(now.ux),
                               output.evaluate(now.uy),
                               output.evaluate(now.uz),
                               output.evaluate(energy));
            }
        }


        double& ref(vector_type& v, index_type idx) const {
            return v(idx[0], idx[1], idx[2]);
        }
    };
}



#endif /* PROBLEMS_ELASTICITY_ELASTICITY_HPP_ */
