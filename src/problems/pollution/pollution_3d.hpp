#ifndef PROBLEMS_POLLLUTION_3D_HPP_
#define PROBLEMS_POLLLUTION_3D_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"


namespace ads {

class pollution_3d : public simulation_3d {
private:
    using Base = simulation_3d;
    vector_type u, u_prev;

    output_manager<3> output;
    galois_executor executor{8};

    lin::band_matrix Kx, Ky, Kz;
    lin::solver_ctx Kx_ctx, Ky_ctx, Kz_ctx;

    int save_every = 1;

    double ambient = 1e-6; // g/m^3

    double Vd = 0.1;
    point_type c_diff{{ 50, 50, 0.5 }}; // m/s^2

    double wind_angle = 2 * M_PI / 4;
    double wind_speed = 5; // m/s

    point_type wind{{ wind_speed * cos(wind_angle), wind_speed * sin(wind_angle), 0.0 }};

    double absorbed = 0.0;


public:
    pollution_3d(const config_3d& config)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , output{ x.B, y.B, z.B, 100 }
    , Kx{x.p, x.p, x.B.dofs()}
    , Ky{y.p, y.p, y.B.dofs()}
    , Kz{z.p, z.p, z.B.dofs()}
    , Kx_ctx{ Kx }
    , Ky_ctx{ Ky }
    , Kz_ctx{ Kz }
    {
        matrix(Kx, x.basis, steps.dt / 3, c_diff[0], wind[0]);
        matrix(Ky, y.basis, steps.dt / 3, c_diff[1], wind[1]);
        matrix(Kz, z.basis, steps.dt / 3, c_diff[2], wind[2]);
    }

private:

    void diffusion_matrix(lin::band_matrix& M, const basis_data& d, double h, double diffusion) {
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto da = d.b[e][q][1][a];
                        auto db = d.b[e][q][1][b];
                        auto diff = diffusion * h * da * db;
                        M(ia, ib) += diff * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void advection_matrix(lin::band_matrix& M, const basis_data& d, double h, double advection) {
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto db = d.b[e][q][1][b];
                        auto adv = advection * h * va * db;

                        M(ia, ib) += adv * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void advection_matrix_conv_by_parts(lin::band_matrix& M, const basis_data& d, double h, double advection) {
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto da = d.b[e][q][1][a];
                        auto vb = d.b[e][q][0][b];
                        auto adv = - advection * h * da * vb;

                        M(ia, ib) += adv * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void matrix(lin::band_matrix& K, const basis_data& d, double h, double diffusion, double advection) {
        gram_matrix_1d(K, d);
        diffusion_matrix(K, d, h, diffusion);
        advection_matrix(K, d, h, advection);
    }

    double emission(double x, double y, double z) {
        double dx = (x - 3000) / 25;
        double dy = (y - 2000) / 25;
        double dz = (z - 2000) / 25;

        double r2 = std::min((dx * dx + dy * dy + dz * dz), 1.0);
        // return 0;
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1); // g/m^3
    }

    void impose_bc(vector_type& v) const {
        // {
        //     // {0} x X
        //     lin::vector buf{{ y.dofs() }};
        //     compute_projection(buf, y.basis, [&](double t) {
        //         return std::max(ambient, (200 - std::abs(t - 1500)) / 200.0);
        //     });
        //     for (int i = 0; i < y.dofs(); ++ i) {
        //         v(0, i) = buf(i);
        //     }
        // }
        // {
        //     // X x {0}
        //     lin::vector buf{{ x.dofs() }};
        //     compute_projection(buf, x.basis, [&](double t) {
        //         return std::max(ambient, (100 - std::abs(t - 1000)) / 100.0);
        //     });
        //     for (int i = 0; i < x.dofs(); ++ i) {
        //         v(i, 0) = buf(i);
        //     }
        // }
        // v(0, 0) = ambient;
    }

    void fix_dof(int k, const dimension& dim, lin::band_matrix& K) {
        int last = dim.dofs() - 1;
        for (int i = clamp(k - dim.p, 0, last); i <= clamp(k + dim.p, 0, last); ++ i) {
            Kx(k, i) = 0;
        }
        Kx(k, k) = 1;
    }

    void prepare_implicit_matrices() {
        Kx.zero();
        Ky.zero();
        Kz.zero();
        auto h = steps.dt / 3;
        matrix(Kx, x.basis, h, c_diff[0], wind[0]);
        matrix(Ky, y.basis, h, c_diff[1], wind[1]);
        matrix(Kz, z.basis, h, c_diff[2], wind[2]);

        // fix_dof(0, x, Kx);
        // fix_dof(0, y, Ky);

        lin::factorize(Kx, Kx_ctx);
        lin::factorize(Ky, Ky_ctx);
        lin::factorize(Kz, Kz_ctx);
    }

    void prepare_matrices() {
        // x.fix_left();
        // y.fix_left();
        Base::prepare_matrices();

        prepare_implicit_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double, double, double) { return ambient; };

        std::cout << "Projecting initial state" << std::endl;
        projection(u, init);
        impose_bc(u);
        solve(u);

        std::cout << "Outputting initial state" << std::endl;
        output.to_file(u, "out_0.vti");
        std::cout << "Done, ready for computation" << std::endl;
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        using std::swap;

        compute_rhs_x();
        impose_bc(u);
        ads_solve(u, buffer, dim_data{Kx, Kx_ctx}, y.data(), z.data());

        swap(u, u_prev);

        compute_rhs_y();
        impose_bc(u);
        ads_solve(u, buffer, x.data(), dim_data{Ky, Ky_ctx}, z.data());

        swap(u, u_prev);

        compute_rhs_z();
        impose_bc(u);
        ads_solve(u, buffer, x.data(), y.data(), dim_data{Kz, Kz_ctx});
    }

    void after_step(int iter, double t) override {
        if ((iter + 1) % save_every == 0) {
            output.to_file(u, "out_%d.vti", (iter + 1) / save_every);
            analyze(iter, t + 0.5 * steps.dt);
        }

        auto s = t / 150;
        auto phase = sin(s) + 0.5 * sin(2.3*s);
        wind_angle = M_PI / 3 * phase + 1.5 * M_PI / 4;
        double vertical = wind_speed / 3 * sin(s);

        wind = { wind_speed * cos(wind_angle), wind_speed * sin(wind_angle), vertical };
        prepare_implicit_matrices();

    }

    point_type local_wind(double x, double y) const {
        // auto dx = x - 2500;
        // auto dy = y - 2500;
        // auto r = std::sqrt(dx * dx + dy * dy);
        // auto a = std::exp(-r / 50) * 1 / (r + 1);
        // return {- a * y, a * x};
        // return {0, sin(4 * M_PI * x / 1000) };
        return {};
    }

    // void boundary_integral(index_type e, vector_type& U, double h) {
    //     if (e[0] == x.B.elements() - 1) {
    //         double absorbed_loc = 0.0;
    //         double J = y.basis.J[e[1]];
    //         for (int q = 0; q < y.basis.quad_order; ++ q) {
    //             double w = y.basis.w[q];
    //             int first = y.basis.first_dof(e[1]);
    //             int last = y.basis.last_dof(e[1]);

    //             double u = 0.0;
    //             for (int a = 0; a + first <= last; ++ a) {
    //                 auto B = y.basis.b[e[1]][q][0][a];
    //                 double c = u_prev(x.basis.last_dof(e[0]), a + first);
    //                 u += c * B;
    //             }

    //             for (int a = 0; a + first <= last; ++ a) {
    //                 auto B = y.basis.b[e[1]][q][0][a];
    //                 auto aa = dof_global_to_local(e, {x.basis.last_dof(e[0]), a + first});
    //                 auto prod = - Vd * u * B;

    //                 U(aa[0], aa[1]) += h * prod * w * J;
    //             }
    //             absorbed_loc += h * Vd * u * w * J;;
    //         }
    //         executor.synchronized([&]{ absorbed += absorbed_loc; });
    //     }
    // }

    void compute_rhs_x() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto h = steps.dt / 3;

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval_fun(u_prev, e, q);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double conv_term = (wind[1] * u.dy + wind[2] * u.dz) * v.val;
                    double gradient_prod = c_diff[1] * u.dy * v.dy + c_diff[2] * u.dz * v.dz;
                    double e = emission(x[0], x[1], x[2]) * v.val;
                    double val = u.val * v.val + h * (- conv_term - gradient_prod + e);

                    U(aa[0], aa[1], aa[2]) += val * w * J;
                }
            }
            // boundary_integral(e, U, h);

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    void compute_rhs_y() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto h = steps.dt / 3;

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval_fun(u_prev, e, q);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double conv_term = (wind[0] * u.dx + wind[2] * u.dz) * v.val;
                    double gradient_prod = c_diff[0] * u.dx * v.dx + c_diff[2] * u.dz * v.dz;
                    double e = emission(x[0], x[1], x[2]) * v.val;
                    double val = u.val * v.val + h * (- conv_term - gradient_prod + e);
                    U(aa[0], aa[1], aa[2]) += val * w * J;
                }
            }
            // boundary_integral(e, U, h);

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    void compute_rhs_z() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto h = steps.dt / 3;

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval_fun(u_prev, e, q);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double conv_term = (wind[0] * u.dx + wind[1] * u.dy) * v.val;
                    double gradient_prod = wind[0] * u.dx * v.dx + wind[1] * u.dy * v.dy;
                    double e = emission(x[0], x[1], x[2]) * v.val;
                    double val = u.val * v.val + h * (- conv_term - gradient_prod + e);
                    U(aa[0], aa[1], aa[2]) += val * w * J;
                }
            }
            // boundary_integral(e, U, h);

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    void analyze(int iter, double t) {
        double total = 0.0;
        double emission_rate = 0.0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval_fun(u_prev, e, q);

                total += u.val * w * J;
                emission_rate += emission(x[0], x[1], x[2]) * w * J;
            }
        }

        using namespace std;

        total /= 1000; // g -> kg
        emission_rate /= 1000; // g -> kg

        auto emitted = t * emission_rate; // kg
        auto volume = 5000.0 * 5000.0 * 5000.0; // m^3
        auto initial = volume * ambient / 1000; // kg
        auto absorbed_kg = absorbed / 1000;
        auto loss = initial + emitted - total - absorbed_kg;
        auto loss_percent = 100 * loss / emitted;

        std::cout << "Step " << (iter + 1) << ":"
                  << "  total " << setw(8) << setprecision(5) << total << " kg "
                  << "  absorbed " << setw(8) << setprecision(5) << absorbed_kg << " kg "
                  << "  (loss "
                  << setw(6) << loss << " kg,  "
                  << setw(5) << loss_percent << "%"
                  << ")" << std::endl;
    }
};

}

#endif /* PROBLEMS_POLLLUTION_3D_HPP_*/
