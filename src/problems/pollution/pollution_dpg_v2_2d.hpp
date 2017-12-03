#ifndef PROBLEMS_POLLUTION_DPG_V2_2D_HPP_
#define PROBLEMS_POLLUTION_DPG_V2_2D_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class pollution_dpg_v2_2d : public simulation_2d {
private:
    using Base = simulation_2d;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;
    lin::dense_matrix Kx_x, Kx_y, Ky_x, Ky_y;

    lin::dense_matrix Kx_x_nf, Kx_y_nf, Ky_x_nf, Ky_y_nf; // non-factorized copies

    lin::solver_ctx Kxx_ctx, Kxy_ctx, Kyx_ctx, Kyy_ctx;

    lin::band_matrix Ax, Ay;
    lin::solver_ctx Ax_ctx, Ay_ctx;
    lin::band_matrix MUx, MUy;

    lin::band_matrix MUVx, MUVy;
    // lin::band_matrix Bx, By;
    lin::dense_matrix Bx, By;
    lin::dense_matrix Tx, Ty;

    vector_type u, u_prev;
    vector_type u_buffer;
    vector_type rhs1, rhs2;

    int save_every = 10;

    double ambient = 1e-6; // g/m^3
    // double Vd = 0.1;
    double alpha = 1e-6; // 0.1
    point_type c_diff{{ alpha, alpha }}; // m/s^2

    double wind_angle = M_PI / 6;
    double wind_speed = 100; // m/s

    point_type wind{{ wind_speed * cos(wind_angle), wind_speed * sin(wind_angle) }};

    double absorbed = 0.0;

    output_manager<2> output;

public:
    pollution_dpg_v2_2d(const config_2d& config, int k)
    // : Base{ higher_order(config, k) }
    // : Base{ repeat_nodes(config, k) }
    // : Base{ increase_elements(config, k) }
    : Base{ higher_order(repeat_nodes(config, k), k) }
    , Ux{ config.x, config.derivatives }
    , Uy{ config.y, config.derivatives }
    , Vx{ x }
    , Vy{ y }
    , Kx_x{Ux.dofs(), Ux.dofs()}
    , Kx_y{Uy.dofs(), Uy.dofs()}
    , Ky_x{Ux.dofs(), Ux.dofs()}
    , Ky_y{Uy.dofs(), Uy.dofs()}
    , Kx_x_nf{Ux.dofs(), Ux.dofs()}
    , Kx_y_nf{Uy.dofs(), Uy.dofs()}
    , Ky_x_nf{Ux.dofs(), Ux.dofs()}
    , Ky_y_nf{Uy.dofs(), Uy.dofs()}
    , Kxx_ctx{ Kx_x }
    , Kxy_ctx{ Kx_y }
    , Kyx_ctx{ Ky_x }
    , Kyy_ctx{ Ky_y }
    , Ax{Vx.p, Vx.p, Vx.dofs()}
    , Ay{y.p, Vy.p, Vy.dofs()}
    , Ax_ctx{ Ax }
    , Ay_ctx{ Ay }
    , MUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , MUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , MUVx{ Vx.p, Ux.p, Vx.dofs(), Ux.dofs() }
    , MUVy{ Vy.p, Uy.p, Vy.dofs(), Uy.dofs() }
    , Bx{ Vx.dofs(), Ux.dofs() }
    , By{ Vy.dofs(), Uy.dofs() }
    , Tx{ Vx.dofs(), Ux.dofs() }
    , Ty{ Vy.dofs(), Uy.dofs() }
    , u{{ Ux.dofs(), Uy.dofs() }}
    , u_prev{{ Ux.dofs(), Uy.dofs() }}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
    , rhs1{{ Vx.dofs(), Uy.dofs() }}, rhs2{{ Ux.dofs(), Vy.dofs() }}
    , output{ Ux.B, Uy.B, 400 }
    {
    }

private:

    static config_2d increase_elements(config_2d cfg, int k) {
        cfg.x = increase_elements(cfg.x, k);
        cfg.y = increase_elements(cfg.y, k);
        return cfg;
    }

    static dim_config increase_elements(dim_config cfg, int k) {
        cfg.elements = k;
        return cfg;
    }

    static config_2d repeat_nodes(config_2d cfg, int k) {
        cfg.x = repeat_nodes(cfg.x, k);
        cfg.y = repeat_nodes(cfg.y, k);
        return cfg;
    }

    static dim_config repeat_nodes(dim_config cfg, int k) {
        cfg.repeated_nodes += k;
        return cfg;
    }

    static dim_config higher_order(dim_config cfg, int k) {
        cfg.p += k;
        return cfg;
    }

    static config_2d higher_order(config_2d cfg, int k) {
        cfg.x = higher_order(cfg.x, k);
        cfg.y = higher_order(cfg.y, k);
        return cfg;
    }

    void prod_V(lin::band_matrix& M, const basis_data& bV) const {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                auto first = bV.first_dof(e);
                auto last  = bV.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = bV.b[e][q][0][a];
                        auto vb = bV.b[e][q][0][b];
                        auto da = bV.b[e][q][1][a];
                        auto db = bV.b[e][q][1][b];
                        auto m = va * vb + da * db;
                        M(ia, ib) += m * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void mass_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto va = bV.b[e][q][0][a];
                        auto vb = bU.b[e][q][0][b];
                        auto diff = va * vb;
                        M(ia, ib) += diff * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void diffusion_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV,
                          double h, double diffusion) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto da = bV.b[e][q][1][a];
                        auto db = bU.b[e][q][1][b];
                        auto diff = diffusion * h * da * db;
                        M(ia, ib) += diff * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void advection_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV,
                          double h, double advection) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto va = bV.b[e][q][0][a];
                        auto db = bU.b[e][q][1][b];
                        auto diff = advection * h * va * db;
                        M(ia, ib) += diff * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void fix_dof(int k, const dimension& dim, lin::dense_matrix& K) {
        int last = dim.dofs() - 1;
        for (int i = clamp(k - dim.p, 0, last); i <= clamp(k + dim.p, 0, last); ++ i) {
            K(k, i) = 0;
        }
        K(k, k) = 1;
    }

    void matrix(lin::dense_matrix& B, const basis_data& bU, const basis_data& bV,
                double h, double diffusion, double advection) {
        mass_matrix(B, bU, bV);
        diffusion_matrix(B, bU, bV, h, diffusion);
        advection_matrix(B, bU, bV, h, advection);
    }

    double emission(double x, double y, double u) {
        double dx = (x - 3000) / 25;
        double dy = (y - 400) / 25;
        double r2 = std::min((dx * dx + dy * dy), 1.0);
        // return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1); // g/m^3
        return 0;
    };

    void prepare_implicit_matrices() {
        // MUVx.zero();
        // MUVy.zero();
        Bx.zero();
        By.zero();
        Ax.zero();
        Ay.zero();
        MUx.zero();
        MUy.zero();
        Kx_x.zero();
        Kx_y.zero();
        Ky_x.zero();
        Ky_y.zero();

        // mass_matrix(MUVx, Ux.basis, x.basis);
        // mass_matrix(MUVy, Uy.basis, y.basis);

        matrix(Bx, Ux.basis, x.basis, steps.dt / 2, c_diff[0], wind[0]);
        matrix(By, Uy.basis, y.basis, steps.dt / 2, c_diff[1], wind[1]);

        prod_V(Ax, x.basis);
        prod_V(Ay, y.basis);
        gram_matrix_1d(MUx, Ux.basis);
        gram_matrix_1d(MUy, Uy.basis);

        lin::factorize(Ax, Ax_ctx);
        lin::factorize(Ay, Ay_ctx);

        // TODO: fill
        // Kx_x = Bx' Ax^-1 Bx
        // to_dense(Bx, Tx);
        Tx = Bx;
        solve_with_factorized(Ax, Tx, Ax_ctx);
        multiply(Bx, Tx, Kx_x, "T");

        // Kx_y = MUVy' MVy^-1 MUVy
        // to_dense(MUy, Ty);
        // solve_with_factorized(Uy.M, Ty, Uy.ctx);
        // multiply(MUy, Ty, Kx_y, "T");
        to_dense(MUy, Kx_y);

        // Ky_x = MUVx' MVx^-1 MUVx
        // to_dense(MUx, Tx);
        // solve_with_factorized(Ux.M, Tx, Ux.ctx);
        // multiply(MUx, Tx, Ky_x, "T");
        to_dense(MUx, Ky_x);

        // Ky_y = By' Ay^-1 By
        // to_dense(By, Ty);
        Ty = By;
        solve_with_factorized(Ay, Ty, Ay_ctx);
        multiply(By, Ty, Ky_y, "T");

        // lin::factorize(MUx, MUx_ctx);
        // lin::factorize(MUy, MUy_ctx);

        fix_dof(0, Uy, Kx_y);
        fix_dof(0, Uy, Ky_y);
        fix_dof(Uy.dofs() - 1, Uy, Kx_y);
        fix_dof(Uy.dofs() - 1, Uy, Ky_y);


        fix_dof(0, Ux, Ky_x);
        fix_dof(0, Ux, Kx_x);
        fix_dof(Ux.dofs() - 1, Ux, Ky_x);
        fix_dof(Ux.dofs() - 1, Ux, Kx_x);


        Kx_x_nf = Kx_x;
        Kx_y_nf = Kx_y;
        Ky_x_nf = Ky_x;
        Ky_y_nf = Ky_y;

        lin::factorize(Kx_x, Kxx_ctx);
        lin::factorize(Kx_y, Kxy_ctx);
        lin::factorize(Ky_x, Kyx_ctx);
        lin::factorize(Ky_y, Kyy_ctx);

    }

    void prepare_matrices() {
        Ux.factorize_matrix();
        Uy.factorize_matrix();
        // Base::prepare_matrices();
        prepare_implicit_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return ambient; };

        // projection(u, init);
        // solve(u);
        zero(u);

        output.to_file(u, "out_0.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs_x();
        ads_solve(rhs1, buffer, dim_data{Ax, Ax_ctx}, Uy.data());

        vector_type rhsx1   {{ Ux.dofs(), Uy.dofs() }};
        vector_type rhsx1_t {{ Uy.dofs(), Ux.dofs() }};

        vector_type u_t     {{ Uy.dofs(), Ux.dofs() }};

        // u = (Bx * MUVy)' rhs
        // =>
        // u = (Bx' (MUVy' rhs)')'
        // =>
        // rhsx = MUVy' rhs
        // u = (Bx' rhsx')'

        multiply(Bx, rhs1, rhsx1, Uy.dofs(), "T");
        lin::cyclic_transpose(rhsx1, rhsx1_t);
        multiply(MUy, rhsx1_t, u_t, Ux.dofs(), "T");
        lin::cyclic_transpose(u_t, u);



        lin::band_matrix MUx_loc{ Ux.p, Ux.p, Ux.dofs() };
        gram_matrix_1d(MUx_loc, Ux.basis);
        lin::solver_ctx ctx_x{ MUx_loc };
        lin::factorize(MUx_loc, ctx_x);

        lin::band_matrix MUy_loc{ Uy.p, Uy.p, Uy.dofs() };
        gram_matrix_1d(MUy_loc, Uy.basis);
        lin::solver_ctx ctx_y{ MUy_loc };
        lin::factorize(MUy_loc, ctx_y);


        lin::vector buf_y0{{ Ux.dofs() }};
        compute_projection(buf_y0, Ux.basis, [&](double t) {
            // return std::sin(t * M_PI);
            // return 1 - t;
            return 1;
        });
        lin::solve_with_factorized(MUx_loc, buf_y0, ctx_x);

        lin::vector buf_y1{{ Ux.dofs() }};
        compute_projection(buf_y1, Ux.basis, [&](double t) {
            // auto a =  std::sin(t * 3 * M_PI);
            // return a * a;
            return 0;
        });
        lin::solve_with_factorized(MUx_loc, buf_y1, ctx_x);

        lin::vector buf_x0{{ Uy.dofs() }};
        compute_projection(buf_x0, Uy.basis, [&](double t) {
            // auto a =  std::sin(t * 4 * M_PI);
            // return a * a;
            return t < 0.5 ? 1 : 0;
            // return t < 0.5 ? 1 - 2*t : 0;
        });
        lin::solve_with_factorized(MUy_loc, buf_x0, ctx_y);

        lin::vector buf_x1{{ Uy.dofs() }};
        compute_projection(buf_x1, Uy.basis, [&](double t) {
            // auto a =  std::sin(t * 2 * M_PI);
            // return a * a;
            return 0;
        });
        lin::solve_with_factorized(MUy_loc, buf_x1, ctx_y);

        for (int i = 0; i < Ux.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Ux.dofs(); ++ j) {
                val += buf_y0(j) * Kx_x_nf(i, j);
            }
            u(i, 0) = val;
        }

        for (int i = 0; i < Ux.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Ux.dofs(); ++ j) {
                val += buf_y1(j) * Kx_x_nf(i, j);
            }
            u(i, Uy.dofs() - 1) = val;
        }

        for (int i = 0; i < Uy.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Uy.dofs(); ++ j) {
                val += buf_x0(j) * Kx_y_nf(i, j);
            }
            u(0, i) = val;
        }

        for (int i = 0; i < Uy.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Uy.dofs(); ++ j) {
                val += buf_x1(j) * Kx_y_nf(i, j);
            }
            u(Ux.dofs() - 1, i) = val;
        }


        // ads_solve(u, u_buffer, dim_data{Kx_x, Kxx_ctx}, dim_data{Kx_y, Kxy_ctx});
        lin::solve_with_factorized(Kx_x, u, Kxx_ctx);
        auto F = lin::cyclic_transpose(u, u_buffer.data());
        lin::solve_with_factorized(Kx_y, F, Kxy_ctx);
        lin::cyclic_transpose(F, u);

        using std::swap;
        swap(u, u_prev);

        compute_rhs_y();
        ads_solve(rhs2, buffer, Ux.data(), dim_data{Ay, Ay_ctx});

        vector_type rhsx2   {{ Ux.dofs(), Vy.dofs() }};
        vector_type rhsx2_t {{ Vy.dofs(), Ux.dofs() }};
        // u = (MUVx * By)' rhs

        multiply(MUx, rhs2, rhsx2, Vy.dofs(), "T");
        lin::cyclic_transpose(rhsx2, rhsx2_t);
        multiply(By, rhsx2_t, u_t, Ux.dofs(), "T");
        lin::cyclic_transpose(u_t, u);




        for (int i = 0; i < Ux.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Ux.dofs(); ++ j) {
                val += buf_y0(j) * Ky_x_nf(i, j);
            }
            u(i, 0) = val;
        }

        for (int i = 0; i < Ux.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Ux.dofs(); ++ j) {
                val += buf_y1(j) * Ky_x_nf(i, j);
            }
            u(i, Uy.dofs() - 1) = val;
        }

        for (int i = 0; i < Uy.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Uy.dofs(); ++ j) {
                val += buf_x0(j) * Ky_y_nf(i, j);
            }
            u(0, i) = val;
        }

        for (int i = 0; i < Uy.dofs(); ++ i) {
            double val = 0;
            for (int j = 0; j < Uy.dofs(); ++ j) {
                val += buf_x1(j) * Ky_y_nf(i, j);
            }
            u(Ux.dofs() - 1, i) = val;
        }


        // ads_solve(u, u_buffer, dim_data{Ky_x, Kyx_ctx}, dim_data{Ky_y, Kyy_ctx});
        lin::solve_with_factorized(Ky_x, u, Kyx_ctx);
        auto F2 = lin::cyclic_transpose(u, u_buffer.data());
        lin::solve_with_factorized(Ky_y, F2, Kyy_ctx);
        lin::cyclic_transpose(F2, u);
    }

    void after_step(int iter, double t) override {
        if ((iter + 1) % save_every == 0) {
            std::cout << "Step " << (iter + 1) << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        //     analyze(iter, t + 0.5 * steps.dt);
        }

        // auto s = t / 150;
        // auto phase = sin(s) + 0.5 * sin(2.3*s);
        // wind_angle = M_PI / 3 * phase + 1.5 * M_PI / 4;

        // wind = { wind_speed * cos(wind_angle), wind_speed * sin(wind_angle) };
        prepare_implicit_matrices();

    }

    double grad_dot(point_type a, value_type u) const {
        return a[0] * u.dx + a[1] * u.dy;
    }

    value_type eval_basis(index_type e, index_type q, index_type a, const dimension& x, const dimension& y) const  {
        auto loc = dof_global_to_local(e, a, x, y);

        const auto& bx = x.basis;
        const auto& by = y.basis;

        double B1  = bx.b[e[0]][q[0]][0][loc[0]];
        double B2  = by.b[e[1]][q[1]][0][loc[1]];
        double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
        double dB2 = by.b[e[1]][q[1]][1][loc[1]];

        double v = B1 * B2;
        double dxv = dB1 *  B2;
        double dyv =  B1 * dB2;

        return { v, dxv, dyv };
    }

    value_type eval(const vector_type& v, index_type e, index_type q, const dimension& x, const dimension& y) const {
        value_type u{};
        for (auto b : dofs_on_element(e, x, y)) {
            double c = v(b[0], b[1]);
            value_type B = eval_basis(e, q, b, x, y);
            u += c * B;
        }
        return u;
    }

    index_range elements(const dimension& x, const dimension& y) const {
        return util::product_range<index_type>(x.element_indices(), y.element_indices());
    }

    index_range quad_points(const dimension& x, const dimension& y) const {
        auto rx = boost::counting_range(0, x.basis.quad_order);
        auto ry = boost::counting_range(0, y.basis.quad_order);
        return util::product_range<index_type>(rx, ry);
    }

    index_range dofs_on_element(index_type e, const dimension& x, const dimension& y) const {
        auto rx = x.basis.dof_range(e[0]);
        auto ry = y.basis.dof_range(e[1]);
        return util::product_range<index_type>(rx, ry);
    }

    index_type dof_global_to_local(index_type e, index_type a, const dimension& x, const dimension& y) const {
        const auto& bx = x.basis;
        const auto& by = y.basis;
        return {{ a[0] - bx.first_dof(e[0]), a[1] - by.first_dof(e[1]) }};
    }

    void update_global_rhs(vector_type& global, const vector_type& local, index_type e,
                           const dimension& x, const dimension& y) const {
        for (auto a : dofs_on_element(e, x, y)) {
            auto loc = dof_global_to_local(e, a, x, y);
            global(a[0], a[1]) += local(loc[0], loc[1]);
        }
    }


    void compute_rhs_x() {
        zero(rhs1);

        executor.for_each(elements(Vx, Uy), [&](index_type e) {
            auto U = vector_type{{ Vx.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};
            auto h = 0.5 * steps.dt;

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Uy)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval(u_prev, e, q, Ux, Uy);

                for (auto a : dofs_on_element(e, Vx, Uy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Uy);
                    value_type v = eval_basis(e, q, a, Vx, Uy);

                    double conv_term = wind[1] * u.dy * v.val;

                    double gradient_prod = u.dy * v.dy;
                    double e = emission(x[0], x[1], u.val) * v.val;
                    double val = u.val * v.val + h * (- conv_term - c_diff[1] * gradient_prod + e);

                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs1, U, e, Vx, Uy);
            });
        });
    }

    void compute_rhs_y() {
        zero(rhs2);

        executor.for_each(elements(Ux, Vy), [&](index_type e) {
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto h = 0.5 * steps.dt;

            double J = jacobian(e);
            for (auto q : quad_points(Ux, Vy)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval(u_prev, e, q, Ux, Uy);

                for (auto a : dofs_on_element(e, Ux, Vy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Vy);
                    value_type v = eval_basis(e, q, a, Ux, Vy);

                    double conv_term = wind[0] * u.dx * v.val;

                    double gradient_prod = u.dx * v.dx;
                    double e = emission(x[0], x[1], u.val) * v.val;
                    double val = u.val * v.val + h * (- conv_term - c_diff[0] * gradient_prod + e);
                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs2, U, e, Ux, Vy);
            });
        });
    }

    void analyze(int iter, double t) {
        double total = 0.0;
        double emission_rate = 0.0;
        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e);
            for (auto q : quad_points(Ux, Ux)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type u = eval(u_prev, e, q, Ux, Uy);

                total += u.val * w * J;
                emission_rate += emission(x[0], x[1], u.val) * w * J;
            }
        }

        using namespace std;

        total /= 1000; // g -> kg
        emission_rate /= 1000; // g -> kg

        auto emitted = t * emission_rate; // kg
        auto area = 5000.0 * 5000.0; // m^2
        auto initial = area * ambient / 1000; // kg
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



#endif /* PROBLEMS_POLLUTION_V2_DPG_2D_HPP_ */
