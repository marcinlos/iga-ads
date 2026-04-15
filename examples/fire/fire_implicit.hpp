// SPDX-FileCopyrightText: 2026 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"
#include "manufactured.hpp"
#include "params.hpp"

using namespace ads;

enum class scheme { FE, BE, CN, peaceman_rachford, strang_BE, strang_CN };

inline double falloff(double r, double R, double t) {
    if (t < r)
        return 1.0;
    if (t > R)
        return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

inline double bump(double r, double R, double x, double y) {
    double dx = x - 50;
    double dy = y - 50;
    double t = std::sqrt(dx * dx + dy * dy) / 100;
    return falloff(r / 200, R / 200, t);
}

class fire_implicit : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;

    // New stuff
    lin::band_matrix MVx, MVy;
    lin::band_matrix MUx, MUy;

    lin::band_matrix KVx, KVy;
    lin::band_matrix KUx, KUy;

    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;

    lin::dense_matrix MUUx, MUUy;
    lin::dense_matrix KUUx, KUUy;
    lin::dense_matrix AUUx, AUUy;

    struct residuum {
        vector_type data;
        const dimension* Vx;
        const dimension* Vy;

        residuum(vector_type data, const dimension* Vx, const dimension* Vy)
        : data{std::move(data)}
        , Vx{Vx}
        , Vy{Vy} { }
    };

    vector_type u;
    residuum r;
    vector_type solver_buffer;
    std::vector<double> full_rhs;

    vector_type fuel, fuel_prev;

    int save_every = 1;

    // Model parameters
    point_type beta{{0, 0}};
    fire_params params;

    mumps::solver solver;

    galois_executor executor;
    output_manager<2> output;

    scheme method;
    solution manufactured;

public:
    fire_implicit(dimension const& trial_x, dimension const& trial_y,  //
                  dimension const& test_x, dimension const& test_y,    //
                  fire_params const& params, int threads, scheme method,
                  timesteps_config const& steps)
    : Base{test_x, test_y, steps}
    , Ux{trial_x}
    , Uy{trial_y}
    , Vx{x}
    , Vy{y}
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , MUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , KUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , MUVx{Vx.dofs(), Ux.dofs()}
    , MUVy{Vy.dofs(), Uy.dofs()}
    , KUVx{Vx.dofs(), Ux.dofs()}
    , KUVy{Vy.dofs(), Uy.dofs()}
    , AUVx{Vx.dofs(), Ux.dofs()}
    , AUVy{Vy.dofs(), Uy.dofs()}
    , MUUx{Ux.dofs(), Ux.dofs()}
    , MUUy{Uy.dofs(), Uy.dofs()}
    , KUUx{Ux.dofs(), Ux.dofs()}
    , KUUy{Uy.dofs(), Uy.dofs()}
    , AUUx{Ux.dofs(), Ux.dofs()}
    , AUUy{Uy.dofs(), Uy.dofs()}
    , u{{Ux.dofs(), Uy.dofs()}}
    , r{vector_type{{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy}
    , solver_buffer{{Ux.dofs(), Uy.dofs()}}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , fuel{{Ux.dofs(), Uy.dofs()}}
    , fuel_prev{{Ux.dofs(), Uy.dofs()}}
    , params{params}
    , executor{threads}
    , output{Ux.B, Uy.B, 300}
    , method{method}
    , manufactured{params, beta[0], beta[1]} { }

private:
    struct matrix_set {
        using band_matrix_ref = lin::band_matrix&;
        using dense_matrix_ref = lin::dense_matrix&;

        band_matrix_ref MVx, MVy, KVx, KVy;
        dense_matrix_ref MUVx, MUVy, KUVx, KUVy, AUVx, AUVy;
    };

    void matrix_1d(lin::band_matrix& M, basis_data const& d, double diff, double b_val, double c) {
        for (element_id e = 0; e < d.elements; ++e) {
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
                        auto val = va * vb + c * (diff * da * db + b_val * va * db);
                        M(ia, ib) += val * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    struct kron_matrix {
        lin::band_matrix Ax;
        lin::band_matrix Ay;

        kron_matrix(int px, int py, int dofs_x, int dofs_y)
        : Ax{px, px, dofs_x}
        , Ay{py, py, dofs_y} { }
    };

    auto assemble_problem_ads(double cx, double cy) -> kron_matrix {
        if (cx != 0.0 && cy != 0) {
            std::cerr << "Invalid scheme" << std::endl;
            std::exit(1);
        }

        auto matrices = kron_matrix{Ux.p, Uy.p, Ux.dofs(), Uy.dofs()};

        auto const diff = params.kappa / (params.rho * params.cp);
        auto const wind = params.cw / params.cp;
        auto const bx = wind * beta[0];
        auto const by = wind * beta[1];

        matrix_1d(matrices.Ax, Ux.basis, diff, bx, cx);
        matrix_1d(matrices.Ay, Uy.basis, diff, by, cy);

        return matrices;
    }

    void assemble_problem(mumps::problem& problem, double cx, double cy, double sx, double sy,
                          const dimension& Vx, const dimension& Vy, const matrix_set& M) {
        auto N = Vx.dofs() * Vy.dofs();

        // Gram matrix
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : overlapping_dofs(i, Vx, Vy)) {
                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                auto const MxMy = kron(M.MVx, M.MVy, i, j);
                auto const KxMy = kron(M.KVx, M.MVy, i, j);
                auto const MxKy = kron(M.MVx, M.KVy, i, j);

                double val = MxMy + sx * KxMy + sy * MxKy;
                problem.add(ii, jj, val);
            }
        }

        auto const diff = params.kappa / (params.rho * params.cp);
        auto const wind = params.cw / params.cp;
        auto const bx = wind * beta[0];
        auto const by = wind * beta[1];

        // B, B^T
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : dofs(Ux, Uy)) {
                auto const MxMy = kron(M.MUVx, M.MUVy, i, j);
                auto const KxMy = kron(M.KUVx, M.MUVy, i, j);
                auto const MxKy = kron(M.MUVx, M.KUVy, i, j);
                auto const AxMy = kron(M.AUVx, M.MUVy, i, j);
                auto const MxAy = kron(M.MUVx, M.AUVy, i, j);

                double Lx = diff * KxMy + bx * AxMy;
                double Ly = diff * MxKy + by * MxAy;
                double val = MxMy + cx * Lx + cy * Ly;

                if (val != 0) {
                    int ii = linear_index(i, Vx, Vy) + 1;
                    int jj = linear_index(j, Ux, Uy) + 1;

                    problem.add(ii, N + jj, -val);
                    problem.add(N + jj, ii, val);
                }
            }
        }
    }

    matrix_set matrices(bool x_refined, bool y_refined) {
        if (x_refined && y_refined) {
            return {MVx, MVy, KVx, KVy, MUVx, MUVy, KUVx, KUVy, AUVx, AUVy};
        } else if (x_refined) {
            return {MVx, MUy, KVx, KUy, MUVx, MUUy, KUVx, KUUy, AUVx, AUUy};
        } else if (y_refined) {
            return {MUx, MVy, KUx, KVy, MUUx, MUVy, KUUx, KUVy, AUUx, AUVy};
        } else {
            return {MUx, MUy, KUx, KUy, MUUx, MUUy, KUUx, KUUy, AUUx, AUUy};
        }
    }

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);

        gram_matrix_1d(MUx, Ux.basis);
        gram_matrix_1d(MUy, Uy.basis);

        gram_matrix_1d(MUUx, Ux.basis, Ux.basis);
        gram_matrix_1d(MUUy, Uy.basis, Uy.basis);

        gram_matrix_1d(MUVx, Ux.basis, Vx.basis);
        gram_matrix_1d(MUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);

        stiffness_matrix_1d(KUx, Ux.basis);
        stiffness_matrix_1d(KUy, Uy.basis);

        stiffness_matrix_1d(KUVx, Ux.basis, Vx.basis);
        stiffness_matrix_1d(KUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KUUx, Ux.basis, Ux.basis);
        stiffness_matrix_1d(KUUy, Uy.basis, Uy.basis);

        advection_matrix_1d(AUVx, Ux.basis, Vx.basis);
        advection_matrix_1d(AUVy, Uy.basis, Vy.basis);

        advection_matrix_1d(AUUx, Ux.basis, Ux.basis);
        advection_matrix_1d(AUUy, Uy.basis, Uy.basis);
    }

    template <typename RHS>
    void const_bc(RHS& u, dimension& Ux, dimension& Uy, double value) {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = value; });
    }

    double init_state(double x, double y) {
        return manufactured.exact(x, y, 0).val;
    }

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();

        zero(r.data);
        zero(u);

        auto init = [this](double x, double y) { return init_state(x, y); };
        compute_projection(u, Ux.basis, Uy.basis, init);
        ads_solve(u, solver_buffer, Ux.data(), Uy.data());
        output.to_file(u, "out_0.data");

        auto fuel_init = [](double /*x*/, double /*y*/) { return 1.0; };
        projection(fuel, fuel_init);
        ads_solve(fuel, solver_buffer, Ux.data(), Uy.data());
        output.to_file(fuel, "fuel_0.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        // no need to swap u, this is handled differently
        swap(fuel, fuel_prev);
    }

    void copy_solution(const vector_view& u_rhs, const vector_view& r_rhs, vector_type& u) {
        for (auto i = 0; i < Ux.dofs(); ++i) {
            for (auto j = 0; j < Uy.dofs(); ++j) {
                u(i, j) = u_rhs(i, j);
            }
        }
        r = residuum{vector_type{{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy};
        for (auto i = 0; i < Vx.dofs(); ++i) {
            for (auto j = 0; j < Vy.dofs(); ++j) {
                r.data(i, j) = r_rhs(i, j);
            }
        }
    }

    template <typename Fun>
    void substep(vector_type& u, bool x_refine, bool y_refine, double Lx_lhs, double Ly_lhs,
                 double Lx_rhs, double Ly_rhs, double dt, Fun&& f) {
        dimension& Vx = x_refine ? this->Vx : Ux;
        dimension& Vy = y_refine ? this->Vy : Uy;

        double sx = x_refine ? 0 : 1;
        double sy = y_refine ? 0 : 1;

        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);
        compute_rhs(Lx_rhs, Ly_rhs, Vx, Vy, r_rhs, u_rhs, dt, std::forward<Fun>(f));

        int size = Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs();
        mumps::problem problem(full_rhs.data(), size);
        assemble_problem(problem, Lx_lhs, Ly_lhs, sx, sy, Vx, Vy, matrices(x_refine, y_refine));
        solver.solve(problem);

        copy_solution(u_rhs, r_rhs, u);
    }

    template <typename Fun>
    void fast_substep(vector_type& u, bool x_refine, bool y_refine, double Lx_lhs, double Ly_lhs,
                      double Lx_rhs, double Ly_rhs, double dt, Fun&& f) {
        vector_type r_rhs{{Ux.dofs(), Uy.dofs()}};
        vector_type u_rhs{{Ux.dofs(), Uy.dofs()}};
        compute_rhs(Lx_rhs, Ly_rhs, Ux, Uy, r_rhs, u_rhs, dt, std::forward<Fun>(f));

        int size = Ux.dofs() * Uy.dofs() + Ux.dofs() * Uy.dofs();
        auto [Ax, Ay] = assemble_problem_ads(Lx_lhs, Ly_lhs);

        lin::solver_ctx ctx_x{Ax};
        lin::solver_ctx ctx_y{Ay};

        lin::factorize(Ax, ctx_x);
        lin::factorize(Ay, ctx_y);
        ads_solve(r_rhs, solver_buffer, dim_data{Ax, ctx_x}, dim_data{Ay, ctx_y});

        for (auto i = 0; i < Ux.dofs(); ++i) {
            for (auto j = 0; j < Uy.dofs(); ++j) {
                u(i, j) = -r_rhs(i, j);
            }
        }
    }

    void step(int /*iter*/, double t) override {
        auto dt = steps.dt;

        auto f = [&](point_type x, double s) { return manufactured.forcing(x[0], x[1], s); };
        auto F = [&](double s) { return [&, s](point_type x) { return f(x, s); }; };
        auto Favg = [&](double s1, double s2) {
            return [=, &f](point_type x) { return 0.5 * (f(x, s1) + f(x, s2)); };
        };
        auto zero = [&](point_type) { return 0; };

        // clang-format off
        if (method == scheme::FE) {
            fast_substep(u, true, true, 0, 0, dt, dt, dt, F(t));
        }
        if (method == scheme::BE) {
            substep(u, true, true, dt, dt, 0, 0, dt, F(t + dt));
        }
        if (method == scheme::CN) {
            substep(u, true, true,   dt/2, dt/2, -dt/2, -dt/2,   dt, Favg(t, t + dt));
        }
        if (method == scheme::peaceman_rachford) {
            fast_substep(u, true, true,   dt/2,    0,     0, -dt/2,   dt/2, F(t + dt/2));
            fast_substep(u, true, true,      0, dt/2, -dt/2,     0,   dt/2, F(t + dt/2));
        }
        if (method == scheme::strang_BE) {
            fast_substep(u, false, true,    dt/2,  0,   0, 0,   dt/2, F(t + dt/2));
            fast_substep(u, true,  false,      0, dt,   0, 0,      0, zero);
            fast_substep(u, false, true,    dt/2,  0,   0, 0,   dt/2, F(t + dt));
        }
        if (method == scheme::strang_CN) {
            fast_substep(u, false,  true,   dt/4,    0,   -dt/4,     0,   dt/2, Favg(t, t + dt/2));
            fast_substep(u, true,  false,      0, dt/2,       0, -dt/2,      0, zero);
            fast_substep(u, false,  true,   dt/4,    0,   -dt/4,     0,   dt/2, Favg(t + dt/2, t + dt));
        }
        // clang-format on

        update_fuel();
    }

    void update_fuel() {
        compute_rhs_fuel();
        ads_solve(fuel, solver_buffer, Ux.data(), Uy.data());
    }

    void after_step(int iter, double t) override {
        auto const i = iter + 1;
        if (i % save_every == 0) {
            // report_errors(i, t + steps.dt);
            // output.to_file(u, "out_%d.data", i);
            // output.to_file(fuel, "fuel_%d.data", i);
        }
    }

    void after() override {
        std::ofstream sol("solution.data");
        for (int i = 0; i < Ux.dofs(); ++i) {
            for (int j = 0; j < Uy.dofs(); ++j) {
                sol << i << " " << j << " " << u(i, j) << std::endl;
            }
        }
        report_errors(steps.step_count, 1.0);
    }

    template <typename VecR, typename VecU, typename Fun>
    void compute_rhs(double cx, double cy, const dimension& Vx, const dimension& Vy, VecR& r_rhs,
                     VecU& u_rhs, double dt, Fun&& F) {
        auto const ch = params.ch;
        auto const Ar = params.Ar;
        auto const rho = params.rho;
        auto const cp = params.cp;
        auto const cw = params.cw;
        auto const kappa = params.kappa;
        auto const xi = params.xi;
        auto const sigma = params.sigma;
        auto const hc = params.hc;
        auto const eps = params.eps;
        auto const delta_x = params.delta_x;
        auto const delta_z = params.delta_z;
        auto const T0 = params.T0;
        auto const Tig = params.Tig;
        auto const Ta = params.Ta;
        auto const M_param = params.M;
        auto const M1 = params.M1;

        auto const inv = 1.0 / (rho * cp);
        auto const diff = kappa * inv;
        auto const wind = cw / cp;
        auto const bx = wind * beta[0];
        auto const by = wind * beta[1];

        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element()}};
            auto U = vector_type{{Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element()}};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                auto x = point(e, q);

                value_type T = eval(u, e, q, Ux, Uy);
                value_type fuel = eval(fuel_prev, e, q, Ux, Uy);
                auto const Fx = F(x);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    double M = T.val * v.val;
                    double Lx = diff * T.dx * v.dx + bx * T.dx * v.val;
                    double Ly = diff * T.dy * v.dy + by * T.dy * v.val;

                    // double delta = T.val > Tig && fuel.val > 0.2 ? 1.0 : 0.0;
                    double delta = T.val > Tig /*&& fuel.val > 0.2*/ ? 1.0 : 0.0;
                    double r = delta * Ar * T.val * std::exp(-Ta / T.val);
                    // double Rc = -1e4 * rho * ch * hc * M_param / M1 * r;
                    double Rc = -1e3 * rho * ch * hc * M_param / M1 * r;
                    double Qw = 0;    // -rho * cw * (bx * T.dx + by * T.dy);
                    double qc = 0;    // -kappa * grad_dot(u, v);
                    double qd = 0.0;  // omitted
                    double qr = -4 * sigma * eps * delta_x * std::pow(T.val, 3) * grad_dot(T, v);
                    double Qconv = xi * (T0 - T.val);
                    double Qrz = sigma * eps / delta_z * (std::pow(T0, 4) - std::pow(T.val, 4));
                    // double rhs = (Rc + Qw + Qconv + Qrz + Fx) * v.val + qc + qd + qr;
                    double rhs = (Rc + Qconv + Qrz + Fx) * v.val + qr;

                    double lv = M + cx * Lx + cy * Ly + dt * inv * rhs;
                    double val = -lv;

                    R(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
    }

    void compute_rhs_fuel() {
        auto& rhs = fuel;
        zero(rhs);

        executor.for_each(elements(Ux, Uy), [&](index_type e) {
            auto F = vector_type{{Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element()}};

            double J = jacobian(e);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weight(q);

                value_type T = eval(u, e, q, Ux, Uy);
                value_type fuel = eval(fuel_prev, e, q, Ux, Uy);

                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type v = eval_basis(e, q, a, Ux, Uy);

                    auto const Ar = params.Ar;
                    auto const Tig = params.Tig;
                    auto const Ta = params.Ta;

                    double delta = T.val > Tig && fuel.val > 0.2 ? 1.0 : 0.0;
                    double r = delta * Ar * T.val * std::exp(-Ta / T.val);

                    double fval = -delta * 3e2 * r * fuel.val * v.val;
                    F(aa[0], aa[1]) += (fuel.val * v.val + steps.dt * fval) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, F, e, Ux, Uy); });
        });
    }

    auto errorL2(vector_type const& u, double t) const -> double {
        auto sol = [&](point_type x) { return manufactured.exact(x[0], x[1], t); };

        return Base::errorL2(u, Ux, Uy, sol);
    }

    auto errorH1(vector_type const& u, double t) const -> double {
        auto sol = [&](point_type x) { return manufactured.exact(x[0], x[1], t); };

        return Base::errorH1(u, Ux, Uy, sol);
    }

    auto rel_errorL2(vector_type const& u, double t) const -> double {
        auto sol = [&](point_type x) { return manufactured.exact(x[0], x[1], t); };
        return Base::errorL2(u, Ux, Uy, sol) / Base::normL2(Ux, Uy, sol) * 100;
    }

    auto rel_errorH1(vector_type const& u, double t) const -> double {
        auto sol = [&](point_type x) { return manufactured.exact(x[0], x[1], t); };
        return Base::errorH1(u, Ux, Uy, sol) / Base::normH1(Ux, Uy, sol) * 100;
    }

    auto report_errors(int iter, double t) const -> void {
        auto const e_L2 = errorL2(u, t);
        auto const e_H1 = errorH1(u, t);
        auto const rel_L2 = rel_errorL2(u, t);
        auto const rel_H1 = rel_errorH1(u, t);

        std::cout << "Step " << iter                    //
                  << "  t: " << t                       //
                  << "  L2: " << e_L2 << " " << rel_L2  //
                  << "  H1: " << e_H1 << " " << rel_H1  //
                  << std::endl;
    }
};
