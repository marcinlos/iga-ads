#ifndef PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_
#define PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_

#include "problems/maxwell/problems.hpp"
#include "problems/maxwell/state.hpp"
#include "ads/executor/galois.hpp"
#include "ads/simulation.hpp"
#include "ads/form_matrix.hpp"
#include "ads/output_manager.hpp"


namespace ads {

class maxwell_galerkin : public simulation_3d {
private:
    using Base = simulation_3d;
    using Problem = maxwell_manufactured1;

    galois_executor executor{6};

    dimension UE1x, UE1y, UE1z;
    dimension UE2x, UE2y, UE2z;
    dimension UE3x, UE3y, UE3z;

    dimension UH1x, UH1y, UH1z;
    dimension UH2x, UH2y, UH2z;
    dimension UH3x, UH3y, UH3z;

    dimension Vx, Vy, Vz;

    state prev, half, now;

    lin::band_matrix Bx, By, Bz;
    lin::solver_ctx Bx_ctx, By_ctx, Bz_ctx;

    Problem problem{1, 1};

    ads::output_manager<3> output;

public:
    maxwell_galerkin(const config_3d& config)
    : Base{config}
    , UE1x{x}, UE1y{y}, UE1z{z}
    , UE2x{x}, UE2y{y}, UE2z{z}
    , UE3x{x}, UE3y{y}, UE3z{z}
    , UH1x{x}, UH1y{y}, UH1z{z}
    , UH2x{x}, UH2y{y}, UH2z{z}
    , UH3x{x}, UH3y{y}, UH3z{z}
    , Vx{x}
    , Vy{y}
    , Vz{z}
    , prev{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
    , half{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
    , now{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
    , Bx{Vx.p, Vx.p, Vx.dofs()}
    , By{Vy.p, Vy.p, Vy.dofs()}
    , Bz{Vz.p, Vz.p, Vz.dofs()}
    , Bx_ctx{Bx}
    , By_ctx{By}
    , Bz_ctx{Bz}
    , output{Vx.B, Vy.B, Vz.B, 50}
    { }

private:

    void fix_dof(int k, const dimension& dim, lin::band_matrix& K) {
        int last = dim.dofs() - 1;
        for (int i = clamp(k - dim.p, 0, last); i <= clamp(k + dim.p, 0, last); ++ i) {
            K(k, i) = 0;
        }
        K(k, k) = 1;
    }

    void prepare_matrices() {
        auto zero = [](auto& dim) { dim.fix_left(); dim.fix_right(); };

        // E x n = 0
        zero(UE1y); zero(UE1z);
        zero(UE2x); zero(UE2z);
        zero(UE3x); zero(UE3y);

        // H * n = 0
        zero(UH1x);
        zero(UH2y);
        zero(UH3z);

        UE1x.factorize_matrix();
        UE1y.factorize_matrix();
        UE1z.factorize_matrix();
        UE2x.factorize_matrix();
        UE2y.factorize_matrix();
        UE2z.factorize_matrix();
        UE3x.factorize_matrix();
        UE3y.factorize_matrix();
        UE3z.factorize_matrix();

        UH1x.factorize_matrix();
        UH1y.factorize_matrix();
        UH1z.factorize_matrix();
        UH2x.factorize_matrix();
        UH2y.factorize_matrix();
        UH2z.factorize_matrix();
        UH3x.factorize_matrix();
        UH3y.factorize_matrix();
        UH3z.factorize_matrix();

        Vx.factorize_matrix();
        Vy.factorize_matrix();
        Vz.factorize_matrix();


        Bx.zero();
        By.zero();
        Bz.zero();

        auto tau = steps.dt;
        auto h = tau * tau / (4 * problem.eps * problem.mu);
        auto form = [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; };

        form_matrix(Bx, Vx.basis, form);
        form_matrix(By, Vy.basis, form);
        form_matrix(Bz, Vz.basis, form);

        fix_dof(0, Vx, Bx); fix_dof(Vx.dofs() - 1, Vx, Bx);
        fix_dof(0, Vy, By); fix_dof(Vy.dofs() - 1, Vy, By);
        fix_dof(0, Vz, Bz); fix_dof(Vz.dofs() - 1, Vz, Bz);

        lin::factorize(Bx, Bx_ctx);
        lin::factorize(By, By_ctx);
        lin::factorize(Bz, Bz_ctx);
    }

    void before() override {
        prepare_matrices();

        auto project = [&](auto& rhs, auto& x, auto& y, auto& z, auto fun) {
            auto f = [&](double x, double y, double z) { return fun({x, y, z}); };
            auto buffer = vector_type{{ x.dofs(), y.dofs(), z.dofs() }};
            compute_projection(rhs, x.basis, y.basis, z.basis, f);
            ads_solve(rhs, buffer, x.data(), y.data(), z.data());
        };

        project(now.E1, UE1x, UE1y, UE1z, problem.E1_val_at(0));
        project(now.E2, UE2x, UE2y, UE2z, problem.E2_val_at(0));
        project(now.E3, UE3x, UE3y, UE3z, problem.E3_val_at(0));

        project(now.H1, UH1x, UH1y, UH1z, problem.H1_val_at(0));
        project(now.H2, UH2x, UH2y, UH2z, problem.H2_val_at(0));
        project(now.H3, UH3x, UH3y, UH3z, problem.H3_val_at(0));

        after_step(-1, -steps.dt);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
    }

    void zero_sides(const char dims[], vector_type& rhs,
            const dimension& Ux, const dimension& Uy, const dimension& Uz) const {
        for (char dim : std::string{dims}) {
            if (dim == 'x') {
                for (int i = 0; i < Uy.dofs(); ++ i) {
                    for (int j = 0; j < Uz.dofs(); ++ j) {
                        rhs(0, i, j) = 0;
                        rhs(Ux.dofs() - 1, i, j) = 0;
                    }
                }
            } else if (dim == 'y') {
                for (int i = 0; i < Ux.dofs(); ++ i) {
                    for (int j = 0; j < Uz.dofs(); ++ j) {
                        rhs(i, 0, j) = 0;
                        rhs(i, Uy.dofs() - 1, j) = 0;
                    }
                }
            } else if (dim == 'z') {
                for (int i = 0; i < Ux.dofs(); ++ i) {
                    for (int j = 0; j < Uy.dofs(); ++ j) {
                        rhs(i, j, 0) = 0;
                        rhs(i, j, Uz.dofs() - 1) = 0;
                    }
                }
            }
        }
    }

    void step(int /*iter*/, double /*t*/) override {
        auto tau = steps.dt;
        auto eps = problem.eps;
        auto mu = problem.mu;
        auto a = tau / (2 * eps);
        auto b = tau * tau / (4 * eps);
        auto c = tau / (2 * mu);

        using shape = std::array<std::size_t, 3>;
        auto basic_shape = shape{Vx.dofs(), Vy.dofs(), Vz.dofs()};
        // Buffer large enough for all the RHS
        auto buffer = vector_type{basic_shape};

        auto shape_E1 = basic_shape;
        auto shape_E2 = basic_shape;
        auto shape_E3 = basic_shape;

        auto shape_H1 = basic_shape;
        auto shape_H2 = basic_shape;
        auto shape_H3 = basic_shape;

        constexpr int X = 0;
        constexpr int Y = 1;
        constexpr int Z = 2;

        // First substep - E
        auto rhs_E1 = vector_type{shape_E1};
        compute_rhs(rhs_E1, prev, prev.E1, prev.E2, prev.E3, UE1x, UE1y, UE1z, [=](auto E, auto, auto H, auto v) {
            return (E[X].val + a * (H[Z].dy - H[Y].dz)) * v.val + b * E[Y].dx * v.dy;
        });
        zero_sides("yz", rhs_E1, UE1x, UE1y, UE1z);
        ads_solve(rhs_E1, buffer, UE1x.data(), dim_data{By, By_ctx}, UE1z.data());

        auto rhs_E2 = vector_type{shape_E2};
        compute_rhs(rhs_E2, prev, prev.E1, prev.E2, prev.E3, UE2x, UE2y, UE2z, [=](auto E, auto, auto H, auto v) {
            return (E[Y].val + a * (H[X].dz - H[Z].dx)) * v.val + b * E[Z].dy * v.dz;
        });
        zero_sides("xz", rhs_E2, UE2x, UE2y, UE2z);
        ads_solve(rhs_E2, buffer, UE2x.data(), UE2y.data(), dim_data{Bz, Bz_ctx});

        auto rhs_E3 = vector_type{shape_E3};
        compute_rhs(rhs_E3, prev, prev.E1, prev.E2, prev.E3, UE3x, UE3y, UE3z, [=](auto E, auto, auto H, auto v) {
            return (E[Z].val + a * (H[Y].dx - H[X].dy)) * v.val + b * E[X].dz * v.dx;
        });
        zero_sides("xy", rhs_E3, UE3x, UE3y, UE3z);
        ads_solve(rhs_E3, buffer, dim_data{Bx, Bx_ctx}, UE3y.data(), UE3z.data());

        // First substep - H
        auto rhs_H1 = vector_type{shape_H1};
        compute_rhs(rhs_H1, prev, rhs_E1, rhs_E2, rhs_E3, UH1x, UH1y, UH1z, [=](auto E, auto En, auto H, auto v) {
            return (H[X].val - c * (E[Z].dy - En[Y].dz)) * v.val;
        });
        zero_sides("x", rhs_H1, UH1x, UH1y, UH1z);
        ads_solve(rhs_H1, buffer, UH1x.data(), UH1y.data(), UH1z.data());

        auto rhs_H2 = vector_type{shape_H2};
        compute_rhs(rhs_H2, prev, rhs_E1, rhs_E2, rhs_E3, UH2x, UH2y, UH2z, [=](auto E, auto En, auto H, auto v) {
            return (H[Y].val - c * (E[X].dz - En[Z].dx)) * v.val;
        });
        zero_sides("y", rhs_H2, UH2x, UH2y, UH2z);
        ads_solve(rhs_H2, buffer, UH2x.data(), UH2y.data(), UH2z.data());

        auto rhs_H3 = vector_type{shape_H3};
        compute_rhs(rhs_H3, prev, rhs_E1, rhs_E2, rhs_E3, UH3x, UH3y, UH3z, [=](auto E, auto En, auto H, auto v) {
            return (H[Z].val - c * (E[Y].dx - En[X].dy)) * v.val;
        });
        zero_sides("z", rhs_H3, UH3x, UH3y, UH3z);
        ads_solve(rhs_H3, buffer, UH3x.data(), UH3y.data(), UH3z.data());

        prev.E1 = rhs_E1;
        prev.E2 = rhs_E2;
        prev.E3 = rhs_E3;

        prev.H1 = rhs_H1;
        prev.H2 = rhs_H2;
        prev.H3 = rhs_H3;

        // Second substep - E
        compute_rhs(now.E1, prev, prev.E1, prev.E2, prev.E3, UE1x, UE1y, UE1z, [=](auto E, auto, auto H, auto v) {
            return (E[X].val + a * (H[Z].dy - H[Y].dz)) * v.val + b * E[Z].dx * v.dz;
        });
        zero_sides("yz", now.E1, UE1x, UE1y, UE1z);
        ads_solve(now.E1, buffer, UE1x.data(), UE1y.data(), dim_data{Bz, Bz_ctx});

        compute_rhs(now.E2, prev, prev.E1, prev.E2, prev.E3, UE2x, UE2y, UE2z, [=](auto E, auto, auto H, auto v) {
            return (E[Y].val + a * (H[X].dz - H[Z].dx)) * v.val + b * E[X].dy * v.dx;
        });
        zero_sides("xz", now.E2, UE2x, UE2y, UE2z);
        ads_solve(now.E2, buffer, dim_data{Bx, Bx_ctx}, UE2y.data(), UE2z.data());

        compute_rhs(now.E3, prev, prev.E1, prev.E2, prev.E3, UE3x, UE3y, UE3z, [=](auto E, auto, auto H, auto v) {
            return (E[Z].val + a * (H[Y].dx - H[X].dy)) * v.val + b * E[Y].dz * v.dy;
        });
        zero_sides("xy", now.E3, UE3x, UE3y, UE3z);
        ads_solve(now.E3, buffer, UE3x.data(), dim_data{By, By_ctx}, UE3z.data());

        // Second substep - H
        compute_rhs(now.H1, prev, now.E1, now.E2, now.E3, UH1x, UH1y, UH1z, [=](auto E, auto En, auto H, auto v) {
            return (H[X].val - c * (En[Z].dy - E[Y].dz)) * v.val;
        });
        zero_sides("x", now.H1, UH1x, UH1y, UH1z);
        ads_solve(now.H1, buffer, UH1x.data(), UH1y.data(), UH1z.data());

        compute_rhs(now.H2, prev, now.E1, now.E2, now.E3, UH2x, UH2y, UH2z, [=](auto E, auto En, auto H, auto v) {
            return (H[Y].val - c * (En[X].dz - E[Z].dx)) * v.val;
        });
        zero_sides("y", now.H2, UH2x, UH2y, UH2z);
        ads_solve(now.H2, buffer, UH2x.data(), UH2y.data(), UH2z.data());

        compute_rhs(now.H3, prev, now.E1, now.E2, now.E3, UH3x, UH3y, UH3z, [=](auto E, auto En, auto H, auto v) {
            return (H[Z].val - c * (En[Y].dx - E[X].dy)) * v.val;
        });
        zero_sides("z", now.H3, UH3x, UH3y, UH3z);
        ads_solve(now.H3, buffer, UH3x.data(), UH3y.data(), UH3z.data());
    }

    template <typename RHS, typename Form>
    void compute_rhs(RHS& rhs, const state& prev,
            const vector_type& E1_new, const vector_type& E2_new, const vector_type& E3_new,
            const dimension& Vx, const dimension& Vy, const dimension& Vz,
            Form&& form) {
        zero(rhs);
        using shape = std::array<std::size_t, 3>;
        auto shape_loc = shape{Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element(), Vz.basis.dofs_per_element()};

        executor.for_each(elements(Vx, Vy, Vz), [&](auto e) {
            auto loc = vector_type{shape_loc};

            auto J = jacobian(e, Vx, Vy, Vz);
            for (auto q : quad_points(Vx, Vy, Vz)) {
                auto W = weigth(q, Vx, Vy, Vz);

                using vec = std::array<value_type, 3>;

                auto E1 = eval(prev.E1, e, q, UE1x, UE1y, UE1z);
                auto E2 = eval(prev.E2, e, q, UE2x, UE2y, UE2z);
                auto E3 = eval(prev.E3, e, q, UE3x, UE3y, UE3z);
                auto E = vec{E1, E2, E3};

                auto E1_n = eval(E1_new, e, q, UE1x, UE1y, UE1z);
                auto E2_n = eval(E2_new, e, q, UE2x, UE2y, UE2z);
                auto E3_n = eval(E3_new, e, q, UE3x, UE3y, UE3z);
                auto E_n = vec{E1_n, E2_n, E3_n};

                auto H1 = eval(prev.H1, e, q, UH1x, UH1y, UH1z);
                auto H2 = eval(prev.H2, e, q, UH2x, UH2y, UH2z);
                auto H3 = eval(prev.H3, e, q, UH3x, UH3y, UH3z);
                auto H = vec{H1, H2, H3};

                for (auto a : dofs_on_element(e, Vx, Vy, Vz)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy, Vz);
                    auto v = eval_basis(e, q, a, Vx, Vy, Vz);

                    auto val = form(E, E_n, H, v);
                    loc(aa[0], aa[1], aa[2]) += val * W * J;
                }
            }
            executor.synchronized([&]{
                update_global_rhs(rhs, loc, e, Vx, Vy, Vz);
            });
        });
    }

    void after_step(int iter, double t) override {
        auto i = iter + 1;
        auto tt = t + steps.dt;

        if (i % 10 == 0) {
            output.to_file("out_%d.vti", i,
                    output.evaluate(now.E1),
                    output.evaluate(now.E2),
                    output.evaluate(now.E3),
                    output.evaluate(now.H1),
                    output.evaluate(now.H2),
                    output.evaluate(now.H3));
        }

        auto E1_norm_L2 = normL2(now.E1, UE1x, UE1y, UE1z);
        auto E2_norm_L2 = normL2(now.E2, UE2x, UE2y, UE2z);
        auto E3_norm_L2 = normL2(now.E3, UE3x, UE3y, UE3z);
        auto E_norm_L2 = std::sqrt(E1_norm_L2 * E1_norm_L2 + E2_norm_L2 * E2_norm_L2 + E3_norm_L2 * E3_norm_L2);

        auto E1_err_L2 = errorL2(now.E1, UE1x, UE1y, UE1z, problem.E1_at(tt));
        auto E2_err_L2 = errorL2(now.E2, UE2x, UE2y, UE2z, problem.E2_at(tt));
        auto E3_err_L2 = errorL2(now.E3, UE3x, UE3y, UE3z, problem.E3_at(tt));
        auto E_err_L2 = std::sqrt(E1_err_L2 * E1_err_L2 + E2_err_L2 * E2_err_L2 + E3_err_L2 * E3_err_L2);

        auto E1_norm_H1 = normH1(now.E1, UE1x, UE1y, UE1z);
        auto E2_norm_H1 = normH1(now.E2, UE2x, UE2y, UE2z);
        auto E3_norm_H1 = normH1(now.E3, UE3x, UE3y, UE3z);
        auto E_norm_H1 = std::sqrt(E1_norm_H1 * E1_norm_H1 + E2_norm_H1 * E2_norm_H1 + E3_norm_H1 * E3_norm_H1);

        auto E1_err_H1 = errorH1(now.E1, UE1x, UE1y, UE1z, problem.E1_at(tt));
        auto E2_err_H1 = errorH1(now.E2, UE2x, UE2y, UE2z, problem.E2_at(tt));
        auto E3_err_H1 = errorH1(now.E3, UE3x, UE3y, UE3z, problem.E3_at(tt));
        auto E_err_H1 = std::sqrt(E1_err_H1 * E1_err_H1 + E2_err_H1 * E2_err_H1 + E3_err_H1 * E3_err_H1);

        auto rot_E = norm_rot(now.E1, now.E2, now.E3, Vx, Vy, Vz);
        auto div_E = norm_div(now.E1, now.E2, now.E3, Vx, Vy, Vz);

        auto H1_norm_L2 = normL2(now.H1, UH1x, UH1y, UH1z);
        auto H2_norm_L2 = normL2(now.H2, UH2x, UH2y, UH2z);
        auto H3_norm_L2 = normL2(now.H3, UH3x, UH3y, UH3z);
        auto H_norm_L2 = std::sqrt(H1_norm_L2 * H1_norm_L2 + H2_norm_L2 * H2_norm_L2 + H3_norm_L2 * H3_norm_L2);

        auto H1_err_L2 = errorL2(now.H1, UH1x, UH1y, UH1z, problem.H1_at(tt));
        auto H2_err_L2 = errorL2(now.H2, UH2x, UH2y, UH2z, problem.H2_at(tt));
        auto H3_err_L2 = errorL2(now.H3, UH3x, UH3y, UH3z, problem.H3_at(tt));
        auto H_err_L2 = std::sqrt(H1_err_L2 * H1_err_L2 + H2_err_L2 * H2_err_L2 + H3_err_L2 * H3_err_L2);

        auto H1_norm_H1 = normH1(now.H1, UH1x, UH1y, UH1z);
        auto H2_norm_H1 = normH1(now.H2, UH2x, UH2y, UH2z);
        auto H3_norm_H1 = normH1(now.H3, UH3x, UH3y, UH3z);
        auto H_norm_H1 = std::sqrt(H1_norm_H1 * H1_norm_H1 + H2_norm_H1 * H2_norm_H1 + H3_norm_H1 * H3_norm_H1);

        auto H1_err_H1 = errorH1(now.H1, UH1x, UH1y, UH1z, problem.H1_at(tt));
        auto H2_err_H1 = errorH1(now.H2, UH2x, UH2y, UH2z, problem.H2_at(tt));
        auto H3_err_H1 = errorH1(now.H3, UH3x, UH3y, UH3z, problem.H3_at(tt));
        auto H_err_H1 = std::sqrt(H1_err_H1 * H1_err_H1 + H2_err_H1 * H2_err_H1 + H3_err_H1 * H3_err_H1);

        auto rot_H = norm_rot(now.H1, now.H2, now.H1, Vx, Vy, Vz);
        auto div_H = norm_div(now.H1, now.H2, now.H1, Vx, Vy, Vz);

        std::cout << "After step " << i << ", t = " << tt << std::endl;
        std::cout << "  |E|     = " << E_norm_L2 << "  " << E_norm_H1 << std::endl;
        std::cout << "  |rot E| = " << rot_E << std::endl;
        std::cout << "  |div E| = " << div_E << std::endl;
        std::cout << "  E err = " << E_err_L2 << "  " << E_err_H1 << std::endl;
        std::cout << "    rel = " << E_err_L2 / E_norm_L2 * 100 << "%  "
                                  << E_err_H1 / E_norm_H1 * 100 << "%" << std::endl;

        std::cout << "  |H| = " << H_norm_L2 << "  " << H_norm_H1 << std::endl;
        std::cout << "  |rot H| = " << rot_H << std::endl;
        std::cout << "  |div H| = " << div_H << std::endl;
        std::cout << "  H err = " << H_err_L2 << "  " << H_err_H1 << std::endl;
        std::cout << "    rel = " << H_err_L2 / H_norm_L2 * 100 << "%  "
                                  << H_err_H1 / H_norm_H1 * 100 << "%" << std::endl;
    }

};

}

#endif // ifndef PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_
