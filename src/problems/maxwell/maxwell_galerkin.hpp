#ifndef PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_
#define PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_

#include "problems/maxwell/problems.hpp"
#include "problems/maxwell/state.hpp"
#include "ads/executor/galois.hpp"
#include "ads/simulation.hpp"
#include "ads/form_matrix.hpp"

namespace ads {

class maxwell_galerkin : public simulation_3d {
private:
    using Base = simulation_3d;
    using Problem = maxwell_manufactured1;

    galois_executor executor{6};

    dimension Ux, Uy, Uz;
    dimension Vx, Vy, Vz;

    state prev, half, now;

    lin::band_matrix Bx, By, Bz;
    lin::solver_ctx Bx_ctx, By_ctx, Bz_ctx;

    Problem problem{1, 1};

public:
    maxwell_galerkin(const config_3d& config)
    : Base{config}
    , Ux{x}
    , Uy{y}
    , Uz{z}
    , Vx{x}
    , Vy{y}
    , Vz{z}
    , prev{{Ux.dofs(), Uy.dofs(), Uz.dofs()}}
    , half{{Ux.dofs(), Uy.dofs(), Uz.dofs()}}
    , now{{Ux.dofs(), Uy.dofs(), Uz.dofs()}}
    , Bx{Ux.p, Ux.p, Ux.dofs()}
    , By{Uy.p, Uy.p, Uy.dofs()}
    , Bz{Uz.p, Uz.p, Uz.dofs()}
    , Bx_ctx{Bx}
    , By_ctx{By}
    , Bz_ctx{Bz}
    { }

private:

    void prepare_matrices() {
        Ux.factorize_matrix();
        Uy.factorize_matrix();
        Uz.factorize_matrix();

        Vx.factorize_matrix();
        Vy.factorize_matrix();
        Vz.factorize_matrix();

        Bx.zero();
        By.zero();
        Bz.zero();

        auto tau = steps.dt;
        auto h = tau * tau / (4 * problem.eps * problem.mu);

        form_matrix(Bx, Ux.basis, [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; });
        form_matrix(By, Uy.basis, [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; });
        form_matrix(Bz, Uz.basis, [h](auto u, auto v) { return u.val * v.val + h * u.dx * v.dx; });

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

        project(now.E1, Ux, Uy, Uz, problem.E1_val_at(0));
        project(now.E2, Ux, Uy, Uz, problem.E2_val_at(0));
        project(now.E3, Ux, Uy, Uz, problem.E3_val_at(0));

        project(now.H1, Ux, Uy, Uz, problem.H1_val_at(0));
        project(now.H2, Ux, Uy, Uz, problem.H2_val_at(0));
        project(now.H3, Ux, Uy, Uz, problem.H3_val_at(0));

        after_step(-1, -steps.dt);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);
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
        compute_rhs(rhs_E1, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[X].val + a * (H[Z].dy - H[Y].dz)) * v.val + b * E[Y].dx * v.dy;
        });
        ads_solve(rhs_E1, buffer, Vx.data(), dim_data{By, By_ctx}, Vz.data());

        auto rhs_E2 = vector_type{shape_E2};
        compute_rhs(rhs_E2, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[Y].val + a * (H[X].dz - H[Z].dx)) * v.val + b * E[Z].dy * v.dz;
        });
        ads_solve(rhs_E2, buffer, Vx.data(), Vy.data(), dim_data{Bz, Bz_ctx});

        auto rhs_E3 = vector_type{shape_E3};
        compute_rhs(rhs_E3, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[Z].val + a * (H[Y].dx - H[X].dy)) * v.val + b * E[X].dz * v.dx;
        });
        ads_solve(rhs_E3, buffer, dim_data{Bx, Bx_ctx}, Vy.data(), Vz.data());

        // First substep - H
        auto rhs_H1 = vector_type{shape_H1};
        compute_rhs(rhs_H1, prev, rhs_E1, rhs_E2, rhs_E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[X].val - c * (E[Z].dy - En[Y].dz)) * v.val;
        });
        ads_solve(rhs_H1, buffer, Vx.data(), Vy.data(), Vz.data());

        auto rhs_H2 = vector_type{shape_H2};
        compute_rhs(rhs_H2, prev, rhs_E1, rhs_E2, rhs_E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[Y].val - c * (E[X].dz - En[Z].dx)) * v.val;
        });
        ads_solve(rhs_H2, buffer, Vx.data(), Vy.data(), Vz.data());

        auto rhs_H3 = vector_type{shape_H3};
        compute_rhs(rhs_H3, prev, rhs_E1, rhs_E2, rhs_E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[Z].val - c * (E[Y].dx - En[X].dy)) * v.val;
        });
        ads_solve(rhs_H3, buffer, Vx.data(), Vy.data(), Vz.data());

        prev.E1 = rhs_E1;
        prev.E2 = rhs_E2;
        prev.E3 = rhs_E3;

        prev.H1 = rhs_H1;
        prev.H2 = rhs_H2;
        prev.H3 = rhs_H3;

        // Second substep - E
        compute_rhs(now.E1, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[X].val + a * (H[Z].dy - H[Y].dz)) * v.val + b * E[Z].dx * v.dz;
        });
        ads_solve(now.E1, buffer, Vx.data(), Vy.data(), dim_data{Bz, Bz_ctx});

        compute_rhs(now.E2, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[Y].val + a * (H[X].dz - H[Z].dx)) * v.val + b * E[X].dy * v.dx;
        });
        ads_solve(now.E2, buffer, dim_data{Bx, Bx_ctx}, Vy.data(), Vz.data());

        compute_rhs(now.E3, prev, prev.E1, prev.E2, prev.E3, Vx, Vy, Vz, [=](auto E, auto, auto H, auto v) {
            return (E[Z].val + a * (H[Y].dx - H[X].dy)) * v.val + b * E[Y].dz * v.dy;
        });
        ads_solve(now.E3, buffer, Vx.data(), dim_data{By, By_ctx}, Vz.data());

        // Second substep - H
        compute_rhs(now.H1, prev, now.E1, now.E2, now.E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[X].val - c * (En[Z].dy - E[Y].dz)) * v.val;
        });
        ads_solve(now.H1, buffer, Vx.data(), Vy.data(), Vz.data());

        compute_rhs(now.H2, prev, now.E1, now.E2, now.E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[Y].val - c * (En[X].dz - E[Z].dx)) * v.val;
        });
        ads_solve(now.H2, buffer, Vx.data(), Vy.data(), Vz.data());

        compute_rhs(now.H3, prev, now.E1, now.E2, now.E3, Vx, Vy, Vz, [=](auto E, auto En, auto H, auto v) {
            return (H[Z].val - c * (En[Y].dx - E[X].dy)) * v.val;
        });
        ads_solve(now.H3, buffer, Vx.data(), Vy.data(), Vz.data());
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

                auto E1 = eval(prev.E1, e, q, Ux, Uy, Uz);
                auto E2 = eval(prev.E2, e, q, Ux, Uy, Uz);
                auto E3 = eval(prev.E3, e, q, Ux, Uy, Uz);
                auto E = vec{E1, E2, E3};

                auto E1_n = eval(E1_new, e, q, Ux, Uy, Uz);
                auto E2_n = eval(E2_new, e, q, Ux, Uy, Uz);
                auto E3_n = eval(E3_new, e, q, Ux, Uy, Uz);
                auto E_n = vec{E1_n, E2_n, E3_n};

                auto H1 = eval(prev.H1, e, q, Ux, Uy, Uz);
                auto H2 = eval(prev.H2, e, q, Ux, Uy, Uz);
                auto H3 = eval(prev.H3, e, q, Ux, Uy, Uz);
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

        auto E1_norm_L2 = normL2(now.E1, Ux, Uy, Uz);
        auto E2_norm_L2 = normL2(now.E2, Ux, Uy, Uz);
        auto E3_norm_L2 = normL2(now.E3, Ux, Uy, Uz);
        auto E_norm_L2 = std::sqrt(E1_norm_L2 * E1_norm_L2 + E2_norm_L2 * E2_norm_L2 + E3_norm_L2 * E3_norm_L2);

        auto E1_err_L2 = errorL2(now.E1, Ux, Uy, Uz, problem.E1_at(tt));
        auto E2_err_L2 = errorL2(now.E2, Ux, Uy, Uz, problem.E2_at(tt));
        auto E3_err_L2 = errorL2(now.E3, Ux, Uy, Uz, problem.E3_at(tt));

        auto rot_E = norm_rot(now.E1, now.E2, now.E3, Ux, Uy, Uz);

        auto H1_norm_L2 = normL2(now.H1, Ux, Uy, Uz);
        auto H2_norm_L2 = normL2(now.H2, Ux, Uy, Uz);
        auto H3_norm_L2 = normL2(now.H3, Ux, Uy, Uz);
        auto H_norm_L2 = std::sqrt(H1_norm_L2 * H1_norm_L2 + H2_norm_L2 * H2_norm_L2 + H3_norm_L2 * H3_norm_L2);

        auto H1_err_L2 = errorL2(now.H1, Ux, Uy, Uz, problem.H1_at(tt));
        auto H2_err_L2 = errorL2(now.H2, Ux, Uy, Uz, problem.H2_at(tt));
        auto H3_err_L2 = errorL2(now.H3, Ux, Uy, Uz, problem.H3_at(tt));

        auto rot_H = norm_rot(now.H1, now.H2, now.H1, Ux, Uy, Uz);

        std::cout << "After step " << i << ", t = " << tt << std::endl;
        std::cout << "  |E|     = " << E_norm_L2 << std::endl;
        std::cout << "  |rot E| = " << rot_E << std::endl;
        std::cout << "  E1 err = " << E1_err_L2 << std::endl;
        std::cout << "  E2 err = " << E2_err_L2 << std::endl;
        std::cout << "  E3 err = " << E3_err_L2 << std::endl;

        std::cout << "  |H| = " << H_norm_L2 << std::endl;
        std::cout << "  |rot H| = " << rot_H << std::endl;
        std::cout << "  H1 err = " << H1_err_L2 << std::endl;
        std::cout << "  H2 err = " << H2_err_L2 << std::endl;
        std::cout << "  H3 err = " << H3_err_L2 << std::endl;
    }

};

}

#endif // ifndef PROBLEMS_MAXWELL_MAXWELL_GALERKIN_HPP_
