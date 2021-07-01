// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_MAXWELL_MAXWELL_HEAD_HPP_
#define PROBLEMS_MAXWELL_MAXWELL_HEAD_HPP_

#include <cstdint>

#include "ads/executor/galois.hpp"
#include "ads/form_matrix.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"
#include "problems.hpp"
#include "state.hpp"


namespace ads {

class maxwell_head_problem {
private:
    maxwell_manufactured1 init{1, 1};

    enum material { AIR, TISSUE, BONE };

    static constexpr double eps_vals[] = { 1.0, 45.8, 16.6 };
    static constexpr double mu_vals[]  = { 1.0, 1.0, 1.0 };

    using byte = std::uint8_t;
    using density_data = lin::tensor<byte, 3>;

    density_data density_map;

public:
    using point_type = std::array<double, 3>;

    maxwell_head_problem(const std::string& data_path)
    : density_map{read_density_data(data_path)}
    { }

    auto init_E1() const { return init.E1_val_at(0); }
    auto init_E2() const { return init.E2_val_at(0); }
    auto init_E3() const { return init.E3_val_at(0); }

    auto init_H1() const { return init.H1_val_at(0); }
    auto init_H2() const { return init.H2_val_at(0); }
    auto init_H3() const { return init.H3_val_at(0); }

    double eps(point_type x) const {
        auto n = density(x);
        return eps_vals[as_material(n)];
    }

    double mu(point_type x) const {
        int n = density(x);
        return mu_vals[as_material(n)];
    }

private:

    static density_data read_density_data(const std::string& path) {
        std::ifstream input{path};

        int nx, ny, nz;
        input >> nz >> ny >> nx;
        auto data = density_data{{nx, ny, nz}};

        for (int i = 0; i < nx * ny * nz; ++ i) {
            int ix, iy, iz, val;
            input >> iz >> iy >> ix >> val;
            data(ix, iy, iz) = static_cast<byte>(val);
        }

        return data;
    }

    byte density(point_type x) const {
        int ix = to_index(x[0], density_map.size(0));
        int iy = to_index(x[1], density_map.size(1));
        int iz = to_index(x[2], density_map.size(2));

        return density_map(ix, iy, iz);
    }

    static int to_index(double t, int n) {
        return std::min(static_cast<int>(t * n), n - 1);
    }

    material as_material(byte n) const {
        if (n <= 1) return AIR;
        else if (n <= 240) return TISSUE;
        else return BONE;
    }
};


class maxwell_head : public simulation_3d {
private:
    using Base = simulation_3d;
    using Problem = maxwell_head_problem;

    galois_executor executor{4};

    dimension UE1x, UE1y, UE1z;
    dimension UE2x, UE2y, UE2z;
    dimension UE3x, UE3y, UE3z;

    dimension UH1x, UH1y, UH1z;
    dimension UH2x, UH2y, UH2z;
    dimension UH3x, UH3y, UH3z;

    dimension Vx, Vy, Vz;

    mumps::problem E1_1, E2_1, E3_1;
    mumps::problem E1_2, E2_2, E3_2;

    state prev, half, now;

    Problem problem{"mri.dat"};
    mumps::solver solver;

    ads::output_manager<3> output;

public:
    maxwell_head(const config_3d& config)
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
    , E1_1{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , E2_1{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , E3_1{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , E1_2{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , E2_2{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , E3_2{nullptr, Vx.dofs() * Vy.dofs() * Vz.dofs()}
    , prev{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
    , half{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
    , now{{Vx.dofs(), Vy.dofs(), Vz.dofs()}}
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

    template <typename BC>
    void matrix(mumps::problem& A, double tau, double cx, double cy, double cz, BC&& bc) {
        using shape = std::array<std::size_t, 6>;
        auto shape_loc = shape{
            Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element(), Vz.basis.dofs_per_element(),
            Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element(), Vz.basis.dofs_per_element()
        };

        for (auto e : elements(Vx, Vy, Vz)) {
            auto J = jacobian(e, Vx, Vy, Vz);
            using buffer_type = lin::tensor<double, 6>;
            auto loc = buffer_type{shape_loc};

            for (auto q : quad_points(Vx, Vy, Vz)) {
                auto W = weigth(q, Vx, Vy, Vz);
                auto x = point(e, q, Vx, Vy, Vz);
                double eps = problem.eps(x);
                double mu  = problem.mu(x);

                for (auto i : dofs_on_element(e, Vx, Vy, Vz)) {
                    auto il = dof_global_to_local(e, i, Vx, Vy, Vz);
                    auto v = eval_basis(e, q, i, Vx, Vy, Vz);
                    for (auto j : dofs_on_element(e, Vx, Vy, Vz)) {
                        auto jl = dof_global_to_local(e, j, Vx, Vy, Vz);
                        auto u = eval_basis(e, q, j, Vx, Vy, Vz);

                        auto M = u.val * v.val;
                        auto S = cx * u.dx * v.dx + cy * u.dy * v.dy + cz * u.dz * v.dz;
                        auto form = M + tau * tau / (4 * eps * mu) * S;
                        loc(il[0], il[1], il[2], jl[0], jl[1], jl[2]) += form * W * J;
                    }
                }
            }
            // Update global matrix
            for (auto i : dofs_on_element(e, Vx, Vy, Vz)) {
                if (bc(i)) continue;
                auto il = dof_global_to_local(e, i, Vx, Vy, Vz);
                int ii = linear_index(i, Vx, Vy, Vz) + 1;
                for (auto j : dofs_on_element(e, Vx, Vy, Vz)) {
                    if (bc(j)) continue;
                    auto jl = dof_global_to_local(e, j, Vx, Vy, Vz);
                    int jj = linear_index(j, Vx, Vy, Vz) + 1;

                    auto val = loc(il[0], il[1], il[2], jl[0], jl[1], jl[2]);
                    A.add(ii, jj, val);
                }
            }
        }
        // Dirichlet BC
        for (auto i : dofs(Vx, Vy, Vz)) {
            if (bc(i)) {
                int ii = linear_index(i, Vx, Vy, Vz) + 1;
                A.add(ii, ii, 1.0);
            }
        }
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


        auto tau = steps.dt;

        // E x n = 0
        matrix(E1_1, tau, 0, 1, 0, [this](auto i) { return is_boundary(i[1], Vy) || is_boundary(i[2], Vz); });
        matrix(E2_1, tau, 0, 0, 1, [this](auto i) { return is_boundary(i[0], Vx) || is_boundary(i[2], Vz); });
        matrix(E3_1, tau, 1, 0, 0, [this](auto i) { return is_boundary(i[0], Vx) || is_boundary(i[1], Vy); });

        matrix(E1_2, tau, 0, 0, 1, [this](auto i) { return is_boundary(i[1], Vy) || is_boundary(i[2], Vz); });
        matrix(E2_2, tau, 1, 0, 0, [this](auto i) { return is_boundary(i[0], Vx) || is_boundary(i[2], Vz); });
        matrix(E3_2, tau, 0, 1, 0, [this](auto i) { return is_boundary(i[0], Vx) || is_boundary(i[1], Vy); });
    }

    void before() override {
        prepare_matrices();

        auto project = [&](auto& rhs, auto& x, auto& y, auto& z, auto fun) {
            auto f = [&](double x, double y, double z) { return fun({x, y, z}); };
            auto buffer = vector_type{{ x.dofs(), y.dofs(), z.dofs() }};
            compute_projection(rhs, x.basis, y.basis, z.basis, f);
            ads_solve(rhs, buffer, x.data(), y.data(), z.data());
        };

        // project(now.H1, UH1x, UH1y, UH1z, [this](point_type x) { return problem.eps(x); });
        // output.to_file("eps.vti", output.evaluate(now.H1));

        project(now.E1, UE1x, UE1y, UE1z, problem.init_E1());
        project(now.E2, UE2x, UE2y, UE2z, problem.init_E2());
        project(now.E3, UE3x, UE3y, UE3z, problem.init_E3());

        project(now.H1, UH1x, UH1y, UH1z, problem.init_H1());
        project(now.H2, UH2x, UH2y, UH2z, problem.init_H2());
        project(now.H3, UH3x, UH3y, UH3z, problem.init_H3());

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
        // auto eps = problem.eps;
        // auto mu = problem.mu;
        auto a = [this,tau](auto x) { return tau / (2 * problem.eps(x)); };
        auto b = [this,tau](auto x) { return tau * tau / (4 * problem.eps(x)); };
        auto c = [this,tau](auto x) { return tau / (2 * problem.mu(x)); };

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
        compute_rhs(rhs_E1, prev, prev.E1, prev.E2, prev.E3, UE1x, UE1y, UE1z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[X].val + a(x) * (H[Z].dy - H[Y].dz)) * v.val + b(x) * E[Y].dx * v.dy;
        });
        zero_sides("yz", rhs_E1, UE1x, UE1y, UE1z);
        // ads_solve(rhs_E1, buffer, UE1x.data(), dim_data{By, By_ctx}, UE1z.data());
        E1_1.rhs(rhs_E1.data());
        solver.solve(E1_1);

        auto rhs_E2 = vector_type{shape_E2};
        compute_rhs(rhs_E2, prev, prev.E1, prev.E2, prev.E3, UE2x, UE2y, UE2z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Y].val + a(x) * (H[X].dz - H[Z].dx)) * v.val + b(x) * E[Z].dy * v.dz;
        });
        zero_sides("xz", rhs_E2, UE2x, UE2y, UE2z);
        // ads_solve(rhs_E2, buffer, UE2x.data(), UE2y.data(), dim_data{Bz, Bz_ctx});
        E2_1.rhs(rhs_E2.data());
        solver.solve(E2_1);

        auto rhs_E3 = vector_type{shape_E3};
        compute_rhs(rhs_E3, prev, prev.E1, prev.E2, prev.E3, UE3x, UE3y, UE3z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Z].val + a(x) * (H[Y].dx - H[X].dy)) * v.val + b(x) * E[X].dz * v.dx;
        });
        zero_sides("xy", rhs_E3, UE3x, UE3y, UE3z);
        // ads_solve(rhs_E3, buffer, dim_data{Bx, Bx_ctx}, UE3y.data(), UE3z.data());
        E3_1.rhs(rhs_E3.data());
        solver.solve(E3_1);

        // First substep - H
        auto rhs_H1 = vector_type{shape_H1};
        compute_rhs(rhs_H1, prev, rhs_E1, rhs_E2, rhs_E3, UH1x, UH1y, UH1z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[X].val - c(x) * (E[Z].dy - En[Y].dz)) * v.val;
        });
        zero_sides("x", rhs_H1, UH1x, UH1y, UH1z);
        ads_solve(rhs_H1, buffer, UH1x.data(), UH1y.data(), UH1z.data());

        auto rhs_H2 = vector_type{shape_H2};
        compute_rhs(rhs_H2, prev, rhs_E1, rhs_E2, rhs_E3, UH2x, UH2y, UH2z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Y].val - c(x) * (E[X].dz - En[Z].dx)) * v.val;
        });
        zero_sides("y", rhs_H2, UH2x, UH2y, UH2z);
        ads_solve(rhs_H2, buffer, UH2x.data(), UH2y.data(), UH2z.data());

        auto rhs_H3 = vector_type{shape_H3};
        compute_rhs(rhs_H3, prev, rhs_E1, rhs_E2, rhs_E3, UH3x, UH3y, UH3z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Z].val - c(x) * (E[Y].dx - En[X].dy)) * v.val;
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
        compute_rhs(now.E1, prev, prev.E1, prev.E2, prev.E3, UE1x, UE1y, UE1z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[X].val + a(x) * (H[Z].dy - H[Y].dz)) * v.val + b(x) * E[Z].dx * v.dz;
        });
        zero_sides("yz", now.E1, UE1x, UE1y, UE1z);
        // ads_solve(now.E1, buffer, UE1x.data(), UE1y.data(), dim_data{Bz, Bz_ctx});
        E1_2.rhs(now.E1.data());
        solver.solve(E1_2);

        compute_rhs(now.E2, prev, prev.E1, prev.E2, prev.E3, UE2x, UE2y, UE2z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Y].val + a(x) * (H[X].dz - H[Z].dx)) * v.val + b(x) * E[X].dy * v.dx;
        });
        zero_sides("xz", now.E2, UE2x, UE2y, UE2z);
        // ads_solve(now.E2, buffer, dim_data{Bx, Bx_ctx}, UE2y.data(), UE2z.data());
        E2_2.rhs(now.E2.data());
        solver.solve(E2_2);

        compute_rhs(now.E3, prev, prev.E1, prev.E2, prev.E3, UE3x, UE3y, UE3z, [=](auto E, auto, auto H, auto v, auto x) {
            return (E[Z].val + a(x) * (H[Y].dx - H[X].dy)) * v.val + b(x) * E[Y].dz * v.dy;
        });
        zero_sides("xy", now.E3, UE3x, UE3y, UE3z);
        // ads_solve(now.E3, buffer, UE3x.data(), dim_data{By, By_ctx}, UE3z.data());
        E3_2.rhs(now.E3.data());
        solver.solve(E3_2);

        // Second substep - H
        compute_rhs(now.H1, prev, now.E1, now.E2, now.E3, UH1x, UH1y, UH1z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[X].val - c(x) * (En[Z].dy - E[Y].dz)) * v.val;
        });
        zero_sides("x", now.H1, UH1x, UH1y, UH1z);
        ads_solve(now.H1, buffer, UH1x.data(), UH1y.data(), UH1z.data());

        compute_rhs(now.H2, prev, now.E1, now.E2, now.E3, UH2x, UH2y, UH2z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Y].val - c(x) * (En[X].dz - E[Z].dx)) * v.val;
        });
        zero_sides("y", now.H2, UH2x, UH2y, UH2z);
        ads_solve(now.H2, buffer, UH2x.data(), UH2y.data(), UH2z.data());

        compute_rhs(now.H3, prev, now.E1, now.E2, now.E3, UH3x, UH3y, UH3z, [=](auto E, auto En, auto H, auto v, auto x) {
            return (H[Z].val - c(x) * (En[Y].dx - E[X].dy)) * v.val;
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
                auto x = point(e, q, Vx, Vy, Vz);

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

                    auto val = form(E, E_n, H, v, x);
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

        auto E1_norm_H1 = normH1(now.E1, UE1x, UE1y, UE1z);
        auto E2_norm_H1 = normH1(now.E2, UE2x, UE2y, UE2z);
        auto E3_norm_H1 = normH1(now.E3, UE3x, UE3y, UE3z);
        auto E_norm_H1 = std::sqrt(E1_norm_H1 * E1_norm_H1 + E2_norm_H1 * E2_norm_H1 + E3_norm_H1 * E3_norm_H1);

        auto rot_E = norm_rot(now.E1, now.E2, now.E3, Vx, Vy, Vz);
        auto div_E = norm_div(now.E1, now.E2, now.E3, Vx, Vy, Vz);

        auto H1_norm_L2 = normL2(now.H1, UH1x, UH1y, UH1z);
        auto H2_norm_L2 = normL2(now.H2, UH2x, UH2y, UH2z);
        auto H3_norm_L2 = normL2(now.H3, UH3x, UH3y, UH3z);
        auto H_norm_L2 = std::sqrt(H1_norm_L2 * H1_norm_L2 + H2_norm_L2 * H2_norm_L2 + H3_norm_L2 * H3_norm_L2);

        auto H1_norm_H1 = normH1(now.H1, UH1x, UH1y, UH1z);
        auto H2_norm_H1 = normH1(now.H2, UH2x, UH2y, UH2z);
        auto H3_norm_H1 = normH1(now.H3, UH3x, UH3y, UH3z);
        auto H_norm_H1 = std::sqrt(H1_norm_H1 * H1_norm_H1 + H2_norm_H1 * H2_norm_H1 + H3_norm_H1 * H3_norm_H1);

        auto rot_H = norm_rot(now.H1, now.H2, now.H3, Vx, Vy, Vz);
        auto div_H = norm_div(now.H1, now.H2, now.H3, Vx, Vy, Vz);

        std::cout << "After step " << i << ", t = " << tt << std::endl;
        std::cout << "  |E|     = " << E_norm_L2 << "  " << E_norm_H1 << std::endl;
        std::cout << "  |rot E| = " << rot_E << std::endl;
        std::cout << "  |div E| = " << div_E << std::endl;

        std::cout << "  |H| = " << H_norm_L2 << "  " << H_norm_H1 << std::endl;
        std::cout << "  |rot H| = " << rot_H << std::endl;
        std::cout << "  |div H| = " << div_H << std::endl;
    }

};

}

#endif // ifndef PROBLEMS_MAXWELL_MAXWELL_HEAD_HPP_
