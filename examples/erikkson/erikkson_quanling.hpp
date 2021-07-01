// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef PROBLEMS_ERIKKSON_ERIKKSON_QUANLING_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_QUANLING_HPP_

#include "ads/executor/galois.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"
#include "ads/lin/tensor/view.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"


namespace ads {

class erikkson_quanling : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;

    galois_executor executor{8};

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
    };

    vector_type u;
    residuum r;
    vector_type u_buffer;
    std::vector<double> full_rhs;

    int save_every = 1;

    double minh = 1 / 15.;//1e-7;
    double minh2 = minh * minh;

    double rho = 0.9;
    double gamma = 10;


    double pecelet = 1e6;
    double epsilon = 1 / pecelet;

    point_type c_diff{{ epsilon, epsilon }};

    double angle = 0;
    double len = 1;

    point_type beta{{ len * cos(angle), len * sin(angle) }};

    mumps::solver solver;

    output_manager<2> output;

public:
    erikkson_quanling(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
    : Base{std::move(test_x), std::move(test_y), steps}
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , Vx{ x }
    , Vy{ y }
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , MUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KUx{Ux.p, Ux.p, Ux.dofs(), Ux.dofs(), 0}
    , KUy{Uy.p, Uy.p, Uy.dofs(), Uy.dofs(), 0}
    , MUVx{ Vx.dofs(), Ux.dofs() }
    , MUVy{ Vy.dofs(), Uy.dofs() }
    , KUVx{ Vx.dofs(), Ux.dofs() }
    , KUVy{ Vy.dofs(), Uy.dofs() }
    , AUVx{ Vx.dofs(), Ux.dofs() }
    , AUVy{ Vy.dofs(), Uy.dofs() }
    , MUUx{ Ux.dofs(), Ux.dofs() }
    , MUUy{ Uy.dofs(), Uy.dofs() }
    , KUUx{ Ux.dofs(), Ux.dofs() }
    , KUUy{ Uy.dofs(), Uy.dofs() }
    , AUUx{ Ux.dofs(), Ux.dofs() }
    , AUUy{ Uy.dofs(), Uy.dofs() }
    , u{{ Ux.dofs(), Uy.dofs() }}
    , r{ {{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , output{ Ux.B, Uy.B, 500 }
    { }

private:

    struct matrix_set {
        using band_matrix_ref = lin::band_matrix&;
        using dense_matrix_ref = lin::dense_matrix&;

        band_matrix_ref MVx, MVy, KVx, KVy;
        dense_matrix_ref MUVx, MUVy, KUVx, KUVy, AUVx, AUVy;
    };

    void assemble_problem(mumps::problem& problem, const dimension& Vx, const dimension& Vy, const matrix_set& M) {
        using std::min;
        using std::max;
        auto N = Vx.dofs() * Vy.dofs();
        vector_type v{{ Vx.dofs(), Vy.dofs() }};

        // Gram matrix - upper left
        for (auto ix = 1; ix < Vx.dofs() - 1; ++ ix) {
            for (auto iy = 1; iy < Vy.dofs() - 1; ++ iy) {
                int i = &v(ix, iy) - &v(0, 0) + 1;
                for (auto jx = max(1, ix - Vx.B.degree); jx < min(Vx.dofs() - 1, ix + Vx.B.degree + 1); ++ jx) {
                    for (auto jy = max(1, iy - Vy.B.degree); jy < min(Vy.dofs() - 1, iy + Vy.B.degree + 1); ++ jy) {
                        int j = &v(jx, jy) - &v(0, 0) + 1;
                        double val = M.MVx(ix, jx) * M.MVy(iy, jy);
                        val += minh2 * M.KVx(ix, jx) * M.MVy(iy, jy);
                        val += minh2 * M.MVx(ix, jx) * M.KVy(iy, jy);
                        // val += minh2 * minh2 * M.KVx(ix, jx) * M.KVy(iy, jy);

                        problem.add(i, j, val);
                    }
                }
            }
        }

        // Dirichlet BC - upper left
        for (auto jx = 0; jx < Vx.dofs(); ++ jx) {
            int j = &v(jx, 0) - &v(0, 0) + 1;
            problem.add(j, j, 1);
            j = &v(jx, Vy.dofs() - 1) - &v(0, 0) + 1;
            problem.add(j, j, 1);

        }
        for (auto jy = 1; jy < Vy.dofs() - 1; ++ jy) {
            int j = &v(0, jy) - &v(0, 0) + 1;
            problem.add(j, j, 1);
            j = &v(Vx.dofs() - 1, jy) - &v(0, 0) + 1;
            problem.add(j, j, 1);
        }

        // B, B^T
        // for (auto ix = 0; ix < Vx.dofs(); ++ ix) {
        //     for (auto jx = 0; jx < Ux.dofs(); ++ jx) {
        //         for (auto iy = 0; iy < Vy.dofs(); ++ iy) {
        //             for (auto jy = 0; jy < Uy.dofs(); ++ jy) {
        //                 int i = &v(ix, iy) - &v(0, 0) + 1;
        //                 int j = &u(jx, jy) - &u(0, 0) + 1;
        //                 double val = 0;
        //                 val += c_diff[0] * M.KUVx(ix, jx) * M.MUVy(iy, jy) + beta[0] * M.AUVx(ix, jx) * M.MUVy(iy, jy);
        //                 val += c_diff[1] * M.MUVx(ix, jx) * M.KUVy(iy, jy) + beta[1] * M.MUVx(ix, jx) * M.AUVy(iy, jy);
        //                 if (val != 0) {
        //                     if (ix != 0 && ix != Vx.dofs() - 1 && iy != 0 && iy != Vy.dofs() - 1
        //                         && jx != 0 && jx != Ux.dofs() - 1 && jy != 0 && jy != Uy.dofs() - 1
        //                         ) {
        //                         problem.add(i, N + j, val);
        //                     }
        //                     if (jx != 0 && jx != Ux.dofs() - 1 && jy != 0 && jy != Uy.dofs() - 1 &&
        //                         ix != 0 && ix != Vx.dofs() - 1 && iy != 0 && iy != Vy.dofs() - 1) {
        //                         problem.add(N + j, i, val);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // G - lower right
        double alpha = 1;
        double c = gamma * std::pow(rho, alpha);

        for (auto ix = 1; ix < Ux.dofs() - 1; ++ ix) {
            for (auto iy = 1; iy < Uy.dofs() - 1; ++ iy) {
                int i = &u(ix, iy) - &u(0, 0) + 1;
                for (auto jx = max(1, ix - Ux.B.degree); jx < min(Ux.dofs() - 1, ix + Ux.B.degree + 1); ++ jx) {
                    for (auto jy = max(1, iy - Uy.B.degree); jy < min(Uy.dofs() - 1, iy + Uy.B.degree + 1); ++ jy) {
                        int j = &u(jx, jy) - &u(0, 0) + 1;
                        double M = MUx(ix, jx) * MUy(iy, jy);
                        double K = MUx(ix, jx) * KUy(iy, jy) + KUx(ix, jx) * MUy(iy, jy);
                        double val = c * (M + minh2 * K);
                        problem.add(N + i, N + j, val);
                    }
                }
            }
        }

        // Dirichlet BC - lower right
        for (auto jx = 0; jx < Ux.dofs(); ++ jx) {
            int j = &u(jx, 0) - &u(0, 0) + 1;
            problem.add(N + j, N + j, 1);
            j = &u(jx, Uy.dofs() - 1) - &u(0, 0) + 1;
            problem.add(N + j, N + j, 1);

        }
        for (auto jy = 1; jy < Uy.dofs() - 1; ++ jy) {
            int j = &u(0, jy) - &u(0, 0) + 1;
            problem.add(N + j, N + j, 1);
            j = &u(Ux.dofs() - 1, jy) - &u(0, 0) + 1;
            problem.add(N + j, N + j, 1);
        }
    }

    void mass(lin::band_matrix& M, const basis_data& d) {
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
                        M(ia, ib) += va * vb * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void diffusion(lin::band_matrix& M, const basis_data& d) {
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
                        M(ia, ib) += da * db * d.w[q] * d.J[e];
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
                        auto val = va * vb;
                        M(ia, ib) += val * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void diffusion_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto da = bV.b[e][q][1][a];
                        auto db = bU.b[e][q][1][b];
                        auto val = da * db;
                        M(ia, ib) += val * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void advection_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto va = bV.b[e][q][0][a];
                        auto db = bU.b[e][q][1][b];
                        auto val = va * db;
                        M(ia, ib) += val * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    matrix_set matrices(bool x_refined, bool y_refined) {
        if (x_refined && y_refined) {
            return { MVx, MVy, KVx, KVy, MUVx, MUVy, KUVx, KUVy, AUVx, AUVy };
        } else if (x_refined) {
            return { MVx, MUy, KVx, KUy, MUVx, MUUy, KUVx, KUUy, AUVx, AUUy };
        } else if (y_refined) {
            return { MUx, MVy, KUx, KVy, MUUx, MUVy, KUUx, KUVy, AUUx, AUVy };
        } else {
            return { MUx, MUy, KUx, KUy, MUUx, MUUy, KUUx, KUUy, AUUx, AUUy };
        }
    }

    void prepare_matrices() {

        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);

        gram_matrix_1d(MUx, Ux.basis);
        gram_matrix_1d(MUy, Uy.basis);

        diffusion(KVx, Vx.basis);
        diffusion(KVy, Vy.basis);
        diffusion(KUx, Ux.basis);
        diffusion(KUy, Uy.basis);

        mass_matrix(MUVx, Ux.basis, Vx.basis);
        mass_matrix(MUVy, Uy.basis, Vy.basis);

        diffusion_matrix(KUVx, Ux.basis, Vx.basis);
        diffusion_matrix(KUVy, Uy.basis, Vy.basis);

        advection_matrix(AUVx, Ux.basis, Vx.basis);
        advection_matrix(AUVy, Uy.basis, Vy.basis);

        mass_matrix(MUUx, Ux.basis, Ux.basis);
        mass_matrix(MUUy, Uy.basis, Uy.basis);

        diffusion_matrix(KUUx, Ux.basis, Ux.basis);
        diffusion_matrix(KUUy, Uy.basis, Uy.basis);

        advection_matrix(AUUx, Ux.basis, Ux.basis);
        advection_matrix(AUUy, Uy.basis, Uy.basis);
    }

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();

        zero(r.data);
        zero(u);
        apply_bc(u);
        // u(5, 5) = 1;

        // output.to_file(u, "out_0.data");
    }

    void apply_bc(vector_type& u_rhs) {
        lin::band_matrix MUy_loc{ Uy.p, Uy.p, Uy.dofs() };
        gram_matrix_1d(MUy_loc, Uy.basis);
        lin::solver_ctx ctx_y{ MUy_loc };
        lin::factorize(MUy_loc, ctx_y);

        lin::vector buf_x0{{ Uy.dofs() }};
        compute_projection(buf_x0, Uy.basis, [&](double t) {
            return std::sin(M_PI * t);
        });
        lin::solve_with_factorized(MUy_loc, buf_x0, ctx_y);

        for (auto j = 0; j < Uy.dofs(); ++ j) {
            u_rhs(0, j) = buf_x0(j);
            u_rhs(Ux.dofs() - 1, j) = 0;
        }
        for (auto j = 1; j < Ux.dofs(); ++ j) {
            u_rhs(j, 0) = 0;
            u_rhs(j, Uy.dofs() - 1) = 0;
        }
    }

    void zero_bc(vector_view& r_rhs, vector_view& u_rhs) {
        for (auto j = 0; j < Uy.dofs(); ++ j) {
            u_rhs(0, j) = 0;
            u_rhs(Ux.dofs() - 1, j) = 0;
        }
        for (auto j = 1; j < Ux.dofs(); ++ j) {
            u_rhs(j, 0) = 0;
            u_rhs(j, Uy.dofs() - 1) = 0;
        }

        for (auto j = 0; j < Vy.dofs(); ++ j) {
            r_rhs(0, j) = 0;
            r_rhs(Vx.dofs() - 1, j) = 0;
        }
        for (auto j = 1; j < Vx.dofs(); ++ j) {
            r_rhs(j, 0) = 0;
            r_rhs(j, Vy.dofs() - 1) = 0;
        }
    }

    void add_solution(const vector_view& u_rhs, const vector_view& r_rhs, const dimension& Vx, const dimension& Vy) {
        for (auto i = 0; i < Ux.dofs(); ++ i) {
            for (auto j = 0; j < Uy.dofs(); ++ j) {
                u(i, j) += u_rhs(i, j);
            }
        }
        // r = residuum{ {{Vx.dofs(), Vy.dofs()}}, &Vx, &Vy };
        for (auto i = 0; i < Vx.dofs(); ++ i) {
            for (auto j = 0; j < Vy.dofs(); ++ j) {
                r.data(i, j) += r_rhs(i, j);
            }
        }
    }

    double norm(const vector_view& u) {
        double norm = 0;
        for (int i = 0; i < Ux.dofs(); ++ i) {
            for (int j = 0; j < Uy.dofs(); ++ j) {
                auto uu = u(i, j);
                norm += uu * uu;
            }
        }
        norm /= (Ux.dofs() * Uy.dofs());
        return std::sqrt(norm);
    }

    double substep() {
        vector_view r_rhs{full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view u_rhs{full_rhs.data() + r_rhs.size(), {Ux.dofs(), Uy.dofs()}};

        std::fill(begin(full_rhs), end(full_rhs), 0);
        compute_rhs(Vx, Vy, r_rhs, u_rhs);
        zero_bc(r_rhs, u_rhs);

        int size = Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs();
        mumps::problem problem(full_rhs.data(), size);
        assemble_problem(problem, Vx, Vy, matrices(true, true));
        solver.solve(problem);

        add_solution(u_rhs, r_rhs, Vx, Vy);

        return norm(u_rhs);
    }

    void step(int iter, double /*t*/) override {
        zero(r.data);
        apply_bc(u);

        std::cout << "Step " << (iter + 1) << std::endl;
        for (int i = 0; ; ++ i) {
            auto norm = substep();
            std::cout << "  substep " << (i + 1) << ": |eta| = " << norm << std::endl;
            if (norm < 1e-8) {
                break;
            }
        }
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % save_every == 0) {
            std::cout << "Step " << (iter + 1) << " : " << errorL2() << " " << errorH1() << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data");
        std::cout << "{ 'L2': '" << errorL2() << "', 'H1': '" << errorH1() << "'}" << std::endl;

        std::ofstream sol("solution.data");
        for (int i = 0; i < Ux.dofs(); ++ i) {
            for (int j = 0; j < Uy.dofs(); ++ j) {
                sol << i << " " << j << " " << u(i, j) << std::endl;
            }
        }
    }

    void plot_middle(const char* filename) {
        std::ofstream out{filename};
        bspline::eval_ctx ctx_x{ Ux.B.degree }, ctx_y{ Uy.B.degree };

        auto print = [&](double xx) {
            auto val = bspline::eval(xx, 0.5, u, Ux.B, Uy.B, ctx_x, ctx_y);
            out << std::setprecision(16) << xx << " " << val << std::endl;
        };

        print(0);
        auto N = Ux.basis.quad_order;
        for (auto e : Ux.element_indices()) {
            std::vector<double> qs(Ux.basis.x[e], Ux.basis.x[e] + N);
            std::sort(begin(qs), end(qs));
            for (auto xx : qs) {
                print(xx);
            }
        }
        print(1);
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

    void update_global_rhs(vector_view& global, const vector_type& local, index_type e,
                           const dimension& x, const dimension& y) const {
        for (auto a : dofs_on_element(e, x, y)) {
            auto loc = dof_global_to_local(e, a, x, y);
            global(a[0], a[1]) += local(loc[0], loc[1]);
        }
    }


    void compute_rhs(const dimension& Vx, const dimension& Vy, vector_view& r_rhs, vector_view& u_rhs) {
        double c = gamma * std::pow(rho, 1);

        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto R = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};
            auto U = vector_type{{ Ux.basis.dofs_per_element(), Uy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weight(q);
                double WJ = W * J;
                // auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);
                value_type rr = eval(r.data, e, q, *r.Vx, *r.Vy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    double lv = 0;//* v.val;
                    double val = lv;
                    // Bu
                    val += c_diff[0] * uu.dx * v.dx + beta[0] * uu.dx * v.val;
                    val += c_diff[1] * uu.dy * v.dy + beta[1] * uu.dy * v.val;
                    // -Aw
                    val -= (rr.val * v.val + minh2 * (rr.dx * v.dx + rr.dy * v.dy));

                    R(aa[0], aa[1]) += val * WJ;
                }
                for (auto a : dofs_on_element(e, Ux, Uy)) {
                    auto aa = dof_global_to_local(e, a, Ux, Uy);
                    value_type w = eval_basis(e, q, a, Ux, Uy);
                    double val = 0;

                    // -B'w
                    val -= (c_diff[0] * w.dx * rr.dx + beta[0] * w.dx * rr.val);
                    val -= (c_diff[1] * w.dy * rr.dy + beta[1] * w.dy * rr.val);

                    U(aa[0], aa[1]) += val * WJ;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(r_rhs, R, e, Vx, Vy);
                update_global_rhs(u_rhs, U, e, Ux, Uy);
            });
        });
    }

    value_type exact(double x, double y, double eps) const {
        using std::exp;
        auto lambda = M_PI * eps;
        auto del = std::sqrt(1 + 4 * lambda * lambda);
        auto r1 = (1 + del) / (2 * eps);
        auto r2 = (1 - del) / (2 * eps);

        auto norm = exp(-r1) - exp(-r2);
        auto alpha = (exp(r1 * (x - 1)) - exp(r2 * (x - 1))) / norm;
        double val = alpha * std::sin(M_PI * y);

        return {
            val,
            (r1 * exp(r1 * (x - 1)) - r2 * exp(r2 * (x - 1))) / norm * std::sin(M_PI * y),
            alpha * M_PI * std::cos(M_PI * y)
        };
    }

    double errorL2() const {
        double error = 0;

        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e);
            for (auto q : quad_points(Ux, Ux)) {
                double w = weight(q);
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - exact(x[0], x[1], epsilon);
                error += d.val * d.val * w * J;
            }
        }
        return std::sqrt(error);
    }

    double errorH1() const {
        double error = 0;
        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e);
            for (auto q : quad_points(Ux, Ux)) {
                double w = weight(q);
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                auto d = uu - exact(x[0], x[1], epsilon);
                error += (d.val * d.val + d.dx * d.dx + d.dy * d.dy) * w * J;
            }
        }
        return std::sqrt(error);
    }
};

}



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_QUANLING<_HPP */
