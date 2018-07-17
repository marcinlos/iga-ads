#ifndef PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"


#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class erikkson_mumps : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;


    // New stuff
    lin::band_matrix MVx, MVy;
    lin::band_matrix KVx, KVy;

    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;


    vector_type u, rhs;
    vector_type u_buffer;
    std::vector<double> full_rhs;


    int save_every = 1;

    double peclet = 1e6;
    double epsilon = 1 / peclet;

    point_type c_diff{{ epsilon, epsilon }};

    double angle = 0;
    double len = 1;

    point_type beta{{ len * cos(angle), len * sin(angle) }};

    mumps::solver solver;

    output_manager<2> output;

public:
    erikkson_mumps(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
    : Base{std::move(test_x), std::move(test_y), steps}
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , Vx{ x }
    , Vy{ y }
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUVx{ Vx.dofs(), Ux.dofs() }
    , MUVy{ Vy.dofs(), Uy.dofs() }
    , KUVx{ Vx.dofs(), Ux.dofs() }
    , KUVy{ Vy.dofs(), Uy.dofs() }
    , AUVx{ Vx.dofs(), Ux.dofs() }
    , AUVy{ Vy.dofs(), Uy.dofs() }
    , u{{ Ux.dofs(), Uy.dofs() }}
    , rhs{{ Vx.dofs(), Vy.dofs() }}
    , u_buffer{{ Ux.dofs(), Uy.dofs() }}
    , full_rhs(Vx.dofs() * Vy.dofs() + Ux.dofs() * Uy.dofs())
    , output{ Ux.B, Uy.B, 500 }
    { }

private:

    void assemble_problem(mumps::problem& problem, double dt) {
        using std::min;
        using std::max;
        auto N = Vx.dofs() * Vy.dofs();
        vector_type v{{ Vx.dofs(), Vy.dofs() }};

        // for (auto ix = 0; ix < Vx.dofs(); ++ ix) {
        //     for (auto jx = 0; jx < Vx.dofs(); ++ jx) {
        //         for (auto iy = 0; iy < Vy.dofs(); ++ iy) {
        //             for (auto jy = 0; jy < Vy.dofs(); ++ jy) {
        //                 int i = &v(ix, iy) - &v(0, 0) + 1;
        //                 int j = &v(jx, jy) - &v(0, 0) + 1;
        //                 double val = MVx(ix, jx) * MVy(iy, jy) + KVx(ix, jx) * MVy(iy, jy) + MVx(ix, jx) * KVy(iy, jy);
        //                 problem.add(i, j, val);
        //             }
        //         }
        //     }
        // }
        for (auto ix = 1; ix < Vx.dofs() - 1; ++ ix) {
            for (auto iy = 1; iy < Vy.dofs() - 1; ++ iy) {
                int i = &v(ix, iy) - &v(0, 0) + 1;
                for (auto jx = max(1, ix - Vx.B.degree); jx < min(Vx.dofs() - 1, ix + Vx.B.degree + 1); ++ jx) {
                    for (auto jy = max(1, iy - Vy.B.degree); jy < min(Vy.dofs() - 1, iy + Vy.B.degree + 1); ++ jy) {
                        int j = &v(jx, jy) - &v(0, 0) + 1;
                        double val = MVx(ix, jx) * MVy(iy, jy)  + KVx(ix, jx) * MVy(iy, jy) + MVx(ix, jx) * KVy(iy, jy);
                        // std::cout << val << " ";
                        problem.add(i, j, val);
                    }
                }
                // std::cout << std::endl;
            }
        }


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

        for (auto ix = 0; ix < Vx.dofs(); ++ ix) {
            for (auto jx = 0; jx < Ux.dofs(); ++ jx) {
                for (auto iy = 0; iy < Vy.dofs(); ++ iy) {
                    for (auto jy = 0; jy < Uy.dofs(); ++ jy) {
                        int i = &v(ix, iy) - &v(0, 0) + 1;
                        int j = &u(jx, jy) - &u(0, 0) + 1;
                        double val =
                            // MUVx(ix, jx) * MUVy(iy, jy)
                            + (c_diff[0] * KUVx(ix, jx) * MUVy(iy, jy) + c_diff[1] * MUVx(ix, jx) * KUVy(iy, jy))
                            - (beta[0] * AUVx(ix, jx) * MUVy(iy, jy) + beta[1] * MUVx(ix, jx) * AUVy(iy, jy));
                        if (val != 0) {
                            if (ix != 0 && ix != Vx.dofs() - 1 && iy != 0 && iy != Vy.dofs() - 1) {
                                problem.add(i, N + j, -val);
                            }
                            if (jx != 0 && jx != Ux.dofs() - 1 && jy != 0 && jy != Uy.dofs() - 1 &&
                                ix != 0 && ix != Vx.dofs() - 1 && iy != 0 && iy != Vy.dofs() - 1) {
                                problem.add(N + j, i, val);
                            }
                        }
                    }
                }
            }
        }
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
                        auto da = bV.b[e][q][1][a];
                        auto vb = bU.b[e][q][0][b];
                        auto val = da * vb;
                        M(ia, ib) += val * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void prepare_matrices() {

        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);
        diffusion(KVx, Vx.basis);
        diffusion(KVy, Vy.basis);

        mass_matrix(MUVx, Ux.basis, Vx.basis);
        mass_matrix(MUVy, Uy.basis, Vy.basis);

        diffusion_matrix(KUVx, Ux.basis, Vx.basis);
        diffusion_matrix(KUVy, Uy.basis, Vy.basis);

        advection_matrix(AUVx, Ux.basis, Vx.basis);
        advection_matrix(AUVy, Uy.basis, Vy.basis);
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min( (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
        // return 0;
    };

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();


        // auto init = [this](double x, double y) { return init_state(x, y); };

        // compute_projection(u, Ux.basis, Uy.basis, init);
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        // zero(u);
        // u(Ux.dofs() / 2, Uy.dofs() / 2) = 1;

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
            u(0, j) = buf_x0(j);
            u(Ux.dofs() - 1, j) = 0;
        }
        for (auto j = 1; j < Ux.dofs(); ++ j) {
            u(j, 0) = 0;
            u(j, Uy.dofs() - 1) = 0;
        }

        output.to_file(u, "out_0.data");
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs();
        // zero(rhs);

        std::fill(begin(full_rhs), end(full_rhs), 0);
        vector_view view_in{ full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view view_out{ full_rhs.data() + view_in.size(), {Ux.dofs(), Uy.dofs()}};


        for (auto i = 0; i < Vx.dofs(); ++ i) {
            for (auto j = 0; j < Vy.dofs(); ++ j) {
                view_in(i, j) = - rhs(i, j);
            }
        }

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
            view_out(0, j) = buf_x0(j);
            view_out(Ux.dofs() - 1, j) = 0;
        }
        for (auto j = 1; j < Ux.dofs(); ++ j) {
            view_out(j, 0) = 0;
            view_out(j, Uy.dofs() - 1) = 0;
        }

        mumps::problem problem(full_rhs.data(), full_rhs.size());
        assemble_problem(problem, steps.dt);
        solver.solve(problem, "igrm");

        for (auto i = 0; i < Ux.dofs(); ++ i) {
            for (auto j = 0; j < Uy.dofs(); ++ j) {
                u(i, j) = view_out(i, j);
            }
        }
        // for (auto i = 0; i < Vx.dofs(); ++ i) {
        //     for (auto j = 0; j < Vy.dofs(); ++ j) {
        //         std::cout << "res(" << i << ", " << j << ") = " << view_in(i, j) << std::endl;
        //     }
        // }

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

    void update_global_rhs(vector_type& global, const vector_type& local, index_type e,
                           const dimension& x, const dimension& y) const {
        for (auto a : dofs_on_element(e, x, y)) {
            auto loc = dof_global_to_local(e, a, x, y);
            global(a[0], a[1]) += local(loc[0], loc[1]);
        }
    }


    void compute_rhs() {
        zero(rhs);
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto U = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double w = weigth(q);
                auto x = point(e, q);
                value_type uu = eval(u, e, q, Ux, Uy);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);
                    double val = 0*1 * v.val;
                    // double val = init_state(x[0], x[1]) * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e, Vx, Vy);
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
                double w = weigth(q);
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
                double w = weigth(q);
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



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP */
