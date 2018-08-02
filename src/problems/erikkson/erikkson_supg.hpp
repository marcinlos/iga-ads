#ifndef PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"


#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class erikkson_supg : public simulation_2d {
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

    double C1 = 4, C2 = 2;

    point_type c_diff{{ epsilon, epsilon }};

    // double angle = 0;
    double angle = M_PI / 6;

    double len = 1;

    point_type beta{{ len * cos(angle), len * sin(angle) }};

    mumps::solver solver;

    output_manager<2> output;

public:
    erikkson_supg(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
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
        auto N = Vx.dofs() * Vy.dofs();

        // Identity matrix in the upper left
        for (auto i : dofs(Vx, Vy)) {
            int ii = linear_index(i, Vx, Vy) + 1;
            problem.add(ii, ii, 1);
        }

        for (auto a : internal_dofs(Ux, Uy)) {
            for (auto b : overlapping_dofs(a, Ux, Uy)) {

                double bwu =
                    + c_diff[0] * kron(KUVx, MUVy, a, b) + c_diff[1] * kron(MUVx, KUVy, a, b)
                    + beta[0] * kron(AUVx, MUVy, a, b) + beta[1] * kron(MUVx, AUVy, a, b);

                double supg = 0;
                for (auto e : elements_supporting_dof(a, Ux, Uy)) {
                    auto xrange = Ux.basis.element_ranges[b[0]];
                    auto yrange = Uy.basis.element_ranges[b[1]];

                    if (e[0] < xrange.first || e[0] > xrange.second || e[1] < yrange.first || e[1] > yrange.second) {
                        continue;
                    }
                    double J = jacobian(e);
                    for (auto q : quad_points(Ux, Uy)) {
                        double w = weigth(q);
                        value_type ww = eval_basis(e, q, a, Ux, Uy);
                        value_type uu = eval_basis(e, q, b, Ux, Uy);

                        double hx = 2 * Ux.basis.J[e[0]];
                        double hy = 2 * Uy.basis.J[e[1]];
                        // double h = hx;
                        double h = std::sqrt(hx * hy);

                        double tau = 1 / (C1 * epsilon / (h * h) + C2 / h);

                        double bDw = beta[0] * ww.dx + beta[1] * ww.dy;

                        double lap = laplacian(e, q, b, Ux, Uy);
                        double res = - epsilon * lap + beta[0] * uu.dx + beta[1] * uu.dy;
                        double v = bDw * res;
                        supg += tau * v * w * J;
                    }
                }

                double val = bwu + supg;

                if (val != 0 && ! is_boundary(a, Ux, Uy)) {
                    int i = linear_index(a, Ux, Uy) + 1;
                    int j = linear_index(b, Ux, Uy) + 1;
                    problem.add(N + i, N + j, val);
                }
            }
        }

        // 1's for Dirichlet BC
        for_boundary_dofs(Ux, Uy, [&](index_type dof) {
            int i = linear_index(dof, Ux, Uy) + 1;
            problem.add(N + i, N + i, 1);
        });
    }

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);
        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);

        gram_matrix_1d(MUVx, Ux.basis, Vx.basis);
        gram_matrix_1d(MUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KUVx, Ux.basis, Vx.basis);
        stiffness_matrix_1d(KUVy, Uy.basis, Vy.basis);

        advection_matrix_1d(AUVx, Ux.basis, Vx.basis);
        advection_matrix_1d(AUVy, Uy.basis, Vy.basis);
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

        apply_bc(u);
        output.to_file(u, "out_0.data");
    }

    template <typename RHS>
    void apply_bc(RHS& u_rhs) {
        dirichlet_bc(u_rhs, boundary::left, Ux, Uy, [](double t) { return t < 0.5 ? 1 : 0; });
        dirichlet_bc(u_rhs, boundary::right, Ux, Uy, 0);
        dirichlet_bc(u_rhs, boundary::bottom, Ux, Uy, 1);
        dirichlet_bc(u_rhs, boundary::top, Ux, Uy, 0);
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

        apply_bc(view_out);

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

    void compute_rhs() {
        zero(rhs);
        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto U = vector_type{{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() }};

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double w = weigth(q);
                auto x = point(e, q);

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



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP */
