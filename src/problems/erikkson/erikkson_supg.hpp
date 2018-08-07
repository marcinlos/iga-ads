#ifndef PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP_

#include "problems/erikkson/erikkson_base.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "mumps.hpp"

#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"


namespace ads {

class erikkson_supg : public erikkson_base {
private:
    using Base = erikkson_base;

    galois_executor executor{8};

    // New stuff
    lin::band_matrix Mx, My;
    lin::band_matrix Kx, Ky;
    lin::band_matrix Ax, Ay;

    vector_type u, rhs;

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
    erikkson_supg(dimension trial_x, dimension trial_y, const timesteps_config& steps)
    : Base{std::move(trial_x), std::move(trial_y), steps}
    , Mx{x.p, x.p, x.dofs(), x.dofs(), 0}
    , My{y.p, y.p, y.dofs(), y.dofs(), 0}
    , Kx{x.p, x.p, x.dofs(), x.dofs(), 0}
    , Ky{y.p, y.p, y.dofs(), y.dofs(), 0}
    , Ax{x.p, x.p, x.dofs(), x.dofs(), 0}
    , Ay{y.p, y.p, y.dofs(), y.dofs(), 0}
    , u{{ x.dofs(), y.dofs() }}
    , rhs{{ x.dofs(), y.dofs() }}
    , output{ x.B, y.B, 500 }
    { }

private:

    double elem_diam(index_type e) const {
        double hx = 2 * x.basis.J[e[0]];
        double hy = 2 * y.basis.J[e[1]];
        // double h = hx;
        return std::sqrt(hx * hy);
    }

    double diffusion(double x, double y) const {
        constexpr double eta = 1e6;

        bool left = x < 0.5, right = !left;
        bool bottom = y < 0.5, top = !bottom;

        if ((bottom && left) || (top && right)) {
            return eta;
        } else {
            return 1;
        }
    }

    void assemble_problem(mumps::problem& problem, double dt) {
        for (auto a : internal_dofs(x, y)) {
            for (auto b : overlapping_dofs(a, x, y)) {

                double val = 0;
                for (auto e : elements_supporting_dof(a, x, y)) {
                    if (! supported_in(b, e, x, y)) continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(x, y)) {
                        double w = weigth(q, x, y);
                        auto pt = point(e, q, x, y);
                        value_type ww = eval_basis(e, q, a, x, y);
                        value_type uu = eval_basis(e, q, b, x, y);

                        auto diff = diffusion(pt[0], pt[1]);
                        double bwu = diff * grad_dot(uu, ww) + beta[0] * uu.dx * ww.val + beta[1] * uu.dy * ww.val;

                        double h = elem_diam(e);
                        double tau = 1 / (C1 * epsilon / (h * h) + C2 / h);

                        double bDw = beta[0] * ww.dx + beta[1] * ww.dy;

                        double lap = laplacian(e, q, b, x, y);
                        double res = - epsilon * lap + beta[0] * uu.dx + beta[1] * uu.dy;
                        double v = bDw * res;
                        val += (bwu + tau * v) * w * J;
                    }
                }

                if (val != 0 && ! is_boundary(a, x, y)) {
                    int i = linear_index(a, x, y) + 1;
                    int j = linear_index(b, x, y) + 1;
                    problem.add(i, j, val);
                }
            }
        }

        // 1's for Dirichlet BC
        for_boundary_dofs(x, y, [&](index_type dof) {
            int i = linear_index(dof, x, y) + 1;
            problem.add(i, i, 1);
        });
    }

    void prepare_matrices() {
        gram_matrix_1d(Mx, x.basis);
        gram_matrix_1d(My, y.basis);

        stiffness_matrix_1d(Kx, x.basis);
        stiffness_matrix_1d(Ky, y.basis);

        advection_matrix_1d(Ax, x.basis);
        advection_matrix_1d(Ay, y.basis);
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
        x.factorize_matrix();
        y.factorize_matrix();

        skew_bc(u, x, y);
        output.to_file(u, "out_0.data");
    }

    void step(int /*iter*/, double /*t*/) override {
        zero(rhs);
        compute_rhs();

        // stationary_bc(rhs, x, y);
        // skew_bc(rhs, x, y);
        zero_bc(rhs, x, y);

        mumps::problem problem(rhs.data(), rhs.size());
        assemble_problem(problem, steps.dt);
        solver.solve(problem);

        u = rhs;
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % save_every == 0) {
            std::cout << "Step " << (iter + 1) << " : " << errorL2() << " " << errorH1() << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data", u, x, y);
        std::cout << "{ 'L2': '" << errorL2() << "', 'H1': '" << errorH1() << "'}" << std::endl;
        print_solution("solution.data", u, x, y);
    }

    void compute_rhs() {
        zero(rhs);
        executor.for_each(elements(x, y), [&](index_type e) {
            auto U = vector_type{{ x.basis.dofs_per_element(), y.basis.dofs_per_element() }};

            double J = jacobian(e, x, y);
            for (auto q : quad_points(x, y)) {
                double w = weigth(q, x, y);
                // auto x = point(e, q);

                for (auto a : dofs_on_element(e, x, y)) {
                    auto aa = dof_global_to_local(e, a, x, y);
                    value_type v = eval_basis(e, q, a, x, y);

                    double F = 1;
                    double val = F * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e, x, y);
            });
        });
    }

    double errorL2() const {
        return Base::errorL2(u, x, y, exact(epsilon));
    }

    double errorH1() const {
        return Base::errorH1(u, x, y, exact(epsilon));
    }

};

}



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_SUPG_HPP */
