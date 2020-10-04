#ifndef PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP_

#include "problems/erikkson/erikkson_base.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"


#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/dense_solve.hpp"

#include <galois/Timer.h>


namespace ads {

class erikkson_mumps : public erikkson_base {
private:
    using Base = erikkson_base;

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

    double h;

    double peclet = 1e2;
    double epsilon = 1 / peclet;

    point_type c_diff{{ epsilon, epsilon }};

    double angle = 0;
    double len = 1;

    // point_type beta{{ len * cos(angle), len * sin(angle) }};
    point_type beta{{ 1, 1 }};

    galois::StatTimer solver_timer{"solver"};
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
    , h{ element_diam(Ux, Uy) }
    , output{ Ux.B, Uy.B, 500 }
    { }


private:

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    void assemble_problem(mumps::problem& problem, double dt) {
        auto N = Vx.dofs() * Vy.dofs();
        auto hh = h * h;

        // Gram matrix - upper left
        for (auto i : internal_dofs(Vx, Vy)) {
            for (auto j : overlapping_internal_dofs(i, Vx, Vy)) {
                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double val = kron(MVx, MVy, i, j) + hh * (kron(KVx, MVy, i, j) + kron(MVx, KVy, i, j));
                problem.add(ii, jj, val);
            }
        }

        // B, B^T
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : dofs(Ux, Uy)) {
                double val = 0;
                // val += kron(MUVx, MUVy, i, j);
                // val += steps.dt * (c_diff[0] * kron(KUVx, MUVy, i, j) + beta[0] * kron(AUVx, MUVy, i, j));
                // val += steps.dt * (c_diff[1] * kron(MUVx, KUVy, i, j) + beta[1] * kron(MUVx, AUVy, i, j));
                val += c_diff[0] * kron(KUVx, MUVy, i, j) + beta[0] * kron(AUVx, MUVy, i, j);
                val += c_diff[1] * kron(MUVx, KUVy, i, j) + beta[1] * kron(MUVx, AUVy, i, j);

                if (val != 0) {
                    int ii = linear_index(i, Vx, Vy) + 1;
                    int jj = linear_index(j, Ux, Uy) + 1;

                    if (! is_boundary(i, Vx, Vy)) {
                        problem.add(ii, N + jj, -val);
                    }
                    if (! is_boundary(i, Vx, Vy) && ! is_boundary(j, Ux, Uy)) {
                        problem.add(N + jj, ii, val);
                    }
                }
            }
        }

        // Dirichlet BC - upper left
        for_boundary_dofs(Vx, Vy, [&](index_type dof) {
            int i = linear_index(dof, Vx, Vy) + 1;
            problem.add(i, i, 1);
        });

        // Dirichlet BC - lower right
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
    }

    void before() override {
        prepare_matrices();
        Ux.factorize_matrix();
        Uy.factorize_matrix();

        // auto init = [this](double x, double y) { return init_state(x, y); };
        // compute_projection(u, Ux.basis, Uy.basis, init);
        // ads_solve(u, u_buffer, Ux.data(), Uy.data());

        // stationary_bc(u, Ux, Uy);
        zero_bc(u, Ux, Uy);

        output.to_file(u, "out_0.data");
    }

    void step(int /*iter*/, double t) override {
        mumps::problem problem(full_rhs.data(), full_rhs.size());

        std::cout << "Assembling matrix" << std::endl;
        assemble_problem(problem, steps.dt);

        std::cout << "Computing RHS" << std::endl;
        compute_rhs(t);
        // zero(rhs);

        std::fill(begin(full_rhs), end(full_rhs), 0);
        vector_view view_in{ full_rhs.data(), {Vx.dofs(), Vy.dofs()}};
        vector_view view_out{ full_rhs.data() + view_in.size(), {Ux.dofs(), Uy.dofs()}};

        for (auto i : dofs(Vx, Vy)) {
            view_in(i[0], i[1]) = -rhs(i[0], i[1]);
        }
        // stationary_bc(view_out, Ux, Uy);
        zero_bc(view_out, Ux, Uy);

        std::cout << "Solving" << std::endl;
        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        for (auto i : dofs(Ux, Uy)) {
            u(i[0], i[1]) = view_out(i[0], i[1]);
        }

        std::cout << "  solver time:       " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
        std::cout << "  assembly    FLOPS: " << solver.flops_assembly() << std::endl;
        std::cout << "  elimination FLOPS: " << solver.flops_elimination() << std::endl;

        std::cout << "Error: L2 = " << errorL2(t) << "%, H1 =  " << errorH1(t) << "%" << std::endl;
    }

    void after_step(int iter, double t) override {
        if ((iter + 1) % save_every == 0) {
            // std::cout << "Step " << (iter + 1) << " : " << errorL2(t) << "% " << errorH1(t) << "%" << std::endl;
            // std::cout << iter << " " << t << " " << errorL2(t) << " " << errorH1(t) << std::endl;
            output.to_file(u, "out_%d.data", (iter + 1) / save_every);
        }
    }

    void after() override {
        plot_middle("final.data", u, Ux, Uy);
        // std::cout << "{ 'L2': '" << errorL2() << "', 'H1': '" << errorH1() << "'}" << std::endl;
        print_solution("solution.data", u, Ux, Uy);
    }

    void compute_rhs(double t) {
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

                    // double F = 0;
                    // double F = erikkson_forcing(x[0], x[1], epsilon, t + steps.dt);
                    double F = erikkson2_forcing(x[0], x[1], epsilon);
                    // double val = (uu.val + steps.dt * F) * v.val;
                    double val = F * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e, Vx, Vy);
            });
        });
    }

    // double errorL2() const {
    //     return Base::errorL2(u, Ux, Uy, exact(epsilon));
    // }

    // double errorH1() const {
    //     return Base::errorH1(u, Ux, Uy, exact(epsilon));
    // }

    double errorL2(double t) const {
        // auto sol = exact(epsilon);
        auto sol = [&](point_type x) { return erikkson2_exact(x[0], x[1], epsilon); };
        // auto sol = [&](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };

        return Base::errorL2(u, Ux, Uy, sol) / normL2(Ux, Uy, sol) * 100;
    }

    double errorH1(double t) const {
        // auto sol = exact(epsilon);
        auto sol = [&](point_type x) { return erikkson2_exact(x[0], x[1], epsilon); };
        // auto sol = [&](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };

        return Base::errorH1(u, Ux, Uy, sol) / normH1(Ux, Uy, sol) * 100;
    }

};

}



#endif /* ADS_PROBLEMS_ERIKKSON_ERIKKSON_MUMPS_HPP */
