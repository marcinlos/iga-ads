#ifndef PROBLEMS_STOKES_STOKES_PROJECTION_HPP_
#define PROBLEMS_STOKES_STOKES_PROJECTION_HPP_

#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "problems/stokes/space_set.hpp"
#include "mumps.hpp"



namespace ads {

class stokes_projection : public simulation_2d {
private:
    using Base = simulation_2d;
    using value_pair = std::array<value_type, 2>;

    galois_executor executor{8};

    space_set trial, test;
    vector_type vx, vy, p;

    mumps::solver solver;
    output_manager<2> outputU1, outputU2, outputP;

public:
    stokes_projection(space_set trial_, space_set test_, const timesteps_config& steps)
    : Base{ test_.Px, test_.Py, steps }
    , trial{ std::move(trial_) }
    , test{ std::move(test_) }
    , vx{{ trial.U1x.dofs(), trial.U1y.dofs() }}
    , vy{{ trial.U2x.dofs(), trial.U2y.dofs() }}
    , p{{ trial.Px.dofs(), trial.Py.dofs() }}
    , outputU1{ trial.U1x.B, trial.U1y.B, 200 }
    , outputU2{ trial.U2x.B, trial.U2y.B, 200 }
    , outputP{ trial.Px.B, trial.Py.B, 200 }
    { }

    void before() override {
        trial.U1x.factorize_matrix();
        trial.U1y.factorize_matrix();
        trial.U2x.factorize_matrix();
        trial.U2y.factorize_matrix();
        trial.Px.factorize_matrix();
        trial.Py.factorize_matrix();

        auto project = [&](auto& rhs, auto& x, auto& y, auto fun) {
            vector_type buffer{{ x.dofs(), y.dofs() }};
            compute_projection(rhs, x.basis, y.basis, [&](double x, double y) { return fun({x, y}); });
            ads_solve(rhs, buffer, x.data(), y.data());
        };
        project(vx, trial.U1x, trial.U1y, [this](point_type x) { return exact_v(x, 0)[0].val; });
        project(vy, trial.U2x, trial.U2y, [this](point_type x) { return exact_v(x, 0)[1].val; });
        project(p, trial.Px, trial.Py, [this](point_type x) { return 0; /*exact_p(x, 0).val;*/ });

        save_to_file(0);
        // output_exact(0);

        // std::cout << "Initial state error:" << std::endl;
        // print_error(0);
    }

    value_type exact_p(point_type p, double t) const {
        auto x = p[0];
        auto et = std::exp(-t);
        return {et * (x * (1 - x) - 1./6), et * (1 - 2 * x), 0.0};
    }

    std::array<value_type, 2> exact_v(point_type p, double t) const {
        auto et = std::exp(-t);

        auto f = [](double x, double y) {
            return x*x * (1 - x) * (1 - x) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        auto dfx = [](double x, double y) {
            return (4 * x*x*x - 6 * x*x + 2 * x) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        auto dfy = [](double x, double y) {
            return x*x * (1 - x) * (1 - x) * (2 - 12 * y + 12 * y*y);
        };

        double x = p[0], y = p[1];
        value_type vx = {f(x, y), dfx(x, y), dfy(x, y)};
        value_type vy = {-f(y, x), -dfy(y, x), -dfx(y, x)};

        return { et * vx, et * vy };
    }

    value_type exact_div(point_type p, double t) const {
        auto et = std::exp(-t);
        auto v = exact_v(p, t);
        auto div = v[0].dx + v[1].dy;

        auto dfxy = [](double x, double y) {
            return (4 * x*x*x - 6 * x*x + 2 * x) * (2 - 12 * y + 12 * y*y);
        };

        auto dfxx = [](double x, double y) {
            return (12 * x*x - 12 * x + 2) * (2 * y - 6 * y*y + 4 * y*y*y);
        };

        double x = p[0], y = p[1];
        double dx = dfxx(x, y) - dfxy(y, x);
        double dy = dfxy(x, y) - dfxx(y, x);

        return { div, et * dx, et * dy };
    }

    void output_exact(double t) {
        auto p = [this,t](point_type x) { return exact_p(x, t).val; };
        auto vx = [this,t](point_type x) { return exact_v(x, t)[0].val; };
        auto vy = [this,t](point_type x) { return exact_v(x, t)[1].val; };

        auto project = [&](auto& x, auto& y, auto fun) {
            vector_type rhs{{ x.dofs(), y.dofs() }};
            vector_type buffer{{ x.dofs(), y.dofs() }};
            compute_projection(rhs, x.basis, y.basis, [&](double x, double y) { return fun({x, y}); });
            ads_solve(rhs, buffer, x.data(), y.data());
            return rhs;
        };

        outputP.to_file(project(trial.Px, trial.Py, p), "pressure_ref.data");
        outputU1.to_file(project(trial.U1x, trial.U1y, vx), "vx_ref.data");
        outputU2.to_file(project(trial.U2x, trial.U2y, vy), "vy_ref.data");
    }

    point_type forcing(point_type p, double t) const {
        double x = p[0], y = p[1];
        auto v = exact_v(p, t);
        auto et = std::exp(-t);

        auto fx =
            (12 - 24 * y) * x*x*x*x +
            (-24 + 48 * y) * x*x*x +
            (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
            (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
            1 - 4 * y + 12 * y*y - 8 * y*y*y;

        auto fy =
            (8 - 48 * y + 48 * y*y) * x*x*x +
            (-12 + 72 * y - 72 * y*y) * x*x +
            (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x -
            12 * y*y + 24 * y*y*y - 12 * y*y*y*y;

        return { et * fx - v[0].val, et * fy - v[1].val };
    }

    auto shifted(int n, int k, mumps::problem& problem) const {
        return [&problem,n,k](int i, int j, double val) {
            problem.add(n + i, k + j, val);
        };
    }


    bool is_pressure_fixed(index_type dof) const {
        return dof[0] == 0 && dof[1] == 0;
    }

    void assemble_matrix(mumps::problem& problem, double cx, double cy, bool bc,
                         const dimension& Ux, const dimension& Uy) const {
        for (auto i : dofs(Ux, Uy)) {
            for (auto j : overlapping_dofs(i, Ux, Uy)) {
                int ii = linear_index(i, Ux, Uy) + 1;
                int jj = linear_index(j, Ux, Uy) + 1;
                bool at_boundary = is_boundary(i, Ux, Uy) || is_boundary(j, Ux, Uy);

                if (!bc || !at_boundary) {
                    auto form = [cx,cy](auto u, auto v) { return u.val * v.val + cx * u.dx * v.dx + cy * u.dy * v.dy; };
                    auto product = integrate(i, j, Ux, Uy, Ux, Uy, form);
                    problem.add(ii, jj, product);
                }
            }
        }
        if (bc) {
            for_boundary_dofs(Ux, Uy, [&](index_type dof) {
                int i = linear_index(dof, Ux, Uy) + 1;
                problem.add(i, i, 1);
            });
        }
    }

    void assemble_matrix_velocity(mumps::problem& problem, double cx, double cy) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();

        auto D = DU1 + DU2;

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DU1, DU1, problem);

        auto trial_vx = shifted(D, D, problem);
        auto trial_vy = shifted(D + dU1, D + dU1, problem);

        // tx, vx -> (\/tx, \/vx)
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : overlapping_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                if (! is_boundary(i, test.U1x, test.U1y) && ! is_boundary(j, test.U1x, test.U1y)) {
                    double val = eval([this,cx,cy](auto tx, auto vx) {
                        return tx.val * vx.val + cx * tx.dx * vx.dx + cy * tx.dy * vx.dy;
                    });
                    test_vx(ii, jj, val);
                }
            }
        }

        // ty, vy -> (\/ty, \/vy)
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : overlapping_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                if (! is_boundary(i, test.U2x, test.U2y) && ! is_boundary(j, test.U2x, test.U2y)) {
                    double val = eval([this,cx,cy](auto ty, auto vy) {
                        return ty.val * vy.val + cx * ty.dx * vy.dx + cy * ty.dy * vy.dy;
                    });
                    test_vy(ii, jj, val);
                }
            }
        }

        // Strong BC
        for_boundary_dofs(test.U1x, test.U1y, [&](index_type dof) {
            int i = linear_index(dof, test.U1x, test.U1y) + 1;
            test_vx(i, i, 1);
        });
        for_boundary_dofs(test.U2x, test.U2y, [&](index_type dof) {
            int i = linear_index(dof, test.U2x, test.U2y) + 1;
            test_vy(i, i, 1);
        });

        // B, B^t
        auto put = [&](int i, int j, int si, int sj, double val, bool fixed_i, bool fixed_j) {
            int ii = i + si;
            int jj = j + sj;

            if (!fixed_i) {
                problem.add(ii, D + jj, val);
            }
            if (!fixed_i && !fixed_j) {
                problem.add(D + jj, ii, val);
            }
        };


        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! overlap(i, test.U1x, test.U1y, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form); };

                bool bd_i = is_boundary(i, test.U1x, test.U1y);
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y);

                double value = eval([cx,cy](auto u, auto v) { return u.val * v.val + cx * u.dx * v.dx + cy * u.dy * v.dy; });
                put(ii, jj, 0, 0, value, bd_i, bd_j);
            }
        }

        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! overlap(i, test.U2x, test.U2y, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form); };

                bool bd_i = is_boundary(i, test.U2x, test.U2y);
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y);

                double value = eval([cx,cy](auto u, auto v) { return u.val * v.val + cx * u.dx * v.dx + cy * u.dy * v.dy; });
                put(ii, jj, DU1, dU1, value, bd_i, bd_j);
            }
        }

        for_boundary_dofs(trial.U1x, trial.U1y, [&](index_type dof) {
            int i = linear_index(dof, trial.U1x, trial.U1y) + 1;
            trial_vx(i, i, 1);
        });
        for_boundary_dofs(trial.U2x, trial.U2y, [&](index_type dof) {
            int i = linear_index(dof, trial.U2x, trial.U2y) + 1;
            trial_vy(i, i, 1);
        });
    }

    void assemble_matrix_pressure(mumps::problem& problem, double cx, double cy) const {
        auto dP = trial.Px.dofs() * trial.Py.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();

        auto test_p = shifted(0, 0, problem);

        // Gram matrix
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : overlapping_dofs(i, test.Px, test.Py)) {
                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, test.Px, test.Py) + 1;

                if (! is_pressure_fixed(i) && ! is_pressure_fixed(j)) {
                    auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };
                    double val = eval([](auto w, auto p) { return w.val * p.val ; });
                    test_p(ii, jj, val);
                }
            }
        }

        // B, B^t
        auto put = [&](int i, int j, int si, int sj, double val) {
            int ii = i + si;
            int jj = j + sj;
            problem.add(ii, DP + jj, val);
            problem.add(DP + jj, ii, val);
        };

        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! overlap(i, test.Px, test.Py, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.Px, trial.Py, form); };

                double value = eval([cx,cy](auto u, auto v) { return u.val * v.val + cx * u.dx * v.dx + cy * u.dy * v.dy; });
                put(ii, jj, 0, 0, value);
            }
        }
    }


    void update_velocity_minev(double t) {
        using namespace std::placeholders;
        auto dt = steps.dt;
        auto f = [&](point_type x, double s) { return forcing(x, s); };
        auto F = [&](double s) { return std::bind(f, _1, s); };

        vector_type rhs_vx{{ trial.U1x.dofs(), trial.U1y.dofs() }};
        vector_type rhs_vy{{ trial.U2x.dofs(), trial.U2y.dofs() }};

        // Velocity - first equation
        compute_rhs(rhs_vx, rhs_vy, vx, vy, vx, vy, p, dt, F(t + dt/2), 0, 0, -dt, -dt, dt, dt);
        zero_bc(rhs_vx, trial.U1x, trial.U1y);
        zero_bc(rhs_vy, trial.U2x, trial.U2y);

        mumps::problem problem_vx1(rhs_vx.data(), rhs_vx.size());
        assemble_matrix(problem_vx1, 0, 0, true, trial.U1x, trial.U1y);
        solver.solve(problem_vx1);

        mumps::problem problem_vy1(rhs_vy.data(), rhs_vy.size());
        assemble_matrix(problem_vy1, 0, 0, true, trial.U2x, trial.U2y);
        solver.solve(problem_vy1);

        // Velocity - step 2
        vector_type rhs_vx2{{ trial.U1x.dofs(), trial.U1y.dofs() }};
        vector_type rhs_vy2{{ trial.U2x.dofs(), trial.U2y.dofs() }};
        compute_rhs(rhs_vx2, rhs_vy2, vx, vy, rhs_vx, rhs_vy, p, dt, F(t + dt/2), dt/2, 0, 0, 0, 0, 0);

        zero_bc(rhs_vx2, trial.U1x, trial.U1y);
        zero_bc(rhs_vy2, trial.U2x, trial.U2y);

        mumps::problem problem_vx2(rhs_vx2.data(), rhs_vx2.size());
        assemble_matrix(problem_vx2, dt/2, 0, true, trial.U1x, trial.U1y);
        solver.solve(problem_vx2);

        mumps::problem problem_vy2(rhs_vy2.data(), rhs_vy2.size());
        assemble_matrix(problem_vy2, dt/2, 0, true, trial.U2x, trial.U2y);
        solver.solve(problem_vy2);

        // Velocity - step 3
        vector_type rhs_vx3{{ trial.U1x.dofs(), trial.U1y.dofs() }};
        vector_type rhs_vy3{{ trial.U2x.dofs(), trial.U2y.dofs() }};
        compute_rhs(rhs_vx3, rhs_vy3, vx, vy, rhs_vx2, rhs_vy2, p, dt, F(t + dt/2), 0, dt/2, 0, 0, 0, 0);

        zero_bc(rhs_vx3, trial.U1x, trial.U1y);
        zero_bc(rhs_vy3, trial.U2x, trial.U2y);

        mumps::problem problem_vx3(rhs_vx3.data(), rhs_vx3.size());
        assemble_matrix(problem_vx3, 0, dt/2, true, trial.U1x, trial.U1y);
        solver.solve(problem_vx3);

        mumps::problem problem_vy3(rhs_vy3.data(), rhs_vy3.size());
        assemble_matrix(problem_vy3, 0, dt/2, true, trial.U2x, trial.U2y);
        solver.solve(problem_vy3);

        vx = rhs_vx3;
        vy = rhs_vy3;
    }

    void update_velocity_galerkin(double t) {
        using namespace std::placeholders;
        auto dt = steps.dt;
        auto f = [&](point_type x, double s) { return forcing(x, s); };
        auto F = [&](double s) { return std::bind(f, _1, s); };

        vector_type rhs_vx{{ trial.U1x.dofs(), trial.U1y.dofs() }};
        vector_type rhs_vy{{ trial.U2x.dofs(), trial.U2y.dofs() }};

        // Step 1
        compute_rhs(rhs_vx, rhs_vy, vx, vy, vx, vy, p, dt, F(t + dt/2), 0, 0, 0, - dt/2, dt/2, dt/2);
        zero_bc(rhs_vx, trial.U1x, trial.U1y);
        zero_bc(rhs_vy, trial.U2x, trial.U2y);

        mumps::problem problem_vx1(rhs_vx.data(), rhs_vx.size());
        assemble_matrix(problem_vx1, dt/2, 0, true, trial.U1x, trial.U1y);
        solver.solve(problem_vx1);

        mumps::problem problem_vy1(rhs_vy.data(), rhs_vy.size());
        assemble_matrix(problem_vy1, dt/2, 0, true, trial.U2x, trial.U2y);
        solver.solve(problem_vy1);

        // Step 2
        vector_type rhs_vx2{{ trial.U1x.dofs(), trial.U1y.dofs() }};
        vector_type rhs_vy2{{ trial.U2x.dofs(), trial.U2y.dofs() }};

        compute_rhs(rhs_vx2, rhs_vy2, vx, vy, rhs_vx, rhs_vy, p, dt, F(t + dt/2), 0, 0, -dt/2, 0, dt/2, dt/2);
        zero_bc(rhs_vx2, trial.U1x, trial.U1y);
        zero_bc(rhs_vy2, trial.U2x, trial.U2y);

        mumps::problem problem_vx2(rhs_vx2.data(), rhs_vx2.size());
        assemble_matrix(problem_vx2, 0, dt/2, true, trial.U1x, trial.U1y);
        solver.solve(problem_vx2);

        mumps::problem problem_vy2(rhs_vy2.data(), rhs_vy2.size());
        assemble_matrix(problem_vy2, 0, dt/2, true, trial.U2x, trial.U2y);
        solver.solve(problem_vy2);

        vx = rhs_vx2;
        vy = rhs_vy2;
    }

    void update_velocity_igrm(double t) {
        using namespace std::placeholders;
        auto dt = steps.dt;
        auto f = [&](point_type x, double s) { return forcing(x, s); };
        auto F = [&](double s) { return std::bind(f, _1, s); };

        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dim_trial = dU1 + dU2;

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto dim_test = DU1 + DU2;

        // Step 1
        std::vector<double> rhs(dim_test + dim_trial);
        vector_view rhs_vx1{rhs.data(),         {test.U1x.dofs(), test.U1y.dofs()}};
        vector_view rhs_vy1{rhs.data() + DU1,   {test.U2x.dofs(), test.U2y.dofs()}};
        vector_view vx1{rhs.data() + dim_test,  {trial.U1x.dofs(), trial.U1y.dofs()}};
        vector_view vy1{vx1.data() + dU1,       {trial.U2x.dofs(), trial.U2y.dofs()}};

        compute_rhs(rhs_vx1, rhs_vy1, vx, vy, vx, vy, p, dt, F(t + dt/2), 0, 0, 0, - dt/2, dt/2, dt/2);

        mumps::problem problem_vx1(rhs.data(), rhs.size());
        assemble_matrix_velocity(problem_vx1, dt/2, 0);
        solver.solve(problem_vx1);

        // Step 2
        std::vector<double> rhs2(dim_test + dim_trial);
        vector_view rhs_vx2{rhs2.data(),        {test.U1x.dofs(), test.U1y.dofs()}};
        vector_view rhs_vy2{rhs2.data() + DU1,  {test.U2x.dofs(), test.U2y.dofs()}};
        vector_view vx2{rhs2.data() + dim_test, {trial.U1x.dofs(), trial.U1y.dofs()}};
        vector_view vy2{vx2.data() + dU1,       {trial.U2x.dofs(), trial.U2y.dofs()}};

        compute_rhs(rhs_vx2, rhs_vy2, vx, vy, vx1, vy1, p, dt, F(t + dt/2), 0, 0, -dt/2, 0, dt/2, dt/2);

        mumps::problem problem_vx2(rhs2.data(), rhs2.size());
        assemble_matrix_velocity(problem_vx2, 0, dt/2);
        solver.solve(problem_vx2);

        // Update
        for (auto i : dofs(trial.U1x, trial.U1y)) { vx(i[0], i[1]) = vx2(i[0], i[1]); } // vx = vx2;
        for (auto i : dofs(trial.U2x, trial.U2y)) { vy(i[0], i[1]) = vy2(i[0], i[1]); } // vy = vy2;
    }

    void update_pressure() {
        vector_type rhs_p{{ trial.Px.dofs(), trial.Py.dofs() }};

        // Pressure
        compute_rhs_pressure_1(rhs_p, vx, vy, trial.Px, trial.Py, steps.dt);
        mumps::problem problem_px(rhs_p.data(), rhs_p.size());
        assemble_matrix(problem_px, 1, 0, false, trial.Px, trial.Py);
        solver.solve(problem_px);

        zero(p);
        compute_rhs_pressure_2(p, rhs_p, trial.Px, trial.Py, steps.dt);
        mumps::problem problem_py(p.data(), rhs_p.size());
        assemble_matrix(problem_py, 0, 1, false, trial.Px, trial.Py);
        solver.solve(problem_py);
    }

    void update_pressure_igrm() {
        auto dim_trial = trial.Px.dofs() * trial.Py.dofs();
        auto dim_test = test.Px.dofs() * test.Py.dofs();

        // Step 1
        std::vector<double> rhs(dim_test + dim_trial);
        vector_view rhs_p1{rhs.data(),         {test.Px.dofs(), test.Py.dofs()}};
        vector_view p1{rhs.data() + dim_test, {trial.Px.dofs(), trial.Py.dofs()}};

        compute_rhs_pressure_1(rhs_p1, vx, vy, test.Px, test.Py, steps.dt);
        mumps::problem problem_px(rhs.data(), rhs.size());
        assemble_matrix_pressure(problem_px, 1, 0);
        solver.solve(problem_px);

        // Step 2
        std::vector<double> rhs2(dim_test + dim_trial);
        vector_view rhs_p2{rhs2.data(),         {test.Px.dofs(), test.Py.dofs()}};
        vector_view p2{rhs2.data() + dim_test, {trial.Px.dofs(), trial.Py.dofs()}};

        compute_rhs_pressure_2(rhs_p2, p1, test.Px, test.Py, steps.dt);
        mumps::problem problem_py(rhs2.data(), rhs2.size());
        assemble_matrix_pressure(problem_py, 0, 1);
        solver.solve(problem_py);

        for (auto i : dofs(trial.Px, trial.Py)) { p(i[0], i[1]) = p2(i[0], i[1]); } // p = p2;
    }

    void step(int iter, double t) override {
        // update_velocity_galerkin(t);
        // update_velocity_minev(t);
        update_velocity_igrm(t);

        // update_pressure();
        update_pressure_igrm();
    }


    // Compute RHS as
    //
    // (v, w) + ax (dv0/dx, dw/dx) + ay (dv0/dy, dw/dy) +
    //          bx (dv/dx, dw/dx) + by (dv/dy, dw/dy) +
    //          c (\/p, w) + d  (f, w)
    template <typename RHS, typename S1, typename S2, typename S3, typename Fun>
    void compute_rhs(RHS& rhsx, RHS& rhsy,
                     const S1& vx0, const S1& vy0,
                     const S2& vx, const S2& vy, const S3& p,
                     double dt, Fun&& forcing,
                     double ax, double ay, double bx, double by, double c, double d) const {
        using shape = std::array<std::size_t, 2>;
        auto u1_shape = shape{ test.U1x.basis.dofs_per_element(), test.U1y.basis.dofs_per_element() };
        auto u2_shape = shape{ test.U2x.basis.dofs_per_element(), test.U2y.basis.dofs_per_element() };

        executor.for_each(elements(trial.Px, trial.Py), [&](index_type e) {
            auto vx_loc = vector_type{ u1_shape };
            auto vy_loc = vector_type{ u2_shape };

            double J = jacobian(e);
            for (auto q : quad_points(trial.Px, trial.Py)) {
                double W = weigth(q);
                auto x = point(e, q);
                auto F = forcing(x);
                value_type vvx0 = eval(vx0, e, q, trial.U1x, trial.U1y);
                value_type vvy0 = eval(vy0, e, q, trial.U2x, trial.U2y);
                value_type vvx = eval(vx, e, q, trial.U1x, trial.U1y);
                value_type vvy = eval(vy, e, q, trial.U2x, trial.U2y);
                value_type pp = eval(p, e, q, trial.Px, trial.Py);

                for (auto a : dofs_on_element(e, test.U1x, test.U1y)) {
                    auto aa = dof_global_to_local(e, a, test.U1x, test.U1y);
                    value_type v = eval_basis(e, q, a, test.U1x, test.U1y);

                    double val = vvx.val * v.val
                        + ax * vvx0.dx * v.dx
                        + ay * vvx0.dy * v.dy
                        + bx * vvx.dx * v.dx
                        + by * vvx.dy * v.dy
                        + c * pp.val * v.dx
                        + d * F[0] * v.val;
                    vx_loc(aa[0], aa[1]) += val * W * J;
                }
                for (auto a : dofs_on_element(e, test.U2x, test.U2y)) {
                    auto aa = dof_global_to_local(e, a, test.U2x, test.U2y);
                    value_type v = eval_basis(e, q, a, test.U2x, test.U2y);

                    double val = vvy.val * v.val
                        + ax * vvy0.dx * v.dx
                        + ay * vvy0.dy * v.dy
                        + bx * vvy.dx * v.dx
                        + by * vvy.dy * v.dy
                        + c * pp.val * v.dy
                        + d * F[1] * v.val;
                    vy_loc(aa[0], aa[1]) += val * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhsx, vx_loc, e, test.U1x, test.U1y);
                update_global_rhs(rhsy, vy_loc, e, test.U2x, test.U2y);
            });
        });
    }

    template <typename RHS, typename Sol>
    void compute_rhs_pressure_1(RHS& rhs, const Sol& vx, const Sol& vy,
                                const dimension& Vx, const dimension& Vy, double dt) const {
        using shape = std::array<std::size_t, 2>;
        auto p_shape = shape{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() };

        executor.for_each(elements(trial.Px, trial.Py), [&](index_type e) {
            auto loc = vector_type{ p_shape };

            double J = jacobian(e);
            for (auto q : quad_points(trial.Px, trial.Py)) {
                double W = weigth(q);
                value_type vvx = eval(vx, e, q, trial.U1x, trial.U1y);
                value_type vvy = eval(vy, e, q, trial.U2x, trial.U2y);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    double val = - 1/dt * (vvx.dx + vvy.dy) * v.val;
                    loc(aa[0], aa[1]) += val * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs, loc, e, Vx, Vy);
            });
        });
    }

    template <typename RHS, typename Sol>
    void compute_rhs_pressure_2(RHS& rhs, const Sol& p,
                                const dimension& Vx, const dimension& Vy, double dt) const {
        using shape = std::array<std::size_t, 2>;
        auto p_shape = shape{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() };

        executor.for_each(elements(trial.Px, trial.Py), [&](index_type e) {
            auto loc = vector_type{ p_shape };

            double J = jacobian(e);
            for (auto q : quad_points(trial.Px, trial.Py)) {
                double W = weigth(q);
                value_type pp = eval(p, e, q, trial.Px, trial.Py);

                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    double val = pp.val * v.val;
                    loc(aa[0], aa[1]) += val * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(rhs, loc, e, Vx, Vy);
            });
        });
    }

    void after_step(int iter, double t) override {
        int i = iter + 1;
        double tt = t + steps.dt;

        if (i % 1 == 0) {
            // std::cout << "Outputting" << std::endl;
            save_to_file(i);
        }
        // auto p_avg = correct_pressure(p);

        auto e_vx = [this,tt](point_type x) { return exact_v(x, tt)[0]; };
        auto e_vy = [this,tt](point_type x) { return exact_v(x, tt)[1]; };
        auto e_p = [this,tt](point_type x) { return exact_p(x, tt); };

        double vxL2 = errorL2(vx, trial.U1x, trial.U1y, e_vx) / normL2(trial.U1x, trial.U1y, e_vx) * 100;
        double vyL2 = errorL2(vy, trial.U2x, trial.U2y, e_vy) / normL2(trial.U2x, trial.U2y, e_vy) * 100;
        double vL2 = std::hypot(vxL2, vyL2) / std::sqrt(2);

        double vxH1 = errorH1(vx, trial.U1x, trial.U1y, e_vx) / normH1(trial.U1x, trial.U1y, e_vx) * 100;
        double vyH1 = errorH1(vy, trial.U2x, trial.U2y, e_vy) / normH1(trial.U2x, trial.U2y, e_vy) * 100;
        double vH1 = std::hypot(vxH1, vyH1) / std::sqrt(2);

        double pL2 = errorL2(p, trial.Px, trial.Py, e_p) / normL2(trial.Px, trial.Py, e_p) * 100;
        double pH1 = errorH1(p, trial.Px, trial.Py, e_p) / normH1(trial.Px, trial.Py, e_p) * 100;

        std::cout << i << " " << tt
                  << "  vx: " << vxL2 << " " << vxH1
                  << "  vy: " << vyL2 << " " << vyH1
                  // << "  p avg: " << p_avg
                  << "  p:  " << pL2 << " " << pH1 << std::endl;
    }

    template <typename Sol>
    double correct_pressure(Sol& pressure) const {
        auto p_avg = average_value(pressure, trial.Px, trial.Py);
        for (auto i : dofs(trial.Px, trial.Py)) {
            pressure(i[0], i[1]) -= p_avg;
        }
        return p_avg;
    }

    template <typename Sol>
    double average_value(const Sol& u, const dimension& Ux, const dimension& Uy) const {
        double val = 0;
        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q, Ux, Uy);
                value_type uu = eval(u, e, q, Ux, Uy);
                val += uu.val * w * J;
            }
        }
        return val;
    }


    void save_to_file(int i) {
        outputP.to_file(p, "pressure_%d.data", i);
        outputU1.to_file(vx, "vx_%d.data", i);
        outputU2.to_file(vy, "vy_%d.data", i);
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) const {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

    template <typename Sol, typename Fun>
    double div_errorL2(const Sol& u, const Sol& v, const space_set& space, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return div_error(u, v, space, L2, fun);
    }

    template <typename Sol, typename Fun>
    double div_errorH1(const Sol& u, const Sol& v, const space_set& space, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy; };
        return div_error(u, v, space, H1, fun);
    }

    template <typename Sol, typename Fun, typename Norm>
    double div_error(const Sol& u, const Sol& v, const space_set& space, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(space.Px, space.Px)) {
            double J = jacobian(e, space.Px, space.Py);
            for (auto q : quad_points(space.Px, space.Py)) {
                double w = weigth(q, space.Px, space.Py);
                auto x = point(e, q, space.Px, space.Py);
                value_type div = divergence(u, v, e, q, space);

                auto d = div - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol>
    value_type divergence(const Sol& u, const Sol& v, index_type e, index_type q, const space_set& space) const {
        value_type div{};
        for (auto b : dofs_on_element(e, space.U1x, space.U1y)) {
            double c = u(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, space.U1x, space.U1y);

            const auto& bx = space.U1x.basis;
            const auto& by = space.U1y.basis;

            double B2  = by.b[e[1]][q[1]][0][loc[1]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];

            double dx = dB1 *  B2;
            double dxx = ddB1 * B2;
            double dxy = dB1 * dB2;

            div.val += c * dx;
            div.dx += c * dxx;
            div.dy += c * dxy;
        }
        for (auto b : dofs_on_element(e, space.U2x, space.U2y)) {
            double d = v(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, space.U2x, space.U2y);

            const auto& bx = space.U2x.basis;
            const auto& by = space.U2y.basis;

            double B1  = bx.b[e[0]][q[0]][0][loc[0]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB2 = by.b[e[1]][q[1]][2][loc[1]];

            double dy =  B1 * dB2;
            double dyy = B1 * ddB2;
            double dxy = dB1 * dB2;

            div.val += d * dy;
            div.dx += d * dxy;
            div.dy += d * dyy;
        }
        return div;
    }

    bool overlap(int a, const dimension& U, int b, const dimension& V) const {
        auto ar = U.basis.element_ranges[a];
        auto br = V.basis.element_ranges[b];
        return (ar.first >= br.first && ar.first <= br.second) || (br.first >= ar.first && br.first <= ar.second);
    }

    bool overlap(index_type a, const dimension& Ux, const dimension& Uy,
                 index_type b, const dimension& Vx, const dimension& Vy) const {
        return overlap(a[0], Ux, b[0], Vx) && overlap(a[1], Uy, b[1], Vy);
    }

    double dot(value_type a, point_type b) const {
        return a.dx * b[0] + a.dy * b[1];
    }

    bool supported_in_1d(int dof, int e, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return e >= xrange.first && e <= xrange.second;
    }

    bool touches_point(int dof, int sx, const dimension& x) const {
        auto xrange = x.basis.element_ranges[dof];
        return sx >= xrange.first && sx <= xrange.second + 1;
    }

    template <typename Form>
    double integrate(index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                     const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;

        for (auto e : elements_supporting_dof(i, Ux, Uy)) {
            if (! supported_in(j, e, Vx, Vy)) continue;

            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q);
                value_type ww = eval_basis(e, q, i, Ux, Uy);
                value_type uu = eval_basis(e, q, j, Vx, Vy);

                double fuw = form(ww, uu);
                val += fuw * w * J;
            }
        }
        return val;
    }

};

}

#endif // PROBLEMS_STOKES_STOKES_PROJECTION_HPP_
