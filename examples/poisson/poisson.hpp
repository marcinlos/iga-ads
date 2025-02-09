// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef POISSON_POISSON_HPP
#define POISSON_POISSON_HPP

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"

namespace ads::problems {

class poisson_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u;

    output_manager<2> output;
    galois_executor executor{4};
    ads::mumps::solver solver;

public:
    explicit poisson_2d(const config_2d& config)
    : Base{config}
    , u{shape()}
    , output{x.B, y.B, 200} { }

private:
    void solve(vector_type& v) {
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();
        Base::prepare_matrices();
    }

    void before() override { prepare_matrices(); }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs();
        dirichlet_bc(u, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(u, boundary::right, x, y, [](double t) { return 0; });
        dirichlet_bc(u, boundary::top, x, y, [](double t) { return 0; });
        dirichlet_bc(u, boundary::bottom, x, y, [](double t) { return 0; });

        ads::mumps::problem problem_u(u.data(), u.size());
        assemble_problem(problem_u);
        solver.solve(problem_u);
    }

    void compute_rhs() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto const [x, y] = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    // auto fval = 16 * 2 * (x * (1 - x) + y * (1 - y));

                    int const n = 8;
                    int const m = 11;
                    auto const pi = M_PI;

                    auto const lambda = pi * pi * (n * n + m * m);
                    auto fval = lambda * std::sin(n * pi * x) * std::sin(m * pi * y);

                    double val = fval * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    void after() override { output.to_file(u, "u.data"); }

    void assemble_problem(ads::mumps::problem& problem) {
        for (auto a : dofs(x, y)) {
            for (auto b : overlapping_dofs(a, x, y)) {
                if (is_fixed(a, x, y))
                    continue;

                double val = 0;
                for (auto e : elements_supporting_dof(a, x, y)) {
                    if (!supported_in(b, e, x, y))
                        continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(x, y)) {
                        double w = weight(q, x, y);
                        value_type ww = eval_basis(e, q, a, x, y);
                        value_type uu = eval_basis(e, q, b, x, y);

                        double bwu = grad_dot(uu, ww);
                        val += bwu * w * J;
                    }
                }

                if (val != 0) {
                    int i = linear_index(a, x, y) + 1;
                    int j = linear_index(b, x, y) + 1;
                    problem.add(i, j, val);
                }
            }
        }

        // 1's for Dirichlet BC
        for_boundary_dofs(x, y, [&](index_type dof) {
            if (is_fixed(dof, x, y)) {
                int i = linear_index(dof, x, y) + 1;
                problem.add(i, i, 1);
            }
        });
    }

    bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
        return dof[0] == 0 || dof[0] == x.dofs() - 1 || dof[1] == 0 || dof[1] == y.dofs() - 1;
    }
};

}  // namespace ads::problems

#endif  // POISSON_POISSON_HPP
