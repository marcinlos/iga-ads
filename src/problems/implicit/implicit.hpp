#ifndef PROBLEMS_IMPLICIT_HPP_
#define PROBLEMS_IMPLICIT_HPP_

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"


namespace ads {

class implicit_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    output_manager<2> output;
    galois_executor executor{4};

    lin::band_matrix Kx, Ky;

public:
    implicit_2d(const config_2d& config)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , output{ x.B, y.B, 200 }
    , Kx{x.p, x.p, x.B.dofs()}
    , Ky{y.p, y.p, y.B.dofs()}
    {
        matrix(Kx, x.basis, steps.dt);
        matrix(Ky, y.basis, steps.dt);
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(12 * (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    };

private:
    void matrix(lin::band_matrix& K, const basis_data& d, double h) {
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
                        auto da = d.b[e][q][1][a];
                        auto db = d.b[e][q][1][b];
                        K(ia, ib) += (va * vb + 0.5 * h * da * db) * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void solve(vector_type& v) {
        // lin::vector buf{{ y.dofs() }};
        // compute_projection(buf, y.basis, [](double y) {
        //     return std::sin(y * M_PI);
        // });
        // for (int i = 0; i < y.dofs(); ++ i) {
        //     v(0, i) = buf(i);
        // }
        Base::solve(v);
    }

    void prepare_matrices() {
        // x.fix_left();
        Base::prepare_matrices();
        lin::factorize(Kx, x.ctx);
        lin::factorize(Ky, y.ctx);
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);

        output.to_file(u, "out_0.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs_1();
        ads_solve(u, buffer, dim_data{Kx, x.ctx}, y.data());

        using std::swap;
        swap(u, u_prev);

        compute_rhs_2();
        ads_solve(u, buffer, x.data(), dim_data{Ky, y.ctx});
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % 1 == 0) {
            output.to_file(u, "out_%d.data", (iter + 1) / 1);
        }
    }

    void compute_rhs_1() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);

                    double gradient_prod = u.dy * v.dy;
                    double val = u.val * v.val - 0.5 * steps.dt * gradient_prod;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    void compute_rhs_2() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type u = eval_fun(u_prev, e, q);

                    double gradient_prod = u.dx * v.dx;
                    double val = u.val * v.val - 0.5 * steps.dt * gradient_prod;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }
};

}

#endif /* PROBLEMS_IMPLICIT_HPP_
 */
