#ifndef ADS_PROBLEMS_STOKES_STOKES_HPP
#define ADS_PROBLEMS_STOKES_STOKES_HPP

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"
#include "ads/lin/tensor/view.hpp"
#include "mumps.hpp"


namespace ads {
namespace problems {

class stokes_split : public simulation_2d {
private:
    using Base = simulation_2d;
    using vector_view = lin::tensor_view<double, 2>;


    struct state {
        vector_type vx, vy;
        vector_type p;

        state(std::array<std::size_t, 2> shape)
        : vx{{ shape[0] + 1, shape[1] }}
        , vy{{ shape[0], shape[1] + 1 }}
        , p{ shape }
        { }
    };

    state now;

    vector_type u, u_prev;
    dimension x_1, y_1; // p + 1

    output_manager<2> output_p, output_vx, output_vy;
    galois_executor executor{4};

    lin::band_matrix Kx, Ky;
    lin::band_matrix Kx_1, Ky_1;

    // Aij - i-th factor has higher degree, j-th factor under derivative
    lin::band_matrix A11x, A12x, A21x, A22x;
    lin::band_matrix A11y, A12y, A21y, A22y;

    mumps::solver solver;

    static constexpr std::size_t RES = 400;

public:
    stokes_split(const config_2d& config)
    : Base{ config }
    , now{ shape() }
    , u{ shape() }
    , u_prev{ shape() }
    , x_1{ higher(config.x), config.derivatives }
    , y_1{ higher(config.y), config.derivatives }
    , output_p{ x.B, y.B, RES }
    , output_vx{ x_1.B, y.B, RES }
    , output_vy{ x.B, y_1.B, RES }
    , Kx{x.p, x.p, x.B.dofs()}
    , Ky{y.p, y.p, y.B.dofs()}
    , Kx_1{x_1.p, x_1.p, x_1.B.dofs()}
    , Ky_1{y_1.p, y_1.p, y_1.B.dofs()}
    , A11x{x_1.p, x.p, x_1.dofs(), x.dofs()}
    , A12x{x_1.p, x.p, x_1.dofs(), x.dofs()}
    , A21x{x.p, x_1.p, x.dofs(), x_1.dofs()}
    , A22x{x.p, x_1.p, x.dofs(), x_1.dofs()}
    , A11y{y_1.p, y.p, y_1.dofs(), y.dofs()}
    , A12y{y_1.p, y.p, y_1.dofs(), y.dofs()}
    , A21y{y.p, y_1.p, y.dofs(), y_1.dofs()}
    , A22y{y.p, y_1.p, y.dofs(), y_1.dofs()}
    {
        stiffness(Kx, x.basis);
        stiffness(Ky, y.basis);
        stiffness(Kx_1, x_1.basis);
        stiffness(Ky_1, y_1.basis);

        mixed1(A11x, x_1.basis, x.basis);
        mixed2(A12x, x_1.basis, x.basis);
        mixed1(A21x, x.basis, x_1.basis);
        mixed2(A22x, x.basis, x_1.basis);

        mixed1(A11y, y_1.basis, y.basis);
        mixed2(A12y, y_1.basis, y.basis);
        mixed1(A21y, y.basis, y_1.basis);
        mixed2(A22y, y.basis, y_1.basis);
    }

    dim_config higher(dim_config cfg) const {
        ++ cfg.p;
        // ++ cfg.quad_order; // TODO: handle this correctly
        return cfg;
    }

    double init_state(double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(12 * (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    };

private:
    void stiffness(lin::band_matrix& K, const basis_data& d) {
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
                        K(ia, ib) += da * db * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void mixed2(lin::band_matrix& K, const basis_data& bV, const basis_data& bU) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto va = bV.b[e][q][0][a];
                        auto db = bU.b[e][q][1][b];
                        K(ia, ib) += va * db * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void mixed1(lin::band_matrix& K, const basis_data& bV, const basis_data& bU) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto da = bV.b[e][q][1][a];
                        auto vb = bU.b[e][q][0][b];
                        K(ia, ib) += da * vb * bV.w[q] * bV.J[e];
                    }
                }
            }
        }
    }

    void advection_matrix(lin::dense_matrix& M, const basis_data& bU, const basis_data& bV,
                          double h, double advection) {
        for (element_id e = 0; e < bV.elements; ++ e) {
            for (int q = 0; q < bV.quad_order; ++ q) {
                for (int a = 0; a + bV.first_dof(e) <= bV.last_dof(e); ++ a) {
                    for (int b = 0; b + bU.first_dof(e) <= bU.last_dof(e); ++ b) {
                        int ia = a + bV.first_dof(e);
                        int ib = b + bU.first_dof(e);
                        auto va = bV.b[e][q][0][a];
                        auto db = bU.b[e][q][1][b];
                        auto diff = advection * h * va * db;
                        M(ia, ib) += diff * bV.w[q] * bV.J[e];
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
    }

    template <typename Rhs>
    void tensor_product(const Rhs& rhs1, const Rhs& rhs2, mumps::problem& p,
                        const lin::band_matrix& A, const lin::band_matrix& B,
                        double a = 1.0, int bx = 0, int by = 0) const {
        using std::min;
        using std::max;

        for (int ix = 0; ix < A.rows; ++ ix) {
            for (int jx = max(0, ix - A.kl); jx < min(A.rows, ix + A.ku + 1); ++ jx) {
                for (int iy = 0; iy < B.rows; ++ iy) {
                    for (int jy = max(0, iy - B.kl); jy < min(B.rows, iy + B.ku + 1); ++ jy) {
                        int i = &rhs1(ix, iy) - &rhs1(0, 0) + 1;
                        int j = &rhs2(jx, jy) - &rhs2(0, 0) + 1;
                        double val = a * A(ix, jx) * B(iy, jy);
                        p.add(bx + i, by + j, val);
                    }
                }
            }
        }
    }

    template <typename Rhs>
    void solve(Rhs& rhs, const lin::band_matrix& A, const lin::band_matrix& B) {
        mumps::problem problem(rhs.data(), rhs.size());
        tensor_product(rhs, rhs, problem, A, B);
        solver.solve(problem);
    }

    void before() override {
        prepare_matrices();

        compute_exact();
        compute();

        double E = energy();
        double L2 = L2norm();
        std::cout << 0 << " " << E << " " << L2 << std::endl;
    }

    void compute_exact() {
        auto e_pressure = [this](double x, double /*y*/) {
            return x * (1 - x);
        };

        auto fun = [](double x, double y) {
            return x * x * (1 - x) * (1 - x) * (2 * y - 6 * y * y + 4 * y * y * y);
        };
        auto e_vx = [&](double x, double y) { return fun(x, y); };
        auto e_vy = [&](double x, double y) { return -fun(y, x); };

        std::vector<double> rhs(now.vx.size() + now.vy.size() + now.p.size());

        vector_view vx{ rhs.data(), now.vx.sizes() };
        vector_view vy{ vx.data() + vx.size(), now.vy.sizes() };
        vector_view p{ vy.data() + vy.size(), now.p.sizes() };

        compute_projection(vx, x_1.basis, y.basis, e_vx);
        compute_projection(vy, x.basis, y_1.basis, e_vy);
        compute_projection(p, x.basis, y.basis, e_pressure);

        mumps::problem problem(rhs.data(), rhs.size());

        auto nvx = vx.data() - rhs.data();
        auto nvy = vy.data() - rhs.data();
        auto np = p.data() - rhs.data();

        tensor_product(vx, vx,problem, x_1.M, y.M);
        tensor_product(vy, vy, problem, x.M, y_1.M, 1, nvy, nvy);
        tensor_product(p, p, problem, x.M, y.M, 1, np, np);

        solver.solve(problem);

        output_p.to_file(p, "pressure_ref.data");
        output_vx.to_file(vx, "vx_ref.data");
        output_vy.to_file(vy, "vy_ref.data");
    }


    void compute() {
        auto fx = [&](double x, double y) {
            return
            (12 - 24 * y) * x*x*x*x +
            (-24 + 48*y) * x*x*x +
            (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
            (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
            1 - 4 * y + 12 * y*y - 8 * y*y*y;
        };
        auto fy = [&](double x, double y) {
            return
            (8 - 48 * y + 48 * y*y) * x*x*x +
            (-12 + 72 * y - 72 * y*y) * x*x +
            (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x -
            12 * y*y + 24 * y*y*y - 12 * y*y*y*y;
        };

        std::vector<double> rhs(now.vx.size() + now.vy.size() + now.p.size());

        vector_view vx{ rhs.data(), now.vx.sizes() };
        vector_view vy{ vx.data() + vx.size(), now.vy.sizes() };
        vector_view p{ vy.data() + vy.size(), now.p.sizes() };

        compute_projection(vx, x_1.basis, y.basis, fx);
        compute_projection(vy, x.basis, y_1.basis, fy);

        mumps::problem problem(rhs.data(), rhs.size());
        double mu = 1;

        auto nvx = vx.data() - rhs.data();
        auto nvy = vy.data() - rhs.data();
        auto np = p.data() - rhs.data();

        // 11 - vx, vx
        tensor_product(vx, vx, problem, Kx_1, y.M, 2 * mu);
        tensor_product(vx, vx, problem, x_1.M, Ky, mu);

        // 12 - vx, vy
        tensor_product(vx, vy, problem, A12x, A21y, mu, nvx, nvy);

        // 13 - vx, p
        tensor_product(vx, p, problem, A11x, y.M, -1, nvx, np);


        // 21 - vy, vx
        tensor_product(vy, vx, problem, A21x, A12y, mu, nvy, nvx);

        // 22 - vy, vy
        tensor_product(vy, vy, problem, x.M, Ky_1, 2 * mu, nvy, nvy);
        tensor_product(vy, vy, problem, Kx, x_1.M, mu, nvy, nvy);

        // 23 - vy, p
        tensor_product(vy, p, problem, x.M, A11y, -1, nvy, np);


        // 31 - p, vx
        tensor_product(p, vx, problem, A22x, y.M, 1, np, nvx);

        // 32 - p, vy
        tensor_product(p, vy, problem, x.M, A22y, 1, np, nvy);



        solver.solve(problem);

        output_p.to_file(p, "pressure.data");
        output_vx.to_file(vx, "vx.data");
        output_vy.to_file(vy, "vy.data");
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs_1();
        solve(u, Kx, y.M);

        using std::swap;
        swap(u, u_prev);

        compute_rhs_2();
        solve(u, x.M, Ky);
    }

    void after_step(int iter, double /*t*/) override {
        if ((iter + 1) % 1 == 0) {
            double E = energy();
            double L2 = L2norm();
            std::cout << (iter + 1) * steps.dt << " " << E << " " << L2 << std::endl;
            // output.to_file(u, "out_%d.data", (iter + 1) / 1);
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

    double energy() const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                value_type vx = eval_fun(u, e, q);
                E += vx.val * w * J;
            }
        }
        return E;
    }

    double L2norm() const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                value_type vx = eval_fun(u, e, q);
                E += vx.val * vx.val * w * J;
            }
        }
        return E;
    }
};

} // problems
} // ads

#endif /* ADS_PROBLEMS_STOKES_STOKES_HPP */
