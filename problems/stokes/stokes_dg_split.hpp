#ifndef PROBLEMS_STOKES_STOKES_DG_SPLIT_HPP_
#define PROBLEMS_STOKES_STOKES_DG_SPLIT_HPP_

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/solver/mumps.hpp"
#include "space_set.hpp"


namespace ads {

enum class scheme { BE, CN, peaceman_rachford, strang_BE, strang_CN };

class stokes_dg_split : public simulation_2d {
private:
    using Base = simulation_2d;
    using value_pair = std::array<value_type, 2>;

    galois_executor executor{8};

    space_set trial, test;
    scheme method;

    vector_type vx, vy, p;

    double h;
    double eta = 100;


    mumps::solver solver;
    output_manager<2> outputU1, outputU2, outputP;

public:
    stokes_dg_split(scheme method, space_set trial_, space_set test_, const timesteps_config& steps)
    : Base{ test_.Px, test_.Py, steps }
    , trial{ std::move(trial_) }
    , test{ std::move(test_) }
    , method{ method }
    , vx{{ trial.U1x.dofs(), trial.U1y.dofs() }}
    , vy{{ trial.U2x.dofs(), trial.U2y.dofs() }}
    , p{{ trial.Px.dofs(), trial.Py.dofs() }}
    , h{ element_diam(trial.Px, trial.Py) }
    , outputU1{ trial.U1x.B, trial.U1y.B, 500 }
    , outputU2{ trial.U2x.B, trial.U2y.B, 500 }
    , outputP{ trial.Px.B, trial.Py.B, 500 }
    { }

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

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
        project(p, trial.Px, trial.Py, [this](point_type x) { return exact_p(x, 0).val; });

        // save_to_file(0);
        // output_exact(0);

        // std::cout << "Initial state error:" << std::endl;
        // print_error(0);
    }

    value_type exact_p(point_type p, double t) const {
        auto x = p[0];
        auto et = std::exp(-t);
        return {et * x * (1 - x), et * (1 - 2 * x), 0.0};
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


    void assemble_matrix(mumps::problem& problem,
                         double cx, double cy,
                         double dt,
                         double sx, double sy) const {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();

        auto D = DU1 + DU2 + DP;

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DU1, DU1, problem);
        auto test_p = shifted(DU1 + DU2, DU1 + DU2, problem);

        // auto test_vxy = shifted(0, DU1, problem);
        // auto test_vyx = shifted(DU1, 0, problem);

        auto trial_vx = shifted(D, D, problem);
        auto trial_vy = shifted(D + dU1, D + dU1, problem);
        auto trial_p = shifted(D + dU1 + dU2, D + dU1 + dU2, problem);

        // Gram matrix
        // G(w, u)
        // w = (tx, ty, w)
        // u = (vx, vy, p)

        // w, p -> (w, p) +  hh ([w], [p])
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : touching_dofs(i, test.Px, test.Py)) {
                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, test.Px, test.Py) + 1;

                if (! is_pressure_fixed(i) && ! is_pressure_fixed(j)) {
                    auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, test.Px, test.Py, form); };
                    double product = eval([](auto w, auto p) { return w.val * p.val; });

                    auto form = [this](auto w, auto p, auto, auto) {
                        return h * jump(w).val * jump(p).val;
                    };
                    double boundary = integrate_over_internal_skeleton(i, j, test.Px, test.Py, test.Px, test.Py, form);
                    test_p(ii, jj, product + boundary);
                }
            }
        }

        // tx, vx -> (\/tx, \/vx) + 1/h ([tx], [vx])
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : touching_dofs(i, test.U1x, test.U1y)) {
                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, test.U1x, test.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form); };

                // Weak BC
                // if (! is_boundary(i[0], test.U1x) && ! is_boundary(j[0], test.U1x)) {
                // Strong BC
                if (! is_boundary(i, test.U1x, test.U1y) && ! is_boundary(j, test.U1x, test.U1y)) {
                    double val = eval([sx,sy](auto tx, auto vx) {
                        return sx * tx.dx * vx.dx + sy * tx.dy * vx.dy;
                    });
                    // skeleton
                    auto form = [this](auto tx, auto vx, auto, auto) {
                        return 1/h * jump(tx).val * jump(vx).val;
                    };
                    val += integrate_over_skeleton(i, j, test.U1x, test.U1y, test.U1x, test.U1y, form);

                    test_vx(ii, jj, val);
                }
            }
        }

        // ty, vy -> (\/ty, \/vy) + 1/h ([ty], [vy])
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : touching_dofs(i, test.U2x, test.U2y)) {
                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, test.U2x, test.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form); };

                // Weak BC
                // if (! is_boundary(i[1], test.U2y) && ! is_boundary(j[1], test.U2y)) {
                // Strong BC
                if (! is_boundary(i, test.U2x, test.U2y) && ! is_boundary(j, test.U2x, test.U2y)) {
                    double val = eval([sx,sy](auto ty, auto vy) {
                        return sx * ty.dx * vy.dx + sy * ty.dy * vy.dy;
                    });
                    // skeleton
                    auto form = [this](auto ty, auto vy, auto, auto) {
                        return 1/h * jump(ty).val * jump(vy).val;
                    };
                    val += integrate_over_skeleton(i, j, test.U2x, test.U2y, test.U2x, test.U2y, form);

                    test_vy(ii, jj, val);
                }
            }
        }

        // Dirichlet BC - test space
        // Weak BC
        // for (auto iy = 0; iy < test.U1y.dofs(); ++ iy) {
        //     int i0 = linear_index({0, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i0, i0, 1);
        //     int i1 = linear_index({test.U1x.dofs() - 1, iy}, test.U1x, test.U1y) + 1;
        //     test_vx(i1, i1, 1);
        // }
        // for (auto ix = 0; ix < test.U2x.dofs(); ++ ix) {
        //     int i0 = linear_index({ix, 0}, test.U2x, test.U2y) + 1;
        //     test_vy(i0, i0, 1);
        //     int i1 = linear_index({ix, test.U2y.dofs() - 1}, test.U2x, test.U2y) + 1;
        //     test_vy(i1, i1, 1);
        // }
        // Strong BC
        for_boundary_dofs(test.U1x, test.U1y, [&](index_type dof) {
            int i = linear_index(dof, test.U1x, test.U1y) + 1;
            test_vx(i, i, 1);
        });
        for_boundary_dofs(test.U2x, test.U2y, [&](index_type dof) {
            int i = linear_index(dof, test.U2x, test.U2y) + 1;
            test_vy(i, i, 1);
        });


        int i = linear_index({0, 0}, test.Px, test.Py) + 1;
        test_p(i, i, 1.0);


        // B, B^t
        auto put = [&](int i, int j, int si, int sj, double val, bool fixed_i, bool fixed_j) {
            int ii = i + si;
            int jj = j + sj;

            if (!fixed_i) {
                problem.add(ii, D + jj, -val);
            }
            if (!fixed_i && !fixed_j) {
                problem.add(D + jj, ii, val);
            }
        };

        // B(u, v)
        // u = (ux, uy, p)  trial
        // v = (vx, vy, q)  test

        // auto pen_term = [this](auto w, auto u) { return - 2 * Cpen / h * w.val * u.val; };

        // ux, vx -> (\/ux, \/vx) - ({\/ux}, [vx])
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! dofs_touch(i, test.U1x, test.U1y, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[0], test.U1x);
                // bool bd_j = is_boundary(j[0], trial.U1x);
                // Strong BC
                bool bd_i = is_boundary(i, test.U1x, test.U1y);
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y);

                double value = eval([cx,cy](auto vx, auto ux) {
                    return vx.val * ux.val + cx * vx.dx * ux.dx + cy * vx.dy * ux.dy;
                });

                // skeleton
                auto form = [this](auto vx, auto ux, auto, auto n) {
                    return eta / h * jump(ux).val * jump(vx).val
                         - dot(average(ux), n) * jump(vx).val
                         - dot(average(vx), n) * jump(ux).val;
                };
                value += integrate_over_skeleton(i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, form);

                // boundary terms
                // top/bottom
                // if (is_boundary(i[1], test.U1y) || is_boundary(j[1], trial.U1y)) {
                //     auto side = i[1] == 0 ? boundary::bottom : boundary::top;
                //     // int n = side == boundary::top ? 1 : -1;
                //     value -= integrate_boundary(side, i, j, test.U1x, test.U1y, trial.U1x, trial.U1y, [&](auto w, auto u, auto) {
                //         return pen_term(w, u);
                //     });
                // }
                put(ii, jj, 0, 0, value, bd_i, bd_j);
            }
        }

        // uy, vy -> (\/uy, \/vy) - ({\/uy}, [vy])
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! dofs_touch(i, test.U2x, test.U2y, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[1], test.U2y);
                // bool bd_j = is_boundary(j[1], trial.U2y);
                // Strong BC
                bool bd_i = is_boundary(i, test.U2x, test.U2y);
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y);

                double value = eval([cx,cy](auto vy, auto uy) {
                    return vy.val * uy.val + cx * vy.dx * uy.dx + cy * vy.dy * uy.dy;
                });

                // skeleton
                auto form = [this](auto vy, auto uy, auto, auto n) {
                    return eta / h * jump(uy).val * jump(vy).val
                         - dot(average(uy), n) * jump(vy).val
                         - dot(average(vy), n) * jump(uy).val;
                };
                value += integrate_over_skeleton(i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, form);

                // boundary terms
                // left/right
                // if (is_boundary(i[0], test.U2x) || is_boundary(j[0], trial.U2x)) {
                //     auto side = i[0] == 0 ? boundary::left : boundary::right;
                //     // int n = side == boundary::right ? 1 : -1;
                //     value -= integrate_boundary(side, i, j, test.U2x, test.U2y, trial.U2x, trial.U2y, [&](auto w, auto u, auto) {
                //         return pen_term(w, u);
                //     });
                // }
                put(ii, jj, DU1, dU1, value, bd_i, bd_j);
            }
        }

        // ux, q -> (ux,x, q)
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U1x, trial.U1y)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.U1x, trial.U1y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U1x, trial.U1y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U1x, trial.U1y, form); };

                bool fixed_i = is_pressure_fixed(i);
                // Weak BC
                // bool bd_j = is_boundary(j[0], trial.U1x);
                // Strong BC
                bool bd_j = is_boundary(j, trial.U1x, trial.U1y);

                double val = eval([](auto q, auto ux) { return q.val * ux.dx; });

                // skeleton
                auto form = [this](auto q, auto ux, auto, auto n) {
                    return - h * average(q).val * n[0] * jump(ux).val;
                };
                val += integrate_over_skeleton(i, j, test.Px, test.Py, trial.U1x, trial.U1y, form);

                put(ii, jj, DU1 + DU2, 0, val, fixed_i, bd_j);
            }
        }

        // uy, q -> (uy,y, q),
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.U2x, trial.U2y)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.U2x, trial.U2y)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.U2x, trial.U2y) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.Px, test.Py, trial.U2x, trial.U2y, form); };

                bool fixed_i = is_pressure_fixed(i);
                // Weak BC
                // bool bd_j = is_boundary(j[1], trial.U2y);
                // Strong BC
                bool bd_j = is_boundary(j, trial.U2x, trial.U2y);

                double val = eval([](auto q, auto uy) { return q.val * uy.dy; });

                // skeleton
                auto form = [this](auto q, auto uy, auto, auto n) {
                    return - h * average(q).val * n[1] * jump(uy).val;
                };
                val += integrate_over_skeleton(i, j, test.Px, test.Py, trial.U2x, trial.U2y, form);

                put(ii, jj, DU1 + DU2, dU1, val, fixed_i, bd_j);
            }
        }

        // p, vx ->  - (p, vx,x) + h [vx] n1 {p}
        for (auto i : dofs(test.U1x, test.U1y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.U1x, test.U1y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U1x, test.U1y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U1x, test.U1y, trial.Px, trial.Py, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[0], test.U1x);
                // Strong BC
                bool bd_i = is_boundary(i, test.U1x, test.U1y);

                bool fixed_j = is_pressure_fixed(j);

                // double val = eval([cx](auto vx, auto p) { return - cx * vx.dx * p.val; });
                double val = eval([dt](auto vx, auto p) { return - dt * vx.dx * p.val; });

                // skeleton
                auto form = [this](auto vx, auto p, auto, auto n) {
                    return h * average(p).val * n[0] * jump(vx).val;
                };
                val += integrate_over_skeleton(i, j, test.U1x, test.U1y, trial.Px, trial.Py, form);

                put(ii, jj, 0, dU1 + dU2, val, bd_i, fixed_j);
            }
        }

        // p, vy ->  - (p, vy,y) + h [vy] n2 {p}
        for (auto i : dofs(test.U2x, test.U2y)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.U2x, test.U2y, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.U2x, test.U2y) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;
                auto eval = [&](auto form) { return integrate(i, j, test.U2x, test.U2y, trial.Px, trial.Py, form); };

                // Weak BC
                // bool bd_i = is_boundary(i[1], test.U2y);
                // Strong BC
                bool bd_i = is_boundary(i, test.U2x, test.U2y);

                bool fixed_j = is_pressure_fixed(j);

                // double val = eval([cy](auto vy, auto p) { return - cy * vy.dy * p.val; });
                double val = eval([dt](auto vy, auto p) { return - dt * vy.dy * p.val; });

                // skeleton
                auto form = [this](auto vy, auto p, auto, auto n) {
                    return h * average(p).val * n[1] * jump(vy).val;
                };
                val += integrate_over_skeleton(i, j, test.U2x, test.U2y, trial.Px, trial.Py, form);

                put(ii, jj, DU1, dU1 + dU2, val, bd_i, fixed_j);
            }
        }

        // Stabilization term
        // S(p, q) = h [p][q]
        for (auto i : dofs(test.Px, test.Py)) {
            for (auto j : dofs(trial.Px, trial.Py)) {
                if (! dofs_touch(i, test.Px, test.Py, j, trial.Px, trial.Py)) continue;

                int ii = linear_index(i, test.Px, test.Py) + 1;
                int jj = linear_index(j, trial.Px, trial.Py) + 1;

                bool fixed_i = is_pressure_fixed(i);
                bool fixed_j = is_pressure_fixed(j);

                // skeleton
                auto form = [this](auto q, auto p, auto, auto) {
                    return h * jump(p).val * jump(q).val;
                };
                double val = integrate_over_internal_skeleton(i, j, test.Px, test.Py, trial.Px, trial.Py, form);

                put(ii, jj, DU1 + DU2, dU1 + dU2, val, fixed_i, fixed_j);
            }
        }

        // Dirichlet BC - trial space
        // Weak BC
        // for (auto iy = 0; iy < trial.U1y.dofs(); ++ iy) {
        //     int i0 = linear_index({0, iy}, trial.U1x, trial.U1y) + 1;
        //     trial_vx(i0, i0, 1);
        //     int i1 = linear_index({trial.U1x.dofs() - 1, iy}, trial.U1x, trial.U1y) + 1;
        //     trial_vx(i1, i1, 1);
        // }
        // for (auto ix = 0; ix < trial.U2x.dofs(); ++ ix) {
        //     int i0 = linear_index({ix, 0}, trial.U2x, trial.U2y) + 1;
        //     trial_vy(i0, i0, 1);
        //     int i1 = linear_index({ix, trial.U2y.dofs() - 1}, trial.U2x, trial.U2y) + 1;
        //     trial_vy(i1, i1, 1);
        // }
        // Strong BC
        for_boundary_dofs(trial.U1x, trial.U1y, [&](index_type dof) {
            int i = linear_index(dof, trial.U1x, trial.U1y) + 1;
            trial_vx(i, i, 1);
        });
        for_boundary_dofs(trial.U2x, trial.U2y, [&](index_type dof) {
            int i = linear_index(dof, trial.U2x, trial.U2y) + 1;
            trial_vy(i, i, 1);
        });

        int ii = linear_index({0, 0}, trial.Px, trial.Py) + 1;
        trial_p(ii, ii, 1.0);
    }

    template <typename Fun>
    void compute_rhs(double cx, double cy,
                     vector_view& Rvx, vector_view& Rvy, vector_view& /*Rp*/,
                     double dt, Fun&& forcing) const {
        using shape = std::array<std::size_t, 2>;
        auto u1_shape = shape{ test.U1x.basis.dofs_per_element(), test.U1y.basis.dofs_per_element() };
        auto u2_shape = shape{ test.U2x.basis.dofs_per_element(), test.U2y.basis.dofs_per_element() };


        executor.for_each(elements(test.Px, test.Py), [&](index_type e) {
            auto vx_loc = vector_type{ u1_shape };
            auto vy_loc = vector_type{ u2_shape };

            double J = jacobian(e);
            for (auto q : quad_points(test.Px, test.Py)) {
                double W = weigth(q);
                auto x = point(e, q);
                auto F = forcing(x);
                value_type vvx = eval(vx, e, q, trial.U1x, trial.U1y);
                value_type vvy = eval(vy, e, q, trial.U2x, trial.U2y);
                // value_type pp = eval(p, e, q, trial.Px, trial.Py);

                for (auto a : dofs_on_element(e, test.U1x, test.U1y)) {
                    auto aa = dof_global_to_local(e, a, test.U1x, test.U1y);
                    value_type v = eval_basis(e, q, a, test.U1x, test.U1y);

                    double M = vvx.val * v.val;
                    // double Lx = vvx.dx * v.dx - pp.val * v.dx;
                    double Lx = vvx.dx * v.dx;
                    double Ly = vvx.dy * v.dy;
                    double val = M + cx * Lx + cy * Ly + dt * F[0] * v.val;
                    vx_loc(aa[0], aa[1]) -= val * W * J;
                }
                for (auto a : dofs_on_element(e, test.U2x, test.U2y)) {
                    auto aa = dof_global_to_local(e, a, test.U2x, test.U2y);
                    value_type v = eval_basis(e, q, a, test.U2x, test.U2y);

                    double M = vvy.val * v.val;
                    double Lx = vvy.dx * v.dx;
                    // double Ly = vvy.dy * v.dy - pp.val * v.dy;
                    double Ly = vvy.dy * v.dy;
                    double val = M + cx * Lx + cy * Ly + dt * F[1] * v.val;
                    vy_loc(aa[0], aa[1]) -= val * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(Rvx, vx_loc, e, test.U1x, test.U1y);
                update_global_rhs(Rvy, vy_loc, e, test.U2x, test.U2y);
            });
        });
    }


    void apply_bc(vector_view& Rvx, vector_view& Rvy, vector_view& Rp) {
        // Strong BC
        zero_bc(Rvx, test.U1x, test.U1y);
        zero_bc(Rvy, test.U2x, test.U2y);

        // Weak BC
        // vx = 0 on left/right edge
        // for (auto i = 0; i < test.U1y.dofs(); ++ i) {
        //     Rvx(0, i) = 0;
        //     Rvx(test.U1x.dofs() - 1, i) = 0;
        // }

        // vy = 0 on top/bottom edge
        // for (auto i = 0; i < test.U2x.dofs(); ++ i) {
        //     Rvy(i, 0) = 0;
        //     Rvy(i, test.U2y.dofs() - 1) = 0;
        // }

        // zero_bc(Rp, test.Px, test.Py);
        int i = linear_index({0, 0}, test.Px, test.Py);
        Rp(i, i) = 0; // fix pressure at a point
    }

    void step(int /*iter*/, double t) override {
        using namespace std::placeholders;
        auto dt = steps.dt;

        auto f = [&](point_type x, double s) { return forcing(x, s); };
        auto F = [&](double s) { return std::bind(f, _1, s); };
        auto Favg = [&](double s1, double s2) {
            return [=,&f](point_type x) {
                auto f1 = f(x, s1);
                auto f2 = f(x, s2);
                return point_type{0.5 * (f1[0] + f2[0]), 0.5 * (f1[1] + f2[1])};
            };
        };
        auto zero = [&](point_type) { return point_type{0, 0}; };

        if (method == scheme::BE) {
            substep(1, 1, dt, dt, 0, 0, dt, F(t + dt));
        }
        if (method == scheme::CN) {
            substep(1, 1, dt/2, dt/2, -dt/2, -dt/2,   dt, Favg(t, t + dt));
        }
        if (method == scheme::peaceman_rachford) {
            substep(1, 1, dt/2,    0,     0, -dt/2,   dt/2, F(t + dt/2));
            substep(1, 1,    0, dt/2, -dt/2,     0,   dt/2, F(t + dt/2));
        }
        if (method == scheme::strang_BE) {
            substep(1, 1, dt/2,  0,   0, 0,   dt/2, F(t + dt/2));
            substep(1, 1,    0, dt,   0, 0,     dt, zero);
            substep(1, 1, dt/2,  0,   0, 0,   dt/2, F(t + dt));
        }
        if (method == scheme::strang_CN) {
            substep(1, 1, dt/4,    0,   -dt/4,     0,   dt/2, Favg(t, t + dt/2));
            substep(1, 1,    0, dt/2,       0, -dt/2,     dt, zero);
            substep(1, 1, dt/4,    0,   -dt/4,     0,   dt/2, Favg(t + dt/2, t + dt));
        }
    }

    void copy_solution(const vector_view& vx_new, const vector_view& vy_new, const vector_view& p_new) {
        for (auto i : dofs(trial.U1x, trial.U1y)) {
            vx(i[0], i[1]) = vx_new(i[0], i[1]);
        }
        for (auto i : dofs(trial.U2x, trial.U2y)) {
            vy(i[0], i[1]) = vy_new(i[0], i[1]);
        }
        for (auto i : dofs(trial.Px, trial.Py)) {
            p(i[0], i[1]) = p_new(i[0], i[1]);
        }
    }

    template <typename Fun>
    void substep(double sx, double sy,
                 double Lx_lhs, double Ly_lhs,
                 double Lx_rhs, double Ly_rhs,
                 double dt, Fun&& f) {
        auto dU1 = trial.U1x.dofs() * trial.U1y.dofs();
        auto dU2 = trial.U2x.dofs() * trial.U2y.dofs();
        auto dP = trial.Px.dofs() * trial.Py.dofs();
        auto dim_trial = dU1 + dU2 + dP;

        auto DU1 = test.U1x.dofs() * test.U1y.dofs();
        auto DU2 = test.U2x.dofs() * test.U2y.dofs();
        auto DP = test.Px.dofs() * test.Py.dofs();
        auto dim_test = DU1 + DU2 + DP;

        std::vector<double> rhs(dim_test + dim_trial);

        vector_view Rvx{rhs.data(), {test.U1x.dofs(), test.U1y.dofs()}};
        vector_view Rvy{Rvx.data() + DU1, {test.U2x.dofs(), test.U2y.dofs()}};
        vector_view Rp{Rvy.data() + DU2, {test.Px.dofs(), test.Py.dofs()}};

        vector_view vx{Rp.data() + DP, {trial.U1x.dofs(), trial.U1y.dofs()}};
        vector_view vy{vx.data() + dU1, {trial.U2x.dofs(), trial.U2y.dofs()}};
        vector_view p{vy.data() + dU2, {trial.Px.dofs(), trial.Py.dofs()}};


        mumps::problem problem(rhs.data(), rhs.size());

        // std::cout << "Assembling matrix" << std::endl;
        assemble_matrix(problem, Lx_lhs, Ly_lhs, dt, sx, sy);

        // std::cout << "Computing RHS" << std::endl;
        compute_rhs(Lx_rhs, Ly_rhs, Rvx, Rvy, Rp, dt, std::forward<Fun>(f));

        apply_bc(Rvx, Rvy, Rp);

        // std::cout << "Solving" << std::endl;
        solver.solve(problem);

        copy_solution(vx, vy, p);
    }

    void after_step(int iter, double t) override {
        int i = iter + 1;
        double tt = t + steps.dt;

        if (i % 1 == 0) {
            // std::cout << "Outputting" << std::endl;
            // save_to_file(i);
        }

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

        std::cout << i << " " << tt << " " << vL2 << " " << vH1 << " " << pL2 << " " << pH1 << std::endl;
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

    bool dofs_touch(int a, const dimension& U, int b, const dimension& V) const {
        auto ar = U.basis.element_ranges[a];
        auto br = V.basis.element_ranges[b];
        return (ar.first >= br.first && ar.first <= br.second + 1) || (br.first >= ar.first && br.first <= ar.second + 1);
    }

    bool dofs_touch(index_type a, const dimension& Ux, const dimension& Uy,
                    index_type b, const dimension& Vx, const dimension& Vy) const {
        return dofs_touch(a[0], Ux, b[0], Vx) && dofs_touch(a[1], Uy, b[1], Vy);
    }


    double dot(value_type a, point_type b) const {
        return a.dx * b[0] + a.dy * b[1];
    }

    value_type average(const value_pair& v) const {
        return 0.5 * (v[0] + v[1]);
    }

    value_type jump(const value_pair& v) const {
        return v[0] - v[1];
    }

    value_type eval_basis_at(point_type p, index_type span, index_type dof,
                             const dimension& x, const dimension& y) const {
        bspline::eval_ders_ctx cx{x.p, 1};
        bspline::eval_ders_ctx cy{y.p, 1};

        double** bvx = cx.basis_vals();
        double** bvy = cy.basis_vals();

        eval_basis_with_derivatives(span[0], p[0], x.B, bvx, 1, cx);
        eval_basis_with_derivatives(span[1], p[1], y.B, bvy, 1, cy);

        int offsetx = span[0] - x.p;
        int offsety = span[1] - y.p;

        int ix = dof[0] - offsetx;
        int iy = dof[1] - offsety;

        if (ix < 0 || ix > x.p || iy < 0 || iy > y.p) return {};

        auto value = bvx[0][ix] * bvy[0][iy];
        auto dx    = bvx[1][ix] * bvy[0][iy];
        auto dy    = bvx[0][ix] * bvy[1][iy];

        return { value, dx, dy };
    }

    value_pair eval_basis_at_edge(point_type p, boundary orientation, index_type dof,
                                  const dimension& x, const dimension& y) const {
        int spanx = bspline::find_span(p[0], x.B);
        int spany = bspline::find_span(p[1], y.B);
        auto span = index_type{spanx, spany};

        auto val1 = eval_basis_at(p, span, dof, x, y);
        auto val0 = value_type{};

        if (orientation == boundary::vertical) {
            int spanx0 = bspline::find_span(p[0] - 1e-10, x.B);
            val0 = eval_basis_at(p, index_type{spanx0, spany}, dof, x, y);
        } else {
            int spany0 = bspline::find_span(p[1] - 1e-10, y.B);
            val0 = eval_basis_at(p, index_type{spanx, spany0}, dof, x, y);
        }
        return {val0, val1};
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

    auto touching_dofs(int dof, int begin, int end, const dimension& x) const {
        using std::min;
        using std::max;

        auto elems = x.basis.element_ranges[dof];
        auto first = x.basis.first_dof(max(elems.first - 1, 0));
        auto last = x.basis.last_dof(min(elems.second + 1, x.elements - 1));

        auto minx = max(begin, first);
        auto maxx = min(end, last);

        return boost::counting_range(minx, maxx + 1);
    }

    index_range touching_dofs(index_type dof, const dimension& x, const dimension& y) const {
        auto rx = touching_dofs(dof[0], 0, x.dofs(), x);
        auto ry = touching_dofs(dof[1], 0, y.dofs(), y);
        return util::product_range<index_type>(rx, ry);
    }


    template <typename Form>
    double integrate_vface(int sx, int ey, index_type i, index_type j,
                           const dimension& Ux, const dimension& Uy,
                           const dimension& Vx, const dimension& Vy,
                           Form&& form) const {

        if (! touches_point(i[0], sx, Ux)   || ! touches_point(j[0], sx, Vx))   return 0;
        if (! supported_in_1d(i[1], ey, Uy) || ! supported_in_1d(j[1], ey, Vy)) return 0;

        double x0 = Ux.basis.points[sx];
        auto n = point_type{1, 0};

        double val = 0;
        for (int q = 0; q < Uy.basis.quad_order; ++ q) {
            double w = Uy.basis.w[q];
            double J = Uy.basis.J[ey];

            auto x = point_type{x0, Uy.basis.x[ey][q]};
            auto uu = eval_basis_at_edge(x, boundary::vertical, i, Ux, Uy);
            auto vv = eval_basis_at_edge(x, boundary::vertical, j, Vx, Vy);

            auto fuv = form(uu, vv, x, n);
            val += fuv * w * J;
        }
        return val;
    }

    template <typename Form>
    double integrate_hface(int sy, int ex, index_type i, index_type j,
                           const dimension& Ux, const dimension& Uy,
                           const dimension& Vx, const dimension& Vy,
                           Form&& form) const {

        if (! supported_in_1d(i[0], ex, Ux) || ! supported_in_1d(j[0], ex, Vx)) return 0;
        if (! touches_point(i[1], sy, Uy)   || ! touches_point(j[1], sy, Vy))   return 0;

        double y0 = Uy.basis.points[sy];
        auto n = point_type{0, 1};

        double val = 0;
        for (int q = 0; q < Ux.basis.quad_order; ++ q) {
            double w = Ux.basis.w[q];
            double J = Ux.basis.J[ex];

            auto x = point_type{Ux.basis.x[ex][q], y0};
            auto uu = eval_basis_at_edge(x, boundary::horizontal, i, Ux, Uy);
            auto vv = eval_basis_at_edge(x, boundary::horizontal, j, Vx, Vy);

            auto fuv = form(uu, vv, x, n);
            val += fuv * w * J;
        }
        return val;
    }

    template <typename Form>
    double integrate_over_skeleton(index_type i, index_type j,
                                   const dimension& Ux, const dimension& Uy,
                                   const dimension& Vx, const dimension& Vy,
                                   Form&& form) const {
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return 0;

        // TODO: restrict the range
        auto rUx = Ux.basis.element_ranges[i[0]];
        auto rVx = Vx.basis.element_ranges[j[0]];
        auto rUy = Uy.basis.element_ranges[i[1]];
        auto rVy = Vy.basis.element_ranges[j[1]];
        auto ex0 = std::min(rUx.first, rVx.first);
        auto ex1 = std::max(rUx.second, rVx.second);
        auto ey0 = std::min(rUy.first, rVy.first);
        auto ey1 = std::max(rUy.second, rVy.second);

        double val = 0;
        // for (int ix = 0; ix <= Ux.elements; ++ ix) {
        for (int ix = ex0; ix <= ex1 + 1; ++ ix) {
            // for (auto ey : Uy.element_indices()) {
            for (auto ey = ey0; ey <= ey1; ++ ey) {
                val += integrate_vface(ix, ey, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        for (auto ex = ex0; ex <= ex1; ++ ex) {
        // for (auto ex : Ux.element_indices()) {
            // for (int iy = 0; iy <= Uy.elements; ++ iy) {
            for (int iy = ey0; iy <= ey1 + 1; ++ iy) {
                val += integrate_hface(iy, ex, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        return val;
    }

    template <typename Form>
    double integrate_over_internal_skeleton(index_type i, index_type j,
                                            const dimension& Ux, const dimension& Uy,
                                            const dimension& Vx, const dimension& Vy,
                                            Form&& form) const {
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return 0;


        double val = 0;
        // TODO: restrict the range
        for (int ix = 1; ix < Ux.elements; ++ ix) {
            for (auto ey : Uy.element_indices()) {
                val += integrate_vface(ix, ey, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        for (auto ex : Ux.element_indices()) {
            for (int iy = 1; iy < Uy.elements; ++ iy) {
                val += integrate_hface(iy, ex, i, j, Ux, Uy, Vx, Vy, form);
            }
        }
        return val;
    }


    template <typename Form>
    double integrate_boundary(boundary side, index_type i, index_type j, const dimension& Ux, const dimension& Uy,
                              const dimension& Vx, const dimension& Vy, Form&& form) const {
        double val = 0;
        bool horizontal = side == boundary::top || side == boundary::bottom;

        if (horizontal) {
            int ey = side == boundary::bottom ? 0 : Uy.elements - 1;
            if (! supported_in_1d(j[1], ey, Vy) || ! supported_in_1d(i[1], ey, Uy)) return 0;

            auto y0 = side == boundary::bottom ? Uy.a : Uy.b;

            for (auto e : Ux.basis.element_range(i[0])) {
                if (! supported_in_1d(j[0], e, Vx)) continue;

                double J = Ux.basis.J[e];

                for (int q = 0; q < Ux.basis.quad_order; ++ q) {
                    double w = Ux.basis.w[q];
                    point_type x{Ux.basis.x[e][q], y0};
                    value_type ww = eval_basis_at(x, {e, ey}, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, {e, ey}, j, Vx, Vy);
                    double fuw = form(ww, uu, x);
                    val += fuw * w * J;
                }
            }
        } else {
            int ex = side == boundary::left ? 0 : Ux.elements - 1;
            if (! supported_in_1d(j[0], ex, Vx) || ! supported_in_1d(i[0], ex, Ux)) return 0;

            auto x0 = side == boundary::left ? Ux.a : Ux.b;

            for (auto e : Uy.basis.element_range(i[1])) {
                if (! supported_in_1d(j[1], e, Vy)) continue;

                double J = Uy.basis.J[e];

                for (int q = 0; q < Uy.basis.quad_order; ++ q) {
                    double w = Uy.basis.w[q];
                    point_type x{x0, Uy.basis.x[e][q]};
                    value_type ww = eval_basis_at(x, {ex, e}, i, Ux, Uy);
                    value_type uu = eval_basis_at(x, {ex, e}, j, Vx, Vy);
                    double fuw = form(ww, uu, x);
                    val += fuw * w * J;
                }
            }
        }
        return val;
    }
};

}

#endif // PROBLEMS_STOKES_STOKES_DG_SPLIT_HPP_
