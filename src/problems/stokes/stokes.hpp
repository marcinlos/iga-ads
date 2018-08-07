#ifndef PROBLEMS_STOKES_STOKES_HPP_
#define PROBLEMS_STOKES_STOKES_HPP_

#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "mumps.hpp"


namespace ads {

class stokes : public simulation_2d {
private:
    using Base = simulation_2d;

    galois_executor executor{8};

    dimension Ux, Uy;
    dimension& Vx;
    dimension& Vy;

    lin::band_matrix MVx, MVy;
    lin::band_matrix KVx, KVy;
    lin::band_matrix AVx, AVy;

    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;

    double h;
    vector_type buffer;

    mumps::solver solver;
    output_manager<2> output;

public:
    stokes(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps)
    : Base{ std::move(test_x), std::move(test_y), steps }
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , Vx{ x }
    , Vy{ y }
    , MVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , MVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , KVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , KVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , AVx{Vx.p, Vx.p, Vx.dofs(), Vx.dofs(), 0}
    , AVy{Vy.p, Vy.p, Vy.dofs(), Vy.dofs(), 0}
    , MUVx{ Vx.dofs(), Ux.dofs() }
    , MUVy{ Vy.dofs(), Uy.dofs() }
    , KUVx{ Vx.dofs(), Ux.dofs() }
    , KUVy{ Vy.dofs(), Uy.dofs() }
    , AUVx{ Vx.dofs(), Ux.dofs() }
    , AUVy{ Vy.dofs(), Uy.dofs() }
    , h{ element_diam(Ux, Uy) }
    , buffer{{ Ux.dofs(), Uy.dofs() }}
    , output{ Ux.B, Uy.B, 500 }
    { }

    double element_diam(const dimension& Ux, const dimension& Uy) const {
        return std::sqrt(max_element_size(Ux) * max_element_size(Uy));
    }

    void prepare_matrices() {
        gram_matrix_1d(MVx, Vx.basis);
        gram_matrix_1d(MVy, Vy.basis);
        gram_matrix_1d(MUVx, Ux.basis, Vx.basis);
        gram_matrix_1d(MUVy, Uy.basis, Vy.basis);

        stiffness_matrix_1d(KVx, Vx.basis);
        stiffness_matrix_1d(KVy, Vy.basis);
        stiffness_matrix_1d(KUVx, Ux.basis, Vx.basis);
        stiffness_matrix_1d(KUVy, Uy.basis, Vy.basis);

        advection_matrix_1d(AVx, Vx.basis);
        advection_matrix_1d(AVy, Vy.basis);
        advection_matrix_1d(AUVx, Ux.basis, Vx.basis);
        advection_matrix_1d(AUVy, Uy.basis, Vy.basis);
    }

    void before() override {
        prepare_matrices();

        Ux.factorize_matrix();
        Uy.factorize_matrix();

        output_exact();
    }

    double exact_p(point_type p) const {
        auto x = p[0];
        return x * (1 - x);
    }

    point_type exact_v(point_type p) const {
        auto f = [](double x, double y) {
            return x * x * (1 - x) * (1 - x) * (2 * y - 6 * y * y + 4 * y * y * y);
        };
        auto vx = f(p[0], p[1]);
        auto vy = -f(p[1], p[0]);
        return { vx, vy };
    }

    void output_exact() {
        auto p = [this](point_type x) { return exact_p(x); };
        auto vx = [this](point_type x) { return exact_v(x)[0]; };
        auto vy = [this](point_type x) { return exact_v(x)[1]; };

        auto project = [&](auto fun) {
            vector_type rhs{{ Ux.dofs(), Uy.dofs() }};
            compute_projection(rhs, Ux.basis, Uy.basis, [&](double x, double y) { return fun({x, y}); });
            ads_solve(rhs, buffer, Ux.data(), Uy.data());
            return rhs;
        };

        output.to_file(project(p), "pressure_ref.data");
        output.to_file(project(vx), "vx_ref.data");
        output.to_file(project(vy), "vy_ref.data");
    }


    point_type forcing(point_type p) const {
        double x = p[0], y = p[1];

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

        return { fx, fy };
    }

    auto shifted(int n, int k, mumps::problem& problem) const {
        return [&problem,n,k](int i, int j, double val) {
            problem.add(n + i, k + j, val);
        };
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

    void assemble_matrix(mumps::problem& problem) const {
        auto DV = Vx.dofs() * Vy.dofs();
        auto DU = Ux.dofs() * Uy.dofs();

        auto test_vx = shifted(0, 0, problem);
        auto test_vy = shifted(DV, DV, problem);
        auto test_p = shifted(2*DV, 2*DV, problem);

        auto test_vxy = shifted(0, DV, problem);
        auto test_vyx = shifted(DV, 0, problem);

        auto trial_vx = shifted(3*DV, 3*DV, problem);
        auto trial_vy = shifted(3*DV + DU, 3*DV + DU, problem);
        auto trial_p = shifted(3*DV + 2*DU, 3*DV + 2*DU, problem);

        auto hh = h * h;

        // Gram matrix
        for (auto i : dofs(Vx, Vy)) {
            for (auto j : overlapping_dofs(i, Vx, Vy)) {
                int ii = linear_index(i, Vx, Vy) + 1;
                int jj = linear_index(j, Vx, Vy) + 1;

                double M = kron(MVx, MVy, i, j);
                double Kx = kron(KVx, MVy, i, j);
                double Ky = kron(MVx, KVy, i, j);
                double Ax = AVx(j[0], i[0]) * AVy(i[1], j[1]);
                double Ay = AVx(i[0], j[0]) * AVy(j[1], j[0]);

                // test_p(ii, jj, M + hh * (Kx + Ky));
                if (! is_boundary(i, Vx, Vy) && ! is_boundary(j, Vx, Vy)) {
                    test_p(ii, jj, M + hh * (Kx + Ky)); // !!!!!!!!!!!

                    // test_vx(ii, jj, M + hh * (Kx + Ky));
                    // test_vy(ii, jj, M + hh * (Kx + Ky));
                    test_vx(ii, jj, M + hh * Kx);
                    test_vy(ii, jj, M + hh * Ky);
                    test_vxy(ii, jj, hh * Ax);
                    test_vyx(ii, jj, hh * Ay);
                }
            }
        }

        // Dirichlet BC - upper left
        for_boundary_dofs(Vx, Vy, [&](index_type dof) {
            int i = linear_index(dof, Vx, Vy) + 1;
            test_p(i, i, 1); // !!!!!!!!

            test_vx(i, i, 1);
            test_vy(i, i, 1);
        });

        // B, B^t
        auto put = [&](index_type i, index_type j, int sx, int sy, double val, bool bc = false) {
            bc = true; // !!!!!!!!!!!!
            int ii = linear_index(i, Vx, Vy) + 1 + sx;
            int jj = linear_index(j, Ux, Uy) + 1 + sy;

            if (! bc || ! is_boundary(i, Vx, Vy)) {
                problem.add(ii, 3*DV + jj, -val);
            }
            if (! bc || (! is_boundary(i, Vx, Vy) && ! is_boundary(j, Ux, Uy))) {
                problem.add(3*DV + jj, ii, val);
            }
        };

        for (auto i : dofs(Vx, Vy)) {
            for (auto j : dofs(Ux, Uy)) {
                if (! overlap(i, Vx, Vy, j, Ux, Uy)) continue;
                double Kx = kron(KUVx, MUVy, i, j);
                double Ky = kron(MUVx, KUVy, i, j);
                double Ax = kron(AUVx, MUVy, i, j);
                double Ay = kron(MUVx, AUVy, i, j);

                // B(w, u)
                // w = (tx, ty, w)
                // u = (vx, vy, p)

                // tx, vx -> (\/tx, \/vx)
                put(i, j, 0, 0, Kx + Ky, true);

                // ty, vy -> (\/ty, \/vy)
                put(i, j, DV, DU, Kx + Ky, true);

                // w, vx -> (w, vx,x)
                put(i, j, 2*DV, 0, Ax);

                // w, vy -> (w, vy,y)
                put(i, j, 2*DV, DU, Ay);

                // tx, p -> (tx, p,x)
                put(i, j, 0, 2*DU, Ax, true);

                // ty, p -> (ty, p,y)
                put(i, j, DV, 2*DU, Ay, true);
            }
        }

        // Dirichlet BC - lower right
        for_boundary_dofs(Ux, Uy, [&](index_type dof) {
            int i = linear_index(dof, Ux, Uy) + 1;
            trial_p(i, i, 1); // !!!!!!!!!!!

            trial_vx(i, i, 1);
            trial_vy(i, i, 1);
        });
    }

    void compute_rhs(vector_view& vx, vector_view& vy, vector_view& p) const {
        using shape = std::array<std::size_t, 2>;
        auto local_shape = shape{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() };

        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto vx_loc = vector_type{ local_shape };
            auto vy_loc = vector_type{ local_shape };
            // auto Lp = vector_type{ local_shape };

            double J = jacobian(e);
            for (auto q : quad_points(Vx, Vy)) {
                double W = weigth(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e, Vx, Vy)) {
                    auto aa = dof_global_to_local(e, a, Vx, Vy);
                    value_type v = eval_basis(e, q, a, Vx, Vy);

                    auto F = forcing(x);
                    double Lvx = F[0] * v.val;
                    double Lvy = F[1] * v.val;

                    vx_loc(aa[0], aa[1]) -= Lvx * W * J;
                    vy_loc(aa[0], aa[1]) -= Lvy * W * J;
                    // Lp(aa[0], aa[1]) += val * W * J;
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(vx, vx_loc, e, Vx, Vy);
                update_global_rhs(vy, vy_loc, e, Vx, Vy);
                // update_global_rhs(p, Lp, e, Vx, Vy);
            });
        });
    }


    void apply_bc(vector_view& Rvx, vector_view& Rvy, vector_view& Rp) const {
        zero_bc(Rvx, Vx, Vy);
        zero_bc(Rvy, Vx, Vy);
        zero_bc(Rp, Vx, Vy);
    }

    void step(int /*iter*/, double /*t*/) override {
        auto trial_dim = Ux.dofs() * Uy.dofs();
        auto test_dim = Vx.dofs() * Vy.dofs();

        using shape = std::array<std::size_t, 2>;
        auto test_shape = shape{ Vx.dofs(), Vy.dofs() };
        auto trial_shape = shape{ Ux.dofs(), Uy.dofs() };


        std::vector<double> rhs(3 * (test_dim + trial_dim));

        vector_view Rvx{rhs.data(), test_shape};
        vector_view Rvy{Rvx.data() + test_dim, test_shape};
        vector_view Rp{Rvy.data() + test_dim, test_shape};

        vector_view vx{Rp.data() + test_dim, trial_shape};
        vector_view vy{vx.data() + trial_dim, trial_shape};
        vector_view p{vy.data() + trial_dim, trial_shape};


        mumps::problem problem(rhs.data(), rhs.size());

        std::cout << "Assembling matrix" << std::endl;
        assemble_matrix(problem);

        std::cout << "Computing RHS" << std::endl;
        compute_rhs(Rvx, Rvy, Rp);
        apply_bc(Rvx, Rvy, Rp);

        std::cout << "Solving" << std::endl;
        solver.solve(problem);

        std::cout << "Outputting" << std::endl;
        output.to_file(p, "pressure.data");
        output.to_file(vx, "vx.data");
        output.to_file(vy, "vy.data");
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) const {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

};

}

#endif // PROBLEMS_STOKES_STOKES_HPP_
