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

    dimension ref_x, ref_y;

    lin::band_matrix MVx, MVy;
    lin::band_matrix KVx, KVy;
    lin::band_matrix AVx, AVy;

    lin::dense_matrix MUVx, MUVy;
    lin::dense_matrix KUVx, KUVy;
    lin::dense_matrix AUVx, AUVy;

    lin::dense_matrix AVUx, AVUy;


    double h;
    vector_type buffer;

    mumps::solver solver;
    Galois::StatTimer solver_timer{"solver"};
    output_manager<2> output;

public:
    stokes(dimension trial_x, dimension trial_y, dimension test_x, dimension test_y, const timesteps_config& steps,
           dimension ref_x, dimension ref_y)
    : Base{ std::move(test_x), std::move(test_y), steps }
    , Ux{ std::move(trial_x) }
    , Uy{ std::move(trial_y) }
    , ref_x{ ref_x }
    , ref_y{ ref_y }
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
    , AVUx{ Ux.dofs(), Vx.dofs() }
    , AVUy{ Uy.dofs(), Vy.dofs() }
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
        advection_matrix_1d(AVUx, Vx.basis, Ux.basis);
        advection_matrix_1d(AVUy, Vy.basis, Uy.basis);
    }

    void before() override {
        prepare_matrices();

        Ux.factorize_matrix();
        Uy.factorize_matrix();

        output_exact();
    }

    value_type exact_p(point_type p) const {
        auto x = p[0];
        return {x * (1 - x), 1 - 2 * x, 0.0};
    }

    std::array<value_type, 2> exact_v(point_type p) const {
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

        return { vx ,vy };
    }

    value_type exact_div(point_type p) const {
        auto v = exact_v(p);
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

        return { div, dx, dy };
    }

    void output_exact() {
        auto p = [this](point_type x) { return exact_p(x).val; };
        auto vx = [this](point_type x) { return exact_v(x)[0].val; };
        auto vy = [this](point_type x) { return exact_v(x)[1].val; };

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

    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p) const {
        auto e_vx = [this](point_type x) { return exact_v(x)[0]; };
        auto e_vy = [this](point_type x) { return exact_v(x)[1]; };
        auto e_p = [this](point_type x) { return exact_p(x); };
        auto div = [this](point_type x) { return exact_div(x); };

        double vxL2 = errorL2(vx, Ux, Uy, e_vx) / normL2(Ux, Uy, e_vx) * 100;
        double vxH1 = errorH1(vx, Ux, Uy, e_vx) / normH1(Ux, Uy, e_vx) * 100;

        double vyL2 = errorL2(vy, Ux, Uy, e_vy) / normL2(Ux, Uy, e_vy) * 100;
        double vyH1 = errorH1(vy, Ux, Uy, e_vy) / normH1(Ux, Uy, e_vy) * 100;

        double vL2 = std::hypot(vxL2, vyL2) / std::sqrt(2);
        double vH1 = std::hypot(vxH1, vyH1) / std::sqrt(2);


        double pL2 = errorL2(p, Ux, Uy, e_p) / normL2(Ux, Uy, e_p) * 100;
        double pH1 = errorH1(p, Ux, Uy, e_p) / normH1(Ux, Uy, e_p) * 100;

        double divL2 = div_errorL2(vx, vy, Ux, Uy, div) * 100;
        double divH1 = div_errorH1(vx, vy, Ux, Uy, div) * 100;

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
        std::cout << "v   : L2 = " << vL2    << "%, H1 = " << vH1    << "%" << std::endl;
        std::cout << "p   : L2 = " << pL2    << "%, H1 = " << pH1    << "%" << std::endl;
        std::cout << "div : L2 = " << divL2  << ", H1 = " << divH1  << std::endl;
    }




    void print_error(const vector_view& vx, const vector_view& vy, const vector_view& p,
                     const vector_type& ref_vx, const vector_type& ref_vy, const vector_type& ref_p) const {
        bspline::eval_ders_ctx cx{ref_x.p, 1};
        bspline::eval_ders_ctx cy{ref_y.p, 1};

        auto e_vx = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_vx, ref_x.B, ref_y.B, cx, cy); };
        auto e_vy = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_vy, ref_x.B, ref_y.B, cx, cy); };
        auto e_p = [&,this](point_type x) { return eval_ders(x[0], x[1], ref_p, ref_x.B, ref_y.B, cx, cy); };
        auto div = [this](point_type x) { return value_type{}; };

        double vxL2 = errorL2(vx, Ux, Uy, e_vx) / normL2(Ux, Uy, e_vx) * 100;
        double vxH1 = errorH1(vx, Ux, Uy, e_vx) / normH1(Ux, Uy, e_vx) * 100;

        double vyL2 = errorL2(vy, Ux, Uy, e_vy) / normL2(Ux, Uy, e_vy) * 100;
        double vyH1 = errorH1(vy, Ux, Uy, e_vy) / normH1(Ux, Uy, e_vy) * 100;

        double vL2 = std::hypot(vxL2, vyL2) / std::sqrt(2);
        double vH1 = std::hypot(vxH1, vyH1) / std::sqrt(2);


        double pL2 = errorL2(p, Ux, Uy, e_p) / normL2(Ux, Uy, e_p) * 100;
        double pH1 = errorH1(p, Ux, Uy, e_p) / normH1(Ux, Uy, e_p) * 100;

        double divL2 = div_errorL2(vx, vy, Ux, Uy, div) * 100;
        double divH1 = div_errorH1(vx, vy, Ux, Uy, div) * 100;

        std::cout.precision(3);
        std::cout << "vx  : L2 = " << vxL2   << "%, H1 = " << vxH1   << "%" << std::endl;
        std::cout << "vy  : L2 = " << vyL2   << "%, H1 = " << vyH1   << "%" << std::endl;
        std::cout << "v   : L2 = " << vL2    << "%, H1 = " << vH1    << "%" << std::endl;
        std::cout << "p   : L2 = " << pL2    << "%, H1 = " << pH1    << "%" << std::endl;
        std::cout << "div : L2 = " << divL2  << ", H1 = " << divH1  << std::endl;
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
        // return { 0, 0 };
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

    bool is_pressure_fixed(index_type dof) const {
        return dof[0] == 0 && dof[1] == 0;
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
                double Ay = AVx(i[0], j[0]) * AVy(j[1], i[1]);

                // G(w, u)
                // w = (tx, ty, w)
                // u = (vx, vy, p)

                if (! is_pressure_fixed(i) && ! is_pressure_fixed(j)) {
                    // w, p -> (w, p) // + hh (\/w, \/p)
                    test_p(ii, jj, M/* + hh * (Kx + Ky)*/);
                }

                if (! is_boundary(i, Vx, Vy) && ! is_boundary(j, Vx, Vy)) {
                    // tx, vx -> (tx, vx) + hh (tx,x, vx,x)
                    // tx, vx -> (\/tx, \/vx)
                    test_vx(ii, jj, /*M + hh */ Kx + Ky);

                    // ty, vy -> (ty, vy) + hh (ty,y, vy,y)
                    // ty, vy -> (\/ty, \/vy)
                    test_vy(ii, jj, /*M + hh */ Kx + Ky);

                    // tx, vy -> hh (tx,x, vy,y)
                    // test_vxy(ii, jj, hh * Ax);

                    // ty, vx -> hh (ty,y, vx,x)
                    // test_vyx(ii, jj, hh * Ay);
                }
            }
        }

        // Dirichlet BC - test space
        for_boundary_dofs(Vx, Vy, [&](index_type dof) {
            int i = linear_index(dof, Vx, Vy) + 1;
            test_vx(i, i, 1);
            test_vy(i, i, 1);
        });
        test_p(1, 1, 1.0);

        // B, B^t
        auto put = [&](index_type i, index_type j, int si, int sj, double val, bool fixed_i, bool fixed_j) {
            int ii = linear_index(i, Vx, Vy) + 1 + si;
            int jj = linear_index(j, Ux, Uy) + 1 + sj;

            if (!fixed_i) {
                problem.add(ii, 3*DV + jj, -val);
            }
            if (!fixed_i && !fixed_j) {
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
                double AxT = AVUx(j[0], i[0]) * MUVy(i[1], j[1]);
                double AyT = MUVx(i[0], j[0]) * AVUy(j[1], i[1]);

                bool bd_i = is_boundary(i, Vx, Vy);
                bool bd_j = is_boundary(j, Ux, Uy);
                bool fixed_i = is_pressure_fixed(i);
                bool fixed_j = is_pressure_fixed(j);

                // B(w, u)
                // w = (tx, ty, w)
                // u = (vx, vy, p)

                // tx, vx -> (\/tx, \/vx)
                put(i, j, 0, 0, Kx + Ky, bd_i, bd_j);

                // ty, vy -> (\/ty, \/vy)
                put(i, j, DV, DU, Kx + Ky, bd_i, bd_j);

                // w, vx -> (w, vx,x)
                put(i, j, 2*DV, 0, Ax, fixed_i, bd_j);

                // w, vy -> (w, vy,y)
                put(i, j, 2*DV, DU, Ay, fixed_i, bd_j);

                // tx, p -> (tx, p,x) = - (tx,x, p)
                // put(i, j, 0, 2*DU, Ax, bd_i, fixed_j);
                put(i, j, 0, 2*DU, -AxT, bd_i, fixed_j);

                // ty, p -> (ty, p,y) = - (ty,y, p)
                // put(i, j, DV, 2*DU, Ay, bd_i, fixed_j);
                put(i, j, DV, 2*DU, -AyT, bd_i, fixed_j);
            }
        }

        // Dirichlet BC - trial space
        for_boundary_dofs(Ux, Uy, [&](index_type dof) {
            int i = linear_index(dof, Ux, Uy) + 1;
            trial_vx(i, i, 1);
            trial_vy(i, i, 1);
        });
        trial_p(1, 1, 1.0);
    }

    void compute_rhs(vector_view& vx, vector_view& vy, vector_view& /*p*/) const {
        using shape = std::array<std::size_t, 2>;
        auto local_shape = shape{ Vx.basis.dofs_per_element(), Vy.basis.dofs_per_element() };

        executor.for_each(elements(Vx, Vy), [&](index_type e) {
            auto vx_loc = vector_type{ local_shape };
            auto vy_loc = vector_type{ local_shape };

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
                }
            }
            executor.synchronized([&]() {
                update_global_rhs(vx, vx_loc, e, Vx, Vy);
                update_global_rhs(vy, vy_loc, e, Vx, Vy);
            });
        });
    }


    void apply_bc(vector_view& Rvx, vector_view& Rvy, vector_view& Rp) {
        zero_bc(Rvx, Ux, Uy);

        // constexpr double eps = 1e-4;
        // auto drop = [](double t) { return t < 1 - eps ? 0 : 1 - (1 - t) / eps; };
        // dirichlet_bc(Rvx, boundary::left, Ux, Uy, drop);
        // dirichlet_bc(Rvx, boundary::right, Ux, Uy, drop);

        // // dirichlet_bc(Rvx, boundary::left,   Ux, Uy, 0);
        // // dirichlet_bc(Rvx, boundary::right,  Ux, Uy, 0);
        // dirichlet_bc(Rvx, boundary::top,    Ux, Uy, 1);
        // dirichlet_bc(Rvx, boundary::bottom, Ux, Uy, 0);

        zero_bc(Rvy, Ux, Uy);
        Rp(0, 0) = 0; // fix pressure at a point
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
        apply_bc(vx, vy, p);


        std::cout << "Solving" << std::endl;
        solver_timer.start();
        solver.solve(problem);
        solver_timer.stop();

        std::cout << "  solver time:       " << static_cast<double>(solver_timer.get()) << " ms" << std::endl;
        std::cout << "  assembly    FLOPS: " << solver.flops_assembly() << std::endl;
        std::cout << "  elimination FLOPS: " << solver.flops_elimination() << std::endl;

        std::cout << "Error:" << std::endl;
        print_error(vx, vy, p);

        std::cout << "Outputting" << std::endl;
        output.to_file(p, "pressure.data");
        output.to_file(vx, "vx.data");
        output.to_file(vy, "vy.data");

        // print_solution("vx.sol", vx, Ux, Uy);
        // print_solution("vy.sol", vy, Ux, Uy);
        // print_solution("p.sol", p, Ux, Uy);

        // auto ref_vx = read_solution("vx.sol", ref_x, ref_y);
        // auto ref_vy = read_solution("vy.sol", ref_x, ref_y);
        // auto ref_p = read_solution("p.sol", ref_x, ref_y);

        // print_error(vx, vy, p, ref_vx, ref_vy, ref_p);
    }

    vector_type read_solution(const char* filename, const dimension& Ux, const dimension& Uy) const {
        auto sol = vector_type{{ Ux.dofs(), Uy.dofs() }};
        std::ifstream in{filename};
        while (in) {
            index_type dof;
            in >> dof[0] >> dof[1];
            in >> sol(dof[0], dof[1]);
        }
        return sol;
    }

    template <typename V>
    void print_solution(const char* filename, const V& u, const dimension& Ux, const dimension& Uy) const {
        std::ofstream out{filename};
        for (auto dof : dofs(Ux, Uy)) {
            out << dof[0] << " " << dof[1] << " " << std::setprecision(16) << u(dof[0], dof[1]) << std::endl;
        }
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) const {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

    template <typename Sol, typename Fun>
    double div_errorL2(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        auto L2 = [](value_type a) { return a.val * a.val; };
        return div_error(u, v, Ux, Uy, L2, fun);
    }

    template <typename Sol, typename Fun>
    double div_errorH1(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Fun&& fun) const {
        auto H1 = [](value_type a) { return a.val * a.val + a.dx * a.dx + a.dy * a.dy; };
        return div_error(u, v, Ux, Uy, H1, fun);
    }

    template <typename Sol, typename Fun, typename Norm>
    double div_error(const Sol& u, const Sol& v, const dimension& Ux, const dimension& Uy, Norm&& norm, Fun&& fun) const {
        double error = 0;

        for (auto e : elements(Ux, Ux)) {
            double J = jacobian(e, Ux, Uy);
            for (auto q : quad_points(Ux, Uy)) {
                double w = weigth(q, Ux, Uy);
                auto x = point(e, q, Ux, Uy);
                value_type div = divergence(u, v, e, q, Ux, Uy);

                auto d = div - fun(x);
                error += norm(d) * w * J;
            }
        }
        return std::sqrt(error);
    }

    template <typename Sol>
    value_type divergence(const Sol& u, const Sol& v, index_type e, index_type q,
                          const dimension& x, const dimension& y) const {
        value_type div{};
        for (auto b : dofs_on_element(e, x, y)) {
            double c = u(b[0], b[1]);
            double d = v(b[0], b[1]);

            auto loc = dof_global_to_local(e, b, x, y);

            const auto& bx = x.basis;
            const auto& by = y.basis;

            double B1  = bx.b[e[0]][q[0]][0][loc[0]];
            double B2  = by.b[e[1]][q[1]][0][loc[1]];
            double dB1 = bx.b[e[0]][q[0]][1][loc[0]];
            double dB2 = by.b[e[1]][q[1]][1][loc[1]];
            double ddB1 = bx.b[e[0]][q[0]][2][loc[0]];
            double ddB2 = by.b[e[1]][q[1]][2][loc[1]];

            double dx = dB1 *  B2;
            double dy =  B1 * dB2;
            double dxx = ddB1 * B2;
            double dyy = B1 * ddB2;
            double dxy = dB1 * dB2;

            div.val += c * dx + d * dy;
            div.dx += c * dxx + d * dxy;
            div.dy += c * dxy + d * dyy;
        }
        return div;
    }

};

}

#endif // PROBLEMS_STOKES_STOKES_HPP_
