// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <iostream>

#include <lyra/lyra.hpp>

#include "ads/experimental/all.hpp"
#include "ads/experimental/horrible_sparse_matrix.hpp"
#include "ads/experimental/space_factory.hpp"
#include "ads/simulation.hpp"

template <typename Fun>
auto solution_to_file(std::string const& path, Fun&& fun) -> void {
    constexpr auto res = 200;
    auto out = fmt::output_file(path);

    out.print("x y z u\n");
    for (auto const y : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto const x : ads::evenly_spaced(0.0, 1.0, res)) {
            out.print("{} {} {}", x, y, 0);
            out.print(" {:.7}\n", fun({x, y}));
        }
    }
}

auto report_cost(ads::mumps::solver const& solver) -> void {
    fmt::print("FLOPS: {:5.3e}\n", solver.flops_elimination());
}

auto neg_part(double x) -> double {
    return (std::abs(x) - x) / 2;
}

struct element_idx_hasher {
    using type = ads::regular_mesh::element_index;
    auto operator()(type const& e) const noexcept -> std::size_t {
        auto const [i, j] = e;
        return i ^ (j + 0x9e3779b9 + (i << 6) + (i >> 2));
    }
};

struct facet_idx_hasher {
    using type = ads::regular_mesh::facet_index;
    auto operator()(type const& e) const noexcept -> std::size_t {
        auto const [i, j, o] = e;
        auto const oi = static_cast<int>(o);
        auto val = i ^ (j + 0x9e3779b9 + (i << 6) + (i >> 2));
        return val ^ (oi + 0x9e3779b9 + (val << 6) + (val >> 2));
    }
};

struct facet_idx_equality {
    using type = ads::regular_mesh::facet_index;
    auto operator()(type const& a, type const& b) const noexcept -> bool {
        return a.ix == b.ix && a.iy == b.iy && a.dir == b.dir;
    }
};

auto mark_for_refinement(std::vector<double> const& local, double eta) -> std::vector<int> {
    auto const n = local.size();
    auto elems = std::vector<int>(n);

    auto indices = std::vector<std::vector<double>::size_type>(n);
    std::iota(begin(indices), end(indices), 0);
    std::sort(begin(indices), end(indices), [&](auto i, auto j) { return local[i] > local[j]; });

    auto const total = std::accumulate(begin(local), end(local), 0.0);
    auto const enough = total * eta;
    auto current = 0.0;

    fmt::print("Refined: ");
    for (auto idx : indices) {
        elems[idx] = 1;
        auto const val = local[idx];
        fmt::print("{} ", idx);
        current += val;
        if (current > enough) {
            break;
        }
    }
    fmt::print("\n");

    return elems;
}

auto refine_single_dimension(ads::partition const& mesh, std::vector<int> mark) -> ads::partition {
    auto new_mesh = ads::partition{};
    for (auto i = 0; i < ads::as_signed(mesh.size()) - 1; ++i) {
        new_mesh.push_back(mesh[i]);
        if (mark[i]) {
            auto const mid = 0.5 * (mesh[i] + mesh[i + 1]);
            new_mesh.push_back(mid);
        }
    }
    new_mesh.push_back(mesh.back());
    return new_mesh;
}

auto refine_single_dimension(ads::partition const& mesh) -> ads::partition {
    auto const n = ads::as_signed(mesh.size()) - 1;
    auto elems = std::vector<int>(n);
    std::fill(begin(elems), end(elems), 1);
    return refine_single_dimension(mesh, elems);
}

auto minimum_diameter(ads::partition const& xs) -> double {
    auto h = std::numeric_limits<double>::max();
    for (int i = 0; i < ads::as_signed(xs.size()) - 1; ++i) {
        auto dx = xs[i + 1] - xs[i];
        h = std::min(h, dx);
    }
    return h;
}

auto ensure_aspect_ratio(ads::partition const& xs, ads::partition const& ys, double ratio)
    -> ads::partition {
    auto const diam_x = minimum_diameter(xs);
    auto const diam_y = minimum_diameter(ys);
    if (diam_y / diam_x <= ratio) {
        return ys;
    } else {
        return refine_single_dimension(ys);
    }
}

template <typename Facets, typename Space, typename Quad, typename Out, typename Form>
auto assemble_facets(const Facets& facets, const Space& space, const Quad& quad, int ders, Out out,
                     Form&& form) -> void {
    auto const& mesh = space.mesh();

    for (auto const f : facets) {
        auto const facet = mesh.facet(f);
        auto const points = quad.coordinates(f);
        auto const eval = space.dof_evaluator(f, points, std::max(ders, 1));

        auto const n = space.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        using dof_idx = typename Space::dof_index;
        using point_idx = typename decltype(points)::point_index;
        using value_type =
            decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>(), facet.normal));

        auto basis_vals = std::vector<value_type>(n);
        auto basis_nders = std::vector<ads::facet_value<double>>(n);

        for (auto const q : points.indices()) {
            auto const [x, w] = points.data(q);
            for (auto const i : space.dofs_on_facet(f)) {
                auto const iloc = space.facet_local_index(i, f);
                basis_vals[iloc] = eval(i, q, facet.normal);
                basis_nders[iloc] = eval.normal_derivative(i, q, ders, facet.normal);
            }
            for (int iloc = 0; iloc < n; ++iloc) {
                for (int jloc = 0; jloc < n; ++jloc) {
                    const auto& u = basis_vals[iloc];
                    const auto& v = basis_vals[jloc];
                    const auto& u_nder = basis_nders[iloc];
                    const auto& v_nder = basis_nders[jloc];
                    M(jloc, iloc) += form(u, v, u_nder, v_nder, x, facet) * w;
                }
            }
        }

        for (auto i : space.dofs_on_facet(f)) {
            const auto iloc = space.facet_local_index(i, f);
            const auto I = space.global_index(i);
            for (auto j : space.dofs_on_facet(f)) {
                const auto jloc = space.facet_local_index(j, f);
                const auto J = space.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Facets, typename Trial, typename Test, typename Quad, typename Out,
          typename Form>
auto assemble_facets(const Facets& facets, const Trial& trial, const Test& test, const Quad& quad,
                     int ders, Out out, Form&& form) -> void {
    const auto& mesh = test.mesh();

    for (auto f : facets) {
        const auto facet = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval_trial = trial.dof_evaluator(f, points, std::max(ders, 1));
        const auto eval_test = test.dof_evaluator(f, points, std::max(ders, 1));

        const auto n_trial = trial.facet_dof_count(f);
        const auto n_test = test.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n_test, n_trial}};

        using point_idx = typename decltype(points)::point_index;
        using trial_dof_idx = typename Trial::dof_index;
        using test_dof_idx = typename Test::dof_index;
        using trial_value = decltype(eval_trial(std::declval<trial_dof_idx>(),
                                                std::declval<point_idx>(), facet.normal));
        using test_value = decltype(eval_test(std::declval<test_dof_idx>(),
                                              std::declval<point_idx>(), facet.normal));

        auto test_vals = std::vector<test_value>(n_test);
        auto trial_vals = std::vector<trial_value>(n_trial);
        auto test_nders = std::vector<ads::facet_value<double>>(n_test);
        auto trial_nders = std::vector<ads::facet_value<double>>(n_trial);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : trial.dofs_on_facet(f)) {
                const auto iloc = trial.facet_local_index(i, f);
                trial_vals[iloc] = eval_trial(i, q, facet.normal);
                trial_nders[iloc] = eval_trial.normal_derivative(i, q, ders, facet.normal);
            }
            for (auto j : test.dofs_on_facet(f)) {
                const auto jloc = test.facet_local_index(j, f);
                test_vals[jloc] = eval_test(j, q, facet.normal);
                test_nders[jloc] = eval_test.normal_derivative(j, q, ders, facet.normal);
            }
            for (int iloc = 0; iloc < n_trial; ++iloc) {
                for (int jloc = 0; jloc < n_test; ++jloc) {
                    const auto& u = trial_vals[iloc];
                    const auto& v = test_vals[jloc];
                    const auto& u_nder = trial_nders[iloc];
                    const auto& v_nder = test_nders[jloc];
                    M(jloc, iloc) += form(u, v, u_nder, v_nder, x, facet) * w;
                }
            }
        }

        for (auto i : trial.dofs_on_facet(f)) {
            const auto iloc = trial.facet_local_index(i, f);
            const auto I = trial.global_index(i);
            for (auto j : test.dofs_on_facet(f)) {
                const auto jloc = test.facet_local_index(j, f);
                const auto J = test.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Coeffs, typename Space, typename Quad, typename Out, typename Norm>
auto local_contribution(Coeffs const& coeffs, Space const& space, Quad const& quad, Out out,
                        Norm&& norm) -> void {
    const auto& mesh = space.mesh();

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        const auto eval = space.dof_evaluator(e, points, 1);

        using dof_idx = typename Space::dof_index;
        using point_idx = typename decltype(points)::point_index;
        using value_type = decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>()));

        auto norm_val = 0.0;

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);

            auto val = value_type{};

            for (auto const i : space.dofs(e)) {
                const auto I = space.global_index(i);
                const auto u = eval(i, q);
                val += coeffs[I] * u;
            }
            norm_val += norm(val, x) * w;
        }
        out(e, norm_val);
    }
}

template <typename Coeffs, typename Facets, typename Space, typename Quad, typename Out,
          typename Norm>
auto local_contribution_facets(Coeffs const& coeffs, Facets const& facets, Space const& space,
                               Quad const& quad, int ders, Out out, Norm&& norm) -> void {
    auto const& mesh = space.mesh();

    for (auto const f : facets) {
        auto const facet = mesh.facet(f);
        auto const points = quad.coordinates(f);
        auto const eval = space.dof_evaluator(f, points, std::max(ders, 1));

        using dof_idx = typename Space::dof_index;
        using point_idx = typename decltype(points)::point_index;
        using value_type =
            decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>(), facet.normal));

        auto norm_val = 0.0;

        for (auto const q : points.indices()) {
            auto const [x, w] = points.data(q);

            auto val = value_type{};
            auto val_nder = ads::facet_value<double>{};

            for (auto const i : space.dofs_on_facet(f)) {
                const auto I = space.global_index(i);
                auto const u = eval(i, q, facet.normal);
                auto const u_nder = eval.normal_derivative(i, q, ders, facet.normal);
                val.avg += coeffs[I] * u.avg;
                val.jump += coeffs[I] * u.jump;
                val_nder.avg += coeffs[I] * u_nder.avg;
                val_nder.jump += coeffs[I] * u_nder.jump;
            }
            norm_val += norm(val, val_nder, x, facet) * w;
        }
        out(f, norm_val);
    }
}

class eriksson_johnson {
public:
    double eps;

    using point = ads::point_t;

    eriksson_johnson(double eps)
    : eps{eps} { }

    auto u(point p) const noexcept -> double { return u_with_grad(p).val; }

    auto u_with_grad(point p) const noexcept -> ads::value_type {
        using std::exp, std::sqrt, std::sin;
        constexpr double pi = M_PI;
        constexpr int n = 1;

        auto const lambda = n * pi * eps;
        auto const del = sqrt(1 + 4 * lambda * lambda);
        auto const r1 = (1 + del) / (2 * eps);
        auto const r2 = (1 - del) / (2 * eps);

        auto const norm = exp(-r1) - exp(-r2);

        auto const [x, y] = p;
        auto const alpha = (exp(r1 * (x - 1)) - exp(r2 * (x - 1))) / norm;
        auto const val = alpha * sin(n * pi * y);

        return {
            val,
            (r1 * exp(r1 * (x - 1)) - r2 * exp(r2 * (x - 1))) / norm * sin(pi * y),
            alpha * pi * cos(pi * y),
        };
    }

    auto u() const noexcept {
        return [this](auto x) { return u(x); };
    }

    auto u_with_grad() const noexcept {
        return [this](auto x) { return u_with_grad(x); };
    }

    auto f(point /*p*/) const noexcept -> double { return 0; }
};

void validate_args(lyra::cli const& cli, lyra::parse_result const& result, bool show_help) {
    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }
    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }
}

struct galerkin_args {
    double gamma = 1.0;
    double eps = 1e-3;
    int nx = 16;
    int ny = 16;
    int p = 2;
    int c = 0;
    bool weak_bc = false;
    double eta = 1.0;
};

auto parse_args(int argc, char* argv[]) {
    auto args = galerkin_args{};

    bool show_help = false;
    auto const cli = lyra::help(show_help)                                                 //
                   | lyra::opt(args.eps, "epsilon")["--eps"]("diffusion coefficient")      //
                   | lyra::opt(args.gamma, "gamma")["--gamma"]("stabilization parameter")  //
                   | lyra::opt(args.nx, "nx")["--nx"]("mesh size (x direction)")           //
                   | lyra::opt(args.ny, "ny")["--ny"]("mesh size (y direction)")           //
                   | lyra::opt(args.p, "p")["-p"]("order of basis functions")              //
                   | lyra::opt(args.c, "c")["-c"]("continuity of basis functions")         //
                   | lyra::opt(args.weak_bc)["--weak-bc"]("should BC be imposed weakly")   //
                   | lyra::opt(args.eta, "eta")["--eta"]("weak BC parameter")              //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

template <typename Rhs>
auto norm_using_matrix(Rhs const& x, ads::horrible_sparse_matrix const& mat) -> double {
    auto value = 0.0;
    using std::begin, std::end;
    using std::cbegin, std::cend;
    for (auto const [dof, a] : mat) {
        value += a * x[dof.i] * x[dof.j];
    }
    return std::sqrt(value);
}

auto adapted_mesh(double eps) -> ads::regular_mesh {
    auto const n = 32;
    auto const d = std::max(0.5, 1 - 10 * eps);
    auto const ys = ads::evenly_spaced(0.0, 1.0, n);
    auto const left = ads::evenly_spaced(0.0, d, n);
    auto const right = ads::evenly_spaced(d, 1.0, n);

    auto xs = left;
    xs.insert(xs.end(), right.begin() + 1, right.end());
    return ads::regular_mesh{xs, ys};
}

auto adapt_mesh(ads::regular_mesh const& mesh, double eps) -> ads::regular_mesh {
    auto const n = 32;
    auto const d = std::max(0.5, 1 - 10 * eps);

    auto xs = mesh.mesh_x().points();
    xs.pop_back();

    auto const last_point = xs.back();
    if (last_point < d) {
        auto const right = ads::evenly_spaced(d, 1.0, n);
        xs.insert(xs.end(), right.begin() + 1, right.end());
    }
    xs.push_back(1.0);

    auto const ys = ads::evenly_spaced(0.0, 1.0, n);
    return ads::regular_mesh{xs, ys};
}

auto mark_boundary_dofs(ads::space const& U) -> std::vector<int> {
    auto is_fixed = std::vector<int>(U.dof_count());
    for (auto const dof : U.dofs()) {
        auto const [ix, iy] = dof;
        bool fix = false;
        // spatial boundary
        fix = fix || ix == 0 || ix == U.space_x().dof_count() - 1;
        fix = fix || iy == 0 || iy == U.space_y().dof_count() - 1;

        if (fix) {
            is_fixed[U.global_index(dof)] = 1;
        }
    }
    return is_fixed;
}

template <typename Rhs, typename Problem>
auto postprocess(Rhs const& F, ads::space const& U, Problem&& problem) -> void {
    auto const u = ads::bspline_function(&U, F);

    auto const mesh = adapted_mesh(problem.eps);
    auto const quad = ads::quadrature{&mesh, 4};

    auto const normL2 = norm(mesh, quad, L2{}, problem.u());
    auto const errorL2 = error(mesh, quad, L2{}, u, problem.u());

    auto const normH1 = norm(mesh, quad, H1{}, problem.u_with_grad());
    auto const errorH1 = error(mesh, quad, H1{}, u.with_grad(), problem.u_with_grad());

    solution_to_file("result.csv", u);

    auto out = fmt::output_file("coeffs.data");
    for (auto const dof : U.dofs()) {
        auto const [ix, iy] = dof;
        auto const val = F[U.global_index(dof)];
        out.print("{:3} {:3} {:7}\n", ix, iy, val);
    }

    fmt::print("Error: {:7.3f}%  {:7.3f}%\n", 100 * errorL2 / normL2, 100 * errorH1 / normH1);
}

template <typename Rhs, typename Quad, typename Problem, typename Norm, typename BdNorm>
auto report_error(Rhs const& F, ads::regular_mesh const& mesh, Quad const& quad,
                  ads::space const& U, Problem&& problem, Norm&& norm, BdNorm&& bd_norm) -> void {
    auto const u = ads::bspline_function(&U, F);

    // auto const fine_mesh = adapted_mesh(problem.eps);
    auto const fine_mesh = adapt_mesh(mesh, problem.eps);
    auto const fine_quad = ads::quadrature{&fine_mesh, 4};

    const auto bd_fun = [&](auto x, auto const& edge) {
        auto const diff = u.with_grad()(x) - problem.u_with_grad()(x);
        auto const v = bd_norm(diff, x, edge);
        return v;
    };

    auto const err_val = error(fine_mesh, fine_quad, norm, u.with_grad(), problem.u_with_grad());
    auto const err_bd2 = integrate_facets(mesh.boundary_facets(), mesh, quad, bd_fun);
    auto const err = std::sqrt(err_val * err_val + err_bd2);

    fmt::print("Error in V norm:      {}\n", err);
}

auto galerkin(int argc, char* argv[]) -> void {
    auto const args = parse_args(argc, argv);

    auto const elems_x = args.nx;
    auto const elems_y = args.ny;
    auto const p = args.p;
    auto const c = args.c;
    auto const eps = args.eps;
    auto const gamma = args.gamma;
    auto const eta = args.eta;
    auto const beta = ads::point_t{1, 0};
    auto const problem = eriksson_johnson{eps};

    auto const xs = ads::evenly_spaced(0.0, 1.0, elems_x);
    auto const ys = ads::evenly_spaced(0.0, 1.0, elems_y);

    auto const bx = ads::make_bspline_basis(xs, p, c);
    auto const by = ads::make_bspline_basis(ys, p, c);

    auto const mesh = ads::regular_mesh{xs, ys};
    auto const quad = ads::quadrature{&mesh, std::max(p + 1, 2)};

    auto spaces = ads::space_factory{};
    auto const U = spaces.next<ads::space>(&mesh, bx, by);

    auto const N = spaces.dim();
    fmt::print("Dimension: {}\n", N);

    auto F = std::vector<double>(N);
    auto linear_system = ads::mumps::problem(F.data(), N);
    auto solver = ads::mumps::solver{};

    auto mat = ads::horrible_sparse_matrix{};
    auto M = [&mat](int row, int col, double val) { mat(row, col) += val; };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot, ads::grad;
    assemble(U, quad, M, [eps, beta](auto u, auto v, auto /*x*/) {
        return eps * dot(grad(u), grad(v)) + dot(grad(u), beta) * v.val;
    });
    auto const stabilizing_form = [=](auto /*u*/, auto /*v*/, auto u_nder, auto v_nder, auto /*x*/,
                                      const auto& edge) {
        auto const& n = edge.normal;
        auto const h = length(edge.span);
        auto const hfac = std::pow(h, 2 * c + 2);
        auto const bn = std::abs(dot(beta, n));
        return gamma * hfac * bn * jump(u_nder) * jump(v_nder);
    };
    assemble_facets(mesh.interior_facets(), U, quad, c + 1, M, stabilizing_form);

    assemble_rhs(U, quad, rhs, [problem](auto v, auto x) { return problem.f(x) * v.val; });

    if (!args.weak_bc) {
        // Strong BC
        auto is_fixed = mark_boundary_dofs(U);

        fmt::print("Applying BC\n");
        for (auto const e : mesh.facets()) {
            for (auto const i : U.dofs_on_facet(e)) {
                auto const I = U.global_index(i);
                if (is_fixed[I] == 1) {
                    for (auto const j : U.dofs_on_facet(e)) {
                        auto const J = U.global_index(j);
                        mat(I, J) = 0;
                    }
                    mat(I, I) = 1;
                    F[I] = 0;
                }
            }
        }

        // Non-zero Dirichlet BC
        auto dimy = ads::dimension{by, p + 1, 1};
        dimy.factorize_matrix();
        auto buf = ads::lin::vector{{dimy.dofs()}};
        ads::compute_projection(buf, dimy.basis, [problem](auto t) { return problem.u({0, t}); });
        ads::lin::solve_with_factorized(dimy.M, buf, dimy.ctx);
        for (int i = 0; i < dimy.dofs(); ++i) {
            F[U.global_index({0, i})] = buf(i);
        }
    } else {
        // Weak BC
        auto const d = 2;
        auto const g = (p + 1) * (p + d) / static_cast<double>(d);

        auto const a_bd = [=](auto u, auto v, auto const& edge) {
            auto const& n = edge.normal;
            auto const h = length(edge.span);
            auto const etaF = g * 4 / h;
            auto const scale = eps * etaF + neg_part(dot(beta, n));
            return eta * scale * u.val * v.val;
        };

        auto const lhs_form = [=](auto u, auto v, auto /*x*/, auto const& edge) {
            return a_bd(avg(u), avg(v), edge);
        };
        assemble_facets(mesh.boundary_facets(), U, quad, M, lhs_form);

        auto const lhs_sym = [=](auto u, auto v, auto /*x*/, auto const& edge) {
            auto const& n = edge.normal;
            auto const uu = avg(u);
            auto const vv = avg(v);
            return -eps * (dot(grad(uu), n) * vv.val + dot(grad(vv), n) * uu.val);
        };
        assemble_facets(mesh.boundary_facets(), U, quad, M, lhs_sym);

        auto const rhs_form = [=](auto v, auto x, auto const& edge) {
            return a_bd(problem.u_with_grad(x), v, edge);
        };
        assemble_rhs(mesh.boundary_facets(), U, quad, rhs, rhs_form);

        auto const rhs_sym = [=](auto v, auto x, auto const& edge) {
            auto const& n = edge.normal;
            return -eps * dot(grad(v), n) * problem.u(x);
        };
        assemble_rhs(mesh.boundary_facets(), U, quad, rhs, rhs_sym);
    }

    // // H1 projection
    // auto const h1 = [](auto u, auto v, auto /*x*/) {
    //     return u.val * v.val + dot(grad(u), grad(v));
    // };

    // assemble(U, quad, M, h1);
    // assemble_rhs(U, quad, rhs, [&](auto v, auto x) {
    //     auto const u = problem.u_with_grad(x);
    //     return h1(u, v, x);
    // });

    mat.mumpsify(linear_system);
    solver.solve(linear_system);
    report_cost(solver);

    postprocess(F.data(), U, problem);
}

struct igrm_args : galerkin_args {
    int P;
    int C;
};

auto parse_args_igrm(int argc, char* argv[]) {
    auto args = igrm_args{};

    bool show_help = false;
    auto const cli = lyra::help(show_help)                                                 //
                   | lyra::opt(args.eps, "epsilon")["--eps"]("diffusion coefficient")      //
                   | lyra::opt(args.gamma, "gamma")["--gamma"]("stabilization parameter")  //
                   | lyra::opt(args.nx, "nx")["--nx"]("mesh size (x direction)")           //
                   | lyra::opt(args.ny, "ny")["--ny"]("mesh size (y direction)")           //
                   | lyra::opt(args.p, "p")["-p"]("order of trial functions")              //
                   | lyra::opt(args.c, "c")["-c"]("continuity of trial functions")         //
                   | lyra::opt(args.P, "P")["-P"]("order of test functions")               //
                   | lyra::opt(args.C, "C")["-C"]("continuity of test functions")          //
                   | lyra::opt(args.weak_bc)["--weak-bc"]("should BC be imposed weakly")   //
                   | lyra::opt(args.eta, "eta")["--eta"]("weak BC parameter")              //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

auto igrm(int argc, char* argv[]) -> void {
    auto const args = parse_args_igrm(argc, argv);

    auto const elems_x = args.nx;
    auto const elems_y = args.ny;
    auto const p = args.p;
    auto const c = args.c;
    auto const eps = args.eps;
    auto const gamma = args.gamma;
    auto const eta = args.eta;
    auto const beta = ads::point_t{1, 0};
    auto const problem = eriksson_johnson{eps};

    auto const p_test = args.P;
    auto const c_test = args.C;

    auto const xs = ads::evenly_spaced(0.0, 1.0, elems_x);
    auto const ys = ads::evenly_spaced(0.0, 1.0, elems_y);

    auto const bx = ads::make_bspline_basis(xs, p, c);
    auto const by = ads::make_bspline_basis(ys, p, c);

    auto const Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto const By = ads::make_bspline_basis(ys, p_test, c_test);

    auto const mesh = ads::regular_mesh{xs, ys};
    auto const quad = ads::quadrature{&mesh, std::max(p_test + 1, 2)};

    auto const U = ads::space{&mesh, bx, by};
    auto const V = ads::space{&mesh, Bx, By};

    auto const n = U.dof_count();
    auto const N = V.dof_count();
    fmt::print("Dimension: {} + {} = {}\n", n, N, n + N);

    auto F = std::vector<double>(N + n);
    auto linear_system = ads::mumps::problem{F.data(), n + N};
    auto solver = ads::mumps::solver{};

    auto mat = ads::horrible_sparse_matrix{};
    auto G = [&mat](int row, int col, double val) { mat(row, col) += val; };
    auto B = [N, &mat](int row, int col, double val) {
        mat(row, N + col) += val;
        mat(N + col, row) += val;
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot, ads::grad;

    auto const der_jump_form = [=](auto /*u*/, auto /*v*/, auto u_nder, auto v_nder, auto /*x*/,
                                   const auto& edge) {
        auto const& n = edge.normal;
        auto const h = length(edge.span);
        auto const hfac = std::pow(h, 2 * c + 2);
        auto const bn = std::abs(dot(beta, n));
        return gamma * hfac * bn * jump(u_nder) * jump(v_nder);
    };

    auto const dg_jump_form = [=](auto u, auto v, auto /*x*/, auto const& edge) {
        auto const& n = edge.normal;
        auto const h = length(edge.span);
        auto const bn = std::abs(dot(beta, n));
        return (gamma / h + bn / 2) * jump(u).val * jump(v).val;
    };

    auto t_before_integration = std::chrono::steady_clock::now();

    // Test scalar product
    assemble(V, quad, G, [=](auto u, auto v, auto /*x*/) {
        auto const beta_inf = 1;
        auto const L = 1;
        auto const h_T = 1.0 / elems_x;
        return eps * dot(grad(u), grad(v))   //
             + beta_inf / L * u.val * v.val  //
             + h_T / beta_inf * dot(beta, grad(u)) * dot(beta, grad(v));
    });
    if (c_test == -1) {
        assemble_facets(mesh.interior_facets(), V, quad, G, dg_jump_form);
    } else {
        assemble_facets(mesh.interior_facets(), V, quad, c_test + 1, G, der_jump_form);
    }

    // Bilinear form of the problem
    assemble(U, V, quad, B, [eps, beta](auto u, auto v, auto /*x*/) {
        return eps * dot(grad(u), grad(v)) + dot(grad(u), beta) * v.val;
    });
    assemble_facets(mesh.interior_facets(), U, V, quad, c_test + 1, B, der_jump_form);

    // RHS
    assemble_rhs(V, quad, rhs, [problem](auto v, auto x) { return problem.f(x) * v.val; });

    if (args.weak_bc) {
        auto const d = 2;
        auto const g = (p + 1) * (p + d) / static_cast<double>(d);

        auto const test_prod_bd_form = [=](auto u, auto v, auto /*x*/, const auto& edge) {
            auto const& n = edge.normal;
            auto const h = length(edge.span);
            auto const etaF = g * 4 / h;
            auto const scale = eps * etaF + 0.5 * abs(dot(beta, n));
            return eta * scale * avg(u).val * avg(v).val;
        };
        assemble_facets(mesh.boundary_facets(), V, quad, G, test_prod_bd_form);

        auto const a_bd = [=](auto u, auto v, auto const& edge) {
            auto const& n = edge.normal;
            auto const h = length(edge.span);
            auto const etaF = g * 4 / h;
            auto const scale = eps * etaF + neg_part(dot(beta, n));
            return eta * scale * u.val * v.val;
        };

        auto const lhs_form = [=](auto u, auto v, auto /*x*/, auto const& edge) {
            return a_bd(avg(u), avg(v), edge);
        };
        assemble_facets(mesh.boundary_facets(), U, V, quad, B, lhs_form);

        auto const lhs_sym = [=](auto u, auto v, auto /*x*/, auto const& edge) {
            auto const& n = edge.normal;
            auto const uu = avg(u);
            auto const vv = avg(v);
            return -eps * (dot(grad(uu), n) * vv.val + dot(grad(vv), n) * uu.val);
        };
        assemble_facets(mesh.boundary_facets(), U, V, quad, B, lhs_sym);

        auto const rhs_form = [=](auto v, auto x, auto const& edge) {
            return a_bd(problem.u_with_grad(x), v, edge);
        };
        assemble_rhs(mesh.boundary_facets(), V, quad, rhs, rhs_form);

        auto const rhs_sym = [=](auto v, auto x, auto const& edge) {
            auto const& n = edge.normal;
            return -eps * dot(grad(v), n) * problem.u(x);
        };
        assemble_rhs(mesh.boundary_facets(), V, quad, rhs, rhs_sym);
    }
    auto t_after_integration = std::chrono::steady_clock::now();

    mat.mumpsify(linear_system);
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(linear_system);
    auto t_after_solver = std::chrono::steady_clock::now();
    postprocess(F.data() + N, U, problem);

    report_cost(solver);
    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Integration:  {:>8%Q %q}\n", as_ms(t_after_integration - t_before_integration));
    fmt::print("Solver    :   {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
}

struct adaptive_igrm_args : igrm_args {
    int steps;
    bool strong_test = false;
};

auto parse_args_adaptive_igrm(int argc, char* argv[]) {
    auto args = adaptive_igrm_args{};

    bool show_help = false;
    auto const cli = lyra::help(show_help)                                                      //
                   | lyra::opt(args.eps, "epsilon")["--eps"]("diffusion coefficient")           //
                   | lyra::opt(args.gamma, "gamma")["--gamma"]("stabilization parameter")       //
                   | lyra::opt(args.nx, "nx")["--nx"]("mesh size (x direction)")                //
                   | lyra::opt(args.ny, "ny")["--ny"]("mesh size (y direction)")                //
                   | lyra::opt(args.p, "p")["-p"]("order of trial functions")                   //
                   | lyra::opt(args.c, "c")["-c"]("continuity of trial functions")              //
                   | lyra::opt(args.P, "P")["-P"]("order of test functions")                    //
                   | lyra::opt(args.C, "C")["-C"]("continuity of test functions")               //
                   | lyra::opt(args.weak_bc)["--weak-bc"]("should BC be imposed weakly")        //
                   | lyra::opt(args.eta, "eta")["--eta"]("weak BC parameter")                   //
                   | lyra::opt(args.steps, "steps")["--steps"]("adaptation steps")              //
                   | lyra::opt(args.strong_test)["--strong-test"]("remove boundary test DOFS")  //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

auto igrm_with_refinement_step(adaptive_igrm_args const& args, ads::partition const& xs,
                               ads::partition const& ys, ads::mumps::solver& solver)
    -> ads::partition {
    // auto const elems_x = args.nx;
    auto const strong_test = args.strong_test;
    // auto const elems_y = args.ny;
    auto const p = args.p;
    auto const c = args.c;
    auto const eps = args.eps;
    auto const gamma = args.gamma;
    auto const eta = args.eta;
    auto const beta = ads::point_t{1, 0};
    auto const problem = eriksson_johnson{eps};

    auto const p_test = args.P;
    auto const c_test = args.C;

    // auto const xs = ads::evenly_spaced(0.0, 1.0, elems_x);
    // auto const ys = ads::evenly_spaced(0.0, 1.0, elems_y);

    auto const bx = ads::make_bspline_basis(xs, p, c);
    auto const by = ads::make_bspline_basis(ys, p, c);

    auto const Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto const By = ads::make_bspline_basis(ys, p_test, c_test);

    auto const mesh = ads::regular_mesh{xs, ys};
    auto const quad = ads::quadrature{&mesh, std::max(p_test + 1, 2)};

    auto const U = ads::space{&mesh, bx, by};
    auto const V = ads::space{&mesh, Bx, By};

    auto const n = U.dof_count();
    auto const N = V.dof_count();
    fmt::print("Dimension: {} + {} = {}\n", n, N, n + N);

    auto F = std::vector<double>(N + n);
    auto linear_system = ads::mumps::problem{F.data(), n + N};

    auto mat = ads::horrible_sparse_matrix{};
    auto test_mat = ads::horrible_sparse_matrix{};

    auto G = [&mat, &test_mat](int row, int col, double val) {
        mat(row, col) += val;
        test_mat(row, col) += val;
    };
    auto B = [N, &mat](int row, int col, double val) {
        mat(row, N + col) += val;
        mat(N + col, row) += val;
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot, ads::grad;

    auto const area = [](auto const& element) {
        return 2 * (length(element.span_x) + length(element.span_y));
    };
    auto const volume = [](auto const& element) {
        return length(element.span_x) * length(element.span_y);
    };
    auto const coeff = [&](auto const& maybe_element) {
        if (maybe_element) {
            auto const& e = maybe_element.value();
            return area(e) / volume(e);
        } else {
            return 0.0;
        }
    };

    auto const eta_e = [&](auto const& edge) {
        auto const d = 2;
        auto const g = (p + 1) * (p + d) / static_cast<double>(d);

        auto const a = coeff(edge.element1);
        auto const b = coeff(edge.element2);

        if (a > 0 && b > 0) {
            return g * (a + b) / 2;
        } else {
            return g * (a + b);
        }
    };

    auto const der_jump_form = [&](auto /*u*/, auto /*v*/, auto u_nder, auto v_nder, auto /*x*/,
                                   const auto& edge) {
        auto const& n = edge.normal;
        auto const eta = eta_e(edge);
        auto const hp1 = std::pow(1.0 / eta, 2 * c + 1);
        auto const hp2 = std::pow(1.0 / eta, 2 * c + 2);
        auto const bn = std::abs(dot(beta, n));
        auto const jumps = jump(u_nder) * jump(v_nder);
        return (0.5 * hp2 * bn + eps * hp1) * jumps;
    };

    auto t_before_integration = std::chrono::steady_clock::now();

    // Test scalar product
    auto const test_product_int_form = [=](auto u, auto v, auto /*x*/) {
        auto const beta_inf = 1;
        auto const L = 1;
        return eps * dot(grad(u), grad(v)) + beta_inf / L * u.val * v.val;
    };
    assemble(V, quad, G, test_product_int_form);

    assemble_facets(mesh.interior_facets(), V, quad, c_test + 1, G, der_jump_form);

    // Bilinear form of the problem
    assemble(U, V, quad, B, [eps, beta](auto u, auto v, auto /*x*/) {
        return eps * dot(grad(u), grad(v)) + dot(grad(u), beta) * v.val;
    });
    assemble_facets(mesh.interior_facets(), U, V, quad, c_test + 1, B, der_jump_form);

    // For discontinuous test only
    auto const missing_form = [=](auto u, auto v, auto /*x*/, auto const& edge) {
        auto const& n = edge.normal;
        auto const uu = avg(u);
        auto const vv = jump(v);
        return -eps * dot(grad(uu), n) * vv.val;
    };
    assemble_facets(mesh.interior_facets(), U, V, quad, B, missing_form);

    // RHS
    assemble_rhs(V, quad, rhs, [problem](auto v, auto x) { return problem.f(x) * v.val; });

    // Boundary forms

    auto const test_prod_bd_form = [=](auto u, auto v, auto /*x*/, const auto& edge) {
        auto const& n = edge.normal;
        auto const eta = eta_e(edge);
        auto const hp1 = 1.0 / eta;
        auto const hp2 = 1.0;
        auto const bn = std::abs(dot(beta, n));
        auto const jumps = u.val * v.val;
        return (0.5 * hp2 * bn + eps * hp1) * jumps;
    };
    auto const lhs_form = [=](auto u, auto v, auto x, auto const& edge) {
        return test_prod_bd_form(avg(u), avg(v), x, edge);
    };

    if (args.weak_bc || !strong_test) {
        assemble_facets(mesh.boundary_facets(), V, quad, G, lhs_form);
        assemble_facets(mesh.boundary_facets(), U, V, quad, B, lhs_form);

        auto const rhs_form = [=](auto v, auto x, auto const& edge) {
            return test_prod_bd_form(problem.u_with_grad(x), v, x, edge);
        };
        assemble_rhs(mesh.boundary_facets(), V, quad, rhs, rhs_form);

        auto const lhs_sym = [=](auto u, auto v, auto /*x*/, auto const& edge) {
            auto const& n = edge.normal;
            auto const uu = avg(u);
            auto const vv = avg(v);
            return -eps * (dot(grad(uu), n) * vv.val + dot(grad(vv), n) * uu.val);
        };
        assemble_facets(mesh.boundary_facets(), U, V, quad, B, lhs_sym);

        auto const rhs_sym = [=](auto v, auto x, auto const& edge) {
            auto const& n = edge.normal;
            return -eps * dot(grad(v), n) * problem.u(x);
        };
        assemble_rhs(mesh.boundary_facets(), V, quad, rhs, rhs_sym);
    }
    // Strong BC
    if (!args.weak_bc) {
        auto is_fixed_test = mark_boundary_dofs(V);
        auto is_fixed_trial = mark_boundary_dofs(U);

        fmt::print("Applying BC\n");
        for (auto const e : mesh.facets()) {
            // Trial space DOFS
            for (auto const i : U.dofs_on_facet(e)) {
                auto const I = U.global_index(i);
                if (is_fixed_trial[I] == 1) {
                    for (auto const j : V.dofs_on_facet(e)) {
                        auto const J = V.global_index(j);
                        mat(N + I, J) = 0;
                    }
                    for (auto const j : U.dofs_on_facet(e)) {
                        auto const J = U.global_index(j);
                        mat(N + I, N + J) = 0;
                    }
                    mat(N + I, N + I) = 1;
                    F[N + I] = 0;
                }
            }
            // Test space DOFS
            if (strong_test) {
                for (auto const i : V.dofs_on_facet(e)) {
                    auto const I = V.global_index(i);
                    if (is_fixed_test[I] == 1) {
                        for (auto const j : V.dofs_on_facet(e)) {
                            auto const J = V.global_index(j);
                            mat(I, J) = 0;
                            test_mat(I, J) = 0;
                        }
                        for (auto const j : U.dofs_on_facet(e)) {
                            auto const J = U.global_index(j);
                            mat(I, N + J) = 0;
                            mat(N + J, I) = 0;
                        }
                        mat(I, I) = 1;
                        F[I] = 0;
                    }
                }
            }

            // Non-zero Dirichlet BC
            auto dimy = ads::dimension{by, p + 1, 1};
            dimy.factorize_matrix();
            auto buf = ads::lin::vector{{dimy.dofs()}};
            ads::compute_projection(buf, dimy.basis, [problem](auto t) {
                return problem.u({0, t});
            });
            ads::lin::solve_with_factorized(dimy.M, buf, dimy.ctx);
            for (int i = 0; i < dimy.dofs(); ++i) {
                F[N + U.global_index({0, i})] = buf(i);
            }
        }
    }

    auto t_after_integration = std::chrono::steady_clock::now();

    mat.mumpsify(linear_system);
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(linear_system);
    auto t_after_solver = std::chrono::steady_clock::now();

    postprocess(F.data() + N, U, problem);
    // Error in the test space norm
    report_error(
        F.data() + N, mesh, quad, U, problem,
        [=](auto const& v) {  //
            return test_product_int_form(v, v, 0);
        },
        [=](auto const& v, auto x, auto const& facet) {
            return test_prod_bd_form(v, v, x, facet);
        });

    auto const res_norm = norm_using_matrix(F.data(), test_mat);
    fmt::print("Residual norm:        {}\n", res_norm);

    // Compute element contributions to the residual norm
    auto elem_contribs =
        std::unordered_map<ads::regular_mesh::element_index, double, element_idx_hasher>{};
    auto facet_contribs = std::unordered_map<ads::regular_mesh::facet_index, double,
                                             facet_idx_hasher, facet_idx_equality>{};

    auto const apply_elem = [&](auto e, auto v) { elem_contribs[e] = v; };
    auto const apply_facet = [&](auto e, auto v) { facet_contribs[e] = v; };

    local_contribution(F.data(), V, quad, apply_elem,
                       [=](auto v, auto x) { return test_product_int_form(v, v, x); });

    local_contribution_facets(F.data(), mesh.interior_facets(), V, quad, c_test + 1, apply_facet,
                              [=](auto v, auto der, auto x, auto const& facet) {
                                  return der_jump_form(v, v, der, der, x, facet);
                              });

    if (args.weak_bc || !strong_test) {
        local_contribution_facets(F.data(), mesh.boundary_facets(), V, quad, c_test + 1,
                                  apply_facet,
                                  [=](auto v, auto /*der*/, auto x, auto const& facet) {
                                      return lhs_form(v, v, x, facet);
                                  });
    }

    // Include edge contributions
    for (auto const& [idx, v] : facet_contribs) {
        auto const [ix, iy, dir] = idx;
        if (dir == ads::orientation::vertical) {
            if (ix == 0) {
                elem_contribs[{ix, iy}] += v;
            } else if (ix == ads::as_signed(xs.size()) - 1) {
                elem_contribs[{ix - 1, iy}] += v;
            } else {
                elem_contribs[{ix - 1, iy}] += 0.5 * v;
                elem_contribs[{ix, iy}] += 0.5 * v;
            }
        } else {
            if (iy == 0) {
                elem_contribs[{ix, iy}] += v;
            } else if (iy == ads::as_signed(ys.size()) - 1) {
                elem_contribs[{ix, iy - 1}] += v;
            } else {
                elem_contribs[{ix, iy - 1}] += 0.5 * v;
                elem_contribs[{ix, iy}] += 0.5 * v;
            }
        }
    }

    // check
    auto total = 0.0;
    for (auto const& [idx, v] : elem_contribs) {
        total += v;
    }
    fmt::print("Sum of contributions: {}\n", std::sqrt(total));

    // Sum contributions vertically
    auto slice_contribs = std::vector<double>(xs.size() - 1);
    for (auto const& [idx, v] : elem_contribs) {
        auto const [ix, iy] = idx;
        slice_contribs[ix] += v;
    }

    auto const marker = mark_for_refinement(slice_contribs, 0.5);
    auto const new_mesh = refine_single_dimension(xs, marker);

    report_cost(solver);
    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Integration:  {:>8%Q %q}\n", as_ms(t_after_integration - t_before_integration));
    fmt::print("Solver    :   {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));

    return new_mesh;
}

auto igrm_with_refinement(int argc, char* argv[]) -> void {
    auto const args = parse_args_adaptive_igrm(argc, argv);
    auto const elems_x = args.nx;
    auto const elems_y = args.ny;
    auto solver = ads::mumps::solver{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems_x);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems_y);

    for (int i = 0; i < args.steps; ++i) {
        fmt::print("After {} refinement steps\n", i);
        xs = igrm_with_refinement_step(args, xs, ys, solver);
        ys = ensure_aspect_ratio(xs, ys, 8);
        fmt::print("Y mesh size: {}\n", ys.size() - 1);

        auto const p = args.p;
        auto const c = args.c;
        auto const p_test = args.P;
        auto const c_test = args.C;

        auto const bx = ads::make_bspline_basis(xs, p, c);
        auto const by = ads::make_bspline_basis(ys, p, c);

        auto const Bx = ads::make_bspline_basis(xs, p_test, c_test);
        auto const By = ads::make_bspline_basis(ys, p_test, c_test);

        auto const mesh = ads::regular_mesh{xs, ys};

        auto const U = ads::space{&mesh, bx, by};
        auto const V = ads::space{&mesh, Bx, By};

        auto const n = U.dof_count();
        auto const N = V.dof_count();
        fmt::print("Dimension: {} + {} = {}\n", n, N, n + N);

        auto dim = n + N;
        if (dim > 100000) {
            fmt::print("Space too large ({}), aborting\n", dim);
            break;
        }
    }
}

auto main(int argc, char* argv[]) -> int {
    // galerkin(argc, argv);
    // igrm(argc, argv);
    igrm_with_refinement(argc, argv);
}
