// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/os.h>
#include <lyra/lyra.hpp>

#include "ads/experimental/all.hpp"
#include "ads/experimental/horrible_sparse_matrix.hpp"
#include "ads/experimental/space_factory.hpp"

template <typename U>
auto save_heat_to_file(std::string const& path, double time, U const& u) -> void {
    constexpr auto res = 50;
    auto extent = fmt::format("0 {0} 0 {0} 0 {0}", res);
    auto spacing = fmt::format("{0} {0} {0}", 1.0 / res);

    auto out = fmt::output_file(path);
    out.print("<?xml version=\"1.0\"?>\n");
    out.print("<VTKFile type=\"ImageData\" version=\"0.1\">\n");
    out.print("  <ImageData WholeExtent=\"{}\" Origin=\"0 0 0\" Spacing=\"{}\">\n", extent,
              spacing);
    out.print("    <Piece Extent=\"{}\">\n", extent);
    out.print("      <PointData Scalars=\"u\">\n");

    out.print("        <DataArray Name=\"u\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"1\">\n");
    for (auto t : ads::evenly_spaced(0.0, time, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
                const auto X = ads::point3_t{x, y, t};
                out.print("{:.7}\n", u(X));
            }
        }
    }
    out.print("        </DataArray>\n");
    out.print("      </PointData>\n");
    out.print("    </Piece>\n");
    out.print("  </ImageData>\n");
    out.print("</VTKFile>\n");
}

constexpr double pi = M_PI;

class heat_transfer {
public:
    using point = ads::point_t;

    auto u(point p, double t) const noexcept -> double {
        auto const [x, y] = p;
        auto const scale = std::exp(-2 * pi * pi * t);
        auto const shape = std::sin(pi * x) * std::sin(pi * y);
        return scale * shape;
    }

    auto u(double t) const noexcept {
        return [this, t](point p) { return u(p, t); };
    };
};

auto galerkin_main(int /*argc*/, char* /*argv*/[]) -> void {
    auto const elems = 32;
    auto const p = 2;
    auto const c = 1;
    auto const T = 1;
    auto const eps = 5 * 1.0e-3;
    auto const s = 0 * 1.0;
    auto const beta_x = s * 0.0;
    auto const beta_y = s * 1.0;

    auto const xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto const ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto const ts = ads::evenly_spaced(0.0, T, 32);
    // auto const ts = ads::evenly_spaced(0.0, T, 32);

    auto const bx = ads::make_bspline_basis(xs, p, c);
    auto const by = ads::make_bspline_basis(ys, p, c);
    auto const bt = ads::make_bspline_basis(ts, p, c);

    auto const mesh = ads::regular_mesh3{xs, ys, ts};
    auto const quad = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto spaces = ads::space_factory{};

    auto const U = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const Vx = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const Vy = spaces.next<ads::space3>(&mesh, bx, by, bt);

    auto const n = spaces.dim();
    fmt::print("Dimension: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto mat = ads::horrible_sparse_matrix{};
    auto M = [&mat](int row, int col, double val) { mat(row, col) += val; };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    fmt::print("Assembling matrix\n");
    assemble(U, quad, M, [](auto u, auto v, auto /*x*/) { return u.dz * v.val; });
    assemble(Vx, quad, M, [](auto sx, auto tx, auto /*x*/) { return sx.val * tx.val; });
    assemble(Vy, quad, M, [](auto sy, auto ty, auto /*x*/) { return sy.val * ty.val; });
    assemble(U, Vx, quad, M,
             [=](auto u, auto tx, auto /*x*/) { return (eps * u.dx - beta_x * u.val) * tx.val; });
    assemble(U, Vy, quad, M,
             [=](auto u, auto ty, auto /*x*/) { return (eps * u.dy - beta_y * u.val) * ty.val; });
    assemble(Vx, U, quad, M, [](auto sx, auto v, auto /*x*/) { return sx.dx * v.val; });
    assemble(Vy, U, quad, M, [](auto sy, auto v, auto /*x*/) { return sy.dy * v.val; });

    fmt::print("Assembling RHS\n");
    assemble_rhs(U, quad, rhs, [](auto v, auto /*x*/) { return 0 * v.val; });

    fmt::print("Collecting BC\n");
    auto is_fixed = std::vector<int>(n);
    for (auto const dof : U.dofs()) {
        auto const [ix, iy, it] = dof;
        bool fix = false;

        // initial condition t = 0
        fix = fix || it == 0;

        // spatial boundary
        fix = fix || ix == 0 || ix == U.space_x().dof_count() - 1;
        fix = fix || iy == 0 || iy == U.space_y().dof_count() - 1;

        if (fix) {
            is_fixed[U.global_index(dof)] = 1;
        }
    }

    fmt::print("Applying BC\n");
    for (auto const e : mesh.elements()) {
        for (auto const i : U.dofs(e)) {
            auto const I = U.global_index(i);
            if (is_fixed[I] == 1) {
                for (auto const j : U.dofs(e)) {
                    auto const J = U.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : Vx.dofs(e)) {
                    auto const J = Vx.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : Vy.dofs(e)) {
                    auto const J = Vy.global_index(j);
                    mat(I, J) = 0;
                }
                mat(I, I) = 1;
                F[I] = 0;
            }
        }
    }

    F[U.global_index({elems / 2, elems / 2, 0})] = 1;
    F[U.global_index({elems / 2 + 1, elems / 2, 0})] = 1;
    F[U.global_index({elems / 2 + 1, elems / 2 + 1, 0})] = 1;
    F[U.global_index({elems / 2, elems / 2 + 1, 0})] = 1;

    mat.mumpsify(problem);
    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());

    fmt::print("Solving\n");
    solver.solve(problem);

    auto u = ads::bspline_function3(&U, F.data());
    fmt::print("Saving\n");
    save_heat_to_file("dup-full.vti", T, u);
    fmt::print("Done\n");
}

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

struct stabilized_args {
    double eps = 1e-1;
    double beta_x = 0;
    double beta_y = 0;
    double T = 1e-2;
    int nx = 16;
    int ny = 16;
    int nt = 16;
    int p = 2;
    int c = 1;
    int P = 2;
    int C = 1;
};

auto parse_args(int argc, char* argv[]) {
    auto args = stabilized_args{};

    bool show_help = false;
    auto const cli = lyra::help(show_help)                                             //
                   | lyra::opt(args.eps, "epsilon")["--eps"]("diffusion coefficient")  //
                   | lyra::opt(args.beta_x, "beta_x")["--bx"]("x advection")           //
                   | lyra::opt(args.beta_y, "beta_y")["--by"]("y advection")           //
                   | lyra::opt(args.T, "T")["--T"]("total time")                       //
                   | lyra::opt(args.nx, "nx")["--nx"]("mesh size (x direction)")       //
                   | lyra::opt(args.ny, "ny")["--ny"]("mesh size (y direction)")       //
                   | lyra::opt(args.nt, "nt")["--nt"]("mesh size (time direction)")    //
                   | lyra::opt(args.p, "p")["-p"]("order of trial functions")          //
                   | lyra::opt(args.c, "c")["-c"]("continuity of trial functions")     //
                   | lyra::opt(args.P, "P")["-P"]("order of test functions")           //
                   | lyra::opt(args.C, "C")["-C"]("continuity of test functions")      //
        ;

    auto const result = cli.parse({argc, argv});
    validate_args(cli, result, show_help);
    return args;
}

auto stabilized_main(int argc, char* argv[]) -> void {
    auto const args = parse_args(argc, argv);

    auto const elems_x = args.nx;
    auto const elems_y = args.ny;
    auto const elems_t = args.nt;
    auto const p = args.p;
    auto const c = args.c;
    auto const p_test = args.P;
    auto const c_test = args.C;
    auto const T = args.T;
    auto const eps = args.eps;
    auto const beta_x = args.beta_x;
    auto const beta_y = args.beta_y;

    auto const xs = ads::evenly_spaced(0.0, 1.0, elems_x);
    auto const ys = ads::evenly_spaced(0.0, 1.0, elems_y);
    auto const ts = ads::evenly_spaced(0.0, T, elems_t);

    auto const mesh = ads::regular_mesh3{xs, ys, ts};
    auto const quad = ads::quadrature3{&mesh, std::max(p_test + 1, 2)};

    auto const bx = ads::make_bspline_basis(xs, p, c);
    auto const by = ads::make_bspline_basis(ys, p, c);
    auto const bt = ads::make_bspline_basis(ts, p, c);

    auto const Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto const By = ads::make_bspline_basis(ys, p_test, c_test);
    auto const Bt = ads::make_bspline_basis(ts, p_test, c_test);

    auto spaces = ads::space_factory{};

    auto const U = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const Vt = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const Vx = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const Vy = spaces.next<ads::space3>(&mesh, bx, by, bt);
    auto const L = spaces.next<ads::space3>(&mesh, Bx, By, Bt);

    auto const n = spaces.dim();
    fmt::print("Dimension: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto mat = ads::horrible_sparse_matrix{};
    auto M = [&mat](int row, int col, double val) { mat(row, col) += val; };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    fmt::print("Assembling matrix\n");
    assemble(U, quad, M, [=](auto u, auto v, auto /*x*/) {
        return u.val * v.val                                                    //
             + (-eps * u.dx + beta_x * u.val) * (-eps * v.dx + beta_x * v.val)  //
             + (-eps * u.dy + beta_y * u.val) * (-eps * v.dy + beta_y * v.val);
    });
    assemble(Vt, quad, M, [](auto st, auto tt, auto /*x*/) { return st.val * tt.val; });
    assemble(Vx, quad, M, [](auto sx, auto tx, auto /*x*/) { return sx.val * tx.val; });
    assemble(Vy, quad, M, [](auto sy, auto ty, auto /*x*/) { return sy.val * ty.val; });

    assemble(U, Vt, quad, M, [=](auto u, auto tt, auto /*x*/) { return -u.val * tt.val; });
    assemble(Vt, U, quad, M, [=](auto st, auto v, auto /*x*/) { return -st.val * v.val; });

    assemble(U, Vx, quad, M,
             [=](auto u, auto tx, auto /*x*/) { return -(-eps * u.dx + beta_x * u.val) * tx.val; });
    assemble(Vx, U, quad, M,
             [=](auto sx, auto v, auto /*x*/) { return -(-eps * v.dx + beta_x * v.val) * sx.val; });

    assemble(U, Vy, quad, M,
             [=](auto u, auto ty, auto /*x*/) { return -(-eps * u.dy + beta_y * u.val) * ty.val; });
    assemble(Vy, U, quad, M,
             [=](auto sy, auto v, auto /*x*/) { return -(-eps * v.dy + beta_y * v.val) * sy.val; });

    assemble(L, Vt, quad, M, [](auto l, auto vt, auto /*x*/) { return l.val * vt.dz; });
    assemble(L, Vx, quad, M, [](auto l, auto vx, auto /*x*/) { return l.val * vx.dx; });
    assemble(L, Vy, quad, M, [](auto l, auto vy, auto /*x*/) { return l.val * vy.dy; });

    assemble(Vt, L, quad, M, [](auto st, auto w, auto /*x*/) { return w.val * st.dz; });
    assemble(Vx, L, quad, M, [](auto sx, auto w, auto /*x*/) { return w.val * sx.dx; });
    assemble(Vy, L, quad, M, [](auto sy, auto w, auto /*x*/) { return w.val * sy.dy; });

    fmt::print("Assembling RHS\n");
    assemble_rhs(U, quad, rhs, [](auto v, auto /*x*/) { return 0 * v.val; });

    fmt::print("Collecting BC\n");
    auto is_fixed = std::vector<int>(n);
    for (auto const dof : U.dofs()) {
        auto const [ix, iy, it] = dof;
        bool fix = false;

        // initial condition t = 0
        fix = fix || it == 0;

        // spatial boundary
        fix = fix || ix == 0 || ix == U.space_x().dof_count() - 1;
        fix = fix || iy == 0 || iy == U.space_y().dof_count() - 1;

        if (fix) {
            is_fixed[U.global_index(dof)] = 1;
        }
    }

    fmt::print("Applying BC\n");
    for (auto const e : mesh.elements()) {
        for (auto const i : U.dofs(e)) {
            auto const I = U.global_index(i);
            if (is_fixed[I] == 1) {
                for (auto const j : U.dofs(e)) {
                    auto const J = U.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : Vt.dofs(e)) {
                    auto const J = Vt.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : Vx.dofs(e)) {
                    auto const J = Vx.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : Vy.dofs(e)) {
                    auto const J = Vy.global_index(j);
                    mat(I, J) = 0;
                }
                for (auto const j : L.dofs(e)) {
                    auto const J = L.global_index(j);
                    mat(I, J) = 0;
                }
                mat(I, I) = 1;
                F[I] = 0;
            }
        }
    }

    // assemble_rhs(L, quad, rhs, [](auto v, auto p) {
    //     auto const [x, y, t] = p;
    //     return v.val * std::sin(x) * std::sin(y);
    // });
    F[U.global_index({elems_x / 2, elems_y / 2, 0})] = 1;
    F[U.global_index({elems_x / 2 + 1, elems_y / 2, 0})] = 1;
    F[U.global_index({elems_x / 2 + 1, elems_y / 2 + 1, 0})] = 1;
    F[U.global_index({elems_x / 2, elems_y / 2 + 1, 0})] = 1;

    mat.mumpsify(problem);
    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());

    fmt::print("Saving matrix to file\n");
    solver.save_to_file(problem, "matrix");

    fmt::print("Solving\n");
    solver.solve(problem);

    auto u = ads::bspline_function3(&U, F.data());
    fmt::print("Saving\n");
    save_heat_to_file("dup-full.vti", T, u);
    fmt::print("Done\n");
}

auto main(int argc, char* argv[]) -> int {
    // galerkin_main(argc, argv);
    stabilized_main(argc, argv);
}
