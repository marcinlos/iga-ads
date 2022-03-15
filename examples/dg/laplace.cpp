// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/experimental/all.hpp"

constexpr double pi = M_PI;
using std::cos;
using std::sin;

/////////////

template <typename Concrete>
class poisson_base {
private:
    auto self() const -> const Concrete* { return static_cast<const Concrete*>(this); }

public:
    auto u() const noexcept {
        return [this](auto x) { return self()->u(x); };
    }

    auto f() const noexcept {
        return [this](auto x) { return self()->f(x); };
    }

    auto g() const noexcept {
        return [this](auto x) { return self()->g(x); };
    }
};

class poisson_type1 : private poisson_base<poisson_type1> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        const auto v = sin(pi * x) * sin(pi * y);
        return v * v;
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        const auto sx = sin(pi * x);
        const auto c2x = cos(2 * pi * x);
        const auto sy = sin(pi * y);
        const auto c2y = cos(2 * pi * y);
        return -2 * pi * pi * (c2x * sy * sy + c2y * sx * sx);
    }

    auto g(point /*X*/) const noexcept -> double { return 0.0; }
};

class poisson_type2 : private poisson_base<poisson_type2> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        return 1 + sin(pi * x) * sin(pi * y);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        return 2 * pi * pi * sin(pi * x) * sin(pi * y);
    }

    auto g(point /*X*/) const noexcept -> double { return 1.0; }
};

class poisson_type3 : private poisson_base<poisson_type3> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        return x * x + 0.5 * y * y + sin(pi * x) * sin(pi * y);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        return -3 + 2 * pi * pi * sin(pi * x) * sin(pi * y);
    }

    auto g(point X) const noexcept -> double {
        const auto [x, y] = X;
        return x * x + 0.5 * y * y;
    }
};

/////////////

void DG_poisson();
void DG_stokes();
void DGiGRM_stokes();

void poisson_3D();
void DG_poisson_3D();
void DG_stokes_3D();
void DGiGRM_stokes_3D();

int main() {
    try {
        // DG_poisson();
        // DG_stokes();
        // DGiGRM_stokes();

        // poisson_3D();
        // DG_poisson_3D();
        // DG_stokes_3D();
        DGiGRM_stokes_3D();
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        std::exit(1);
    }
}

void DG_poisson() {
    auto elems = 128;
    auto p = 3;
    auto c = -1;
    // auto eta = 1000000.0;  // good for Nitsche BC
    auto eta = 10.0;

    auto poisson = poisson_type1{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);

    auto mesh = ads::regular_mesh{xs, ys};
    auto space = ads::space{&mesh, bx, by};
    auto quad = ads::quadrature{&mesh, p + 1};

    auto n = space.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto out = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int J, double val) { F[J] += val; };
    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(space, quad, out, [](auto u, auto v, auto /*x*/) { return dot(grad(u), grad(v)); });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    auto form = [eta](auto u, auto v, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto h = length(edge.span);
        // clang-format off
        return - dot(grad(avg(v)), n) * jump(u).val
               - dot(grad(avg(u)), n) * jump(v).val
               + eta / h * jump(u).val * jump(v).val;
        // clang-format on
    };
    assemble_facets(mesh.facets(), space, quad, out, form);
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(space, quad, rhs, [&poisson](auto v, auto x) {  //
        return v.val * poisson.f(x);
    });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    auto bd_form = [eta, &poisson](auto v, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto h = length(edge.span);
        const auto g = poisson.g(x);
        // clang-format off
        return - dot(grad(v), n) * g
               + eta/h * g * v.val;
        // clang-format on
    };
    assemble_rhs(mesh.boundary_facets(), space, quad, rhs, bd_form);
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    fmt::print("Computing error\n");

    auto u = ads::bspline_function(&space, F.data());

    auto t_before_err = std::chrono::steady_clock::now();
    auto err = error(mesh, quad, L2{}, u, poisson.u());
    auto t_after_err = std::chrono::steady_clock::now();

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", u);
    auto t_after_output = std::chrono::steady_clock::now();

    fmt::print("error = {:.6}\n", err);

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

template <typename Concrete>
class stokes_base {
private:
    auto self() const -> const Concrete* { return static_cast<const Concrete*>(this); }

public:
    auto vx() const noexcept {
        return [this](auto x) { return self()->vx(x); };
    }

    auto vy() const noexcept {
        return [this](auto x) { return self()->vy(x); };
    }

    auto vz() const noexcept {
        return [this](auto x) { return self()->vz(x); };
    }

    auto p(double avg = 0) const noexcept {
        return [this, avg](auto x) { return avg + self()->p(x); };
    }

    auto fx() const noexcept {
        return [this](auto x) { return self()->fx(x); };
    }

    auto fy() const noexcept {
        return [this](auto x) { return self()->fy(x); };
    }

    auto fz() const noexcept {
        return [this](auto x) { return self()->fz(x); };
    }
};

class stokes_polynomial : private stokes_base<stokes_polynomial> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p(point p) const noexcept -> double {
        const auto [x, y] = p;
        return x * (1 - x) - 1. / 6;
    }

    auto vx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return x * x * (1 - x) * (1 - x) * (2 * y - 6 * y * y + 4 * y * y * y);
    }

    auto vy(point p) const noexcept -> double {
        const auto [x, y] = p;
        return -y * y * (1 - y) * (1 - y) * (2 * x - 6 * x * x + 4 * x * x * x);
    }

    auto fx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return (12 - 24 * y) * x * x * x * x                         //
             + (-24 + 48 * y) * x * x * x                            //
             + (-48 * y + 72 * y * y - 48 * y * y * y + 12) * x * x  //
             + (-2 + 24 * y - 72 * y * y + 48 * y * y * y) * x       //
             + 1 - 4 * y + 12 * y * y - 8 * y * y * y;
    }

    auto fy(point p) const noexcept -> double {
        const auto [x, y] = p;
        return (8 - 48 * y + 48 * y * y) * x * x * x                                //
             + (-12 + 72 * y - 72 * y * y) * x * x                                  //
             + (4 - 24 * y + 48 * y * y - 48 * y * y * y + 24 * y * y * y * y) * x  //
             - 12 * y * y + 24 * y * y * y - 12 * y * y * y * y;
    }
};

class stokes_nonpoly : private stokes_base<stokes_nonpoly> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto xx = x * x;
        const auto yy = y * y;
        const auto ex = exp(x);
        const auto e = exp(1);

        const auto a = 456                               //
                     + xx * (228 - 5 * (yy - y))         //
                     + 2 * x * (-228 + (yy - y))         //
                     + 2 * pow(x, 3) * (-36 + (yy - y))  //
                     + pow(x, 4) * (12 + (yy - y));

        return -424 + 156 * e + (yy - y) * (-456 + ex * a);
    };

    auto vx(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto ex = exp(x);

        return 2 * ex * pow(-1 + x, 2) * x * x * (y * y - y) * (-1 + 2 * y);
    }

    auto vy(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto ex = exp(x);

        return -ex * (-1 + x) * x * (-2 + x * (3 + x)) * pow(-1 + y, 2) * y * y;
    }

    auto fx(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto x2 = x * x;
        const auto x3 = pow(x, 3);
        const auto x4 = pow(x, 4);
        const auto y2 = y * y;
        const auto y3 = pow(x, y);
        const auto ex = exp(x);

        const auto px = ex * (y - 1) * y
                      * (                            //
                            x4 * (y2 - y + 12)       //
                            + 6 * x3 * (y2 - y - 4)  //
                            + x2 * (y2 - y + 12)     //
                            - 8 * x * (y - 1) * y    //
                            + 2 * (y - 1) * y        //
                      );

        const auto Lux = 2 * ex
                       * (                                             //
                             x4 * (2 * y3 - 3 * y2 + 13 * y - 6)       //
                             + 6 * x3 * (2 * y3 - 3 * y2 - 3 * y + 2)  //
                             + x2 * (2 * y3 - 3 * y2 + 13 * y - 6)     //
                             - 8 * x * y * (2 * y2 - 3 * y + 1)        //
                             + 2 * y * (2 * y2 - 3 * y + 1)            //
                       );

        return -Lux + px;
    }

    auto fy(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto x2 = x * x;
        const auto x3 = pow(x, 3);
        const auto x4 = pow(x, 4);
        const auto y2 = y * y;
        const auto y3 = pow(x, y);
        const auto y4 = pow(x, y);
        const auto ex = exp(x);

        const auto py = 2 * (2 * y - 1)
                      * (ex
                             * (                                 //
                                 x4 * (y2 - y + 6)               //
                                 + 2 * x3 * (y2 - y - 18)        //
                                 + x2 * (-5 * y2 + 5 * y + 114)  //
                                 + 2 * x * (y2 - y - 114) + 228  //
                                 )
                         - 228);

        const auto Luy = -ex
                       * (                                                       //
                           x4 * (y4 - 2 * y3 + 13 * y2 - 12 * y + 2)             //
                           + 2 * x3 * (5 * y4 - 10 * y3 + 17 * y2 - 12 * y + 2)  //
                           + x2 * (19 * y4 - 38 * y3 - 41 * y2 + 60 * y - 10)    //
                           + x * (-6 * y4 + 12 * y3 + 18 * y2 - 24 * y + 4)      //
                           - 6 * pow(y - 1, 2) * y2                              //
                       );

        return -Luy + py;
    }
};

class stokes_cavity : private stokes_base<stokes_cavity> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p([[maybe_unused]] point p) const noexcept -> double { return 0.0; }

    auto vx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return std::abs(y - 1) < 1e-5 ? 1.0 : 0.0;
    }

    auto vy([[maybe_unused]] point p) const noexcept -> double { return 0.0; }

    auto fx([[maybe_unused]] point p) const noexcept -> double { return 0.0; }

    auto fy([[maybe_unused]] point p) const noexcept -> double { return 0.0; }
};

class space_factory {
private:
    ads::global_dof offset_ = 0;

public:
    template <typename Space, typename... Args>
    auto next(Args&&... args) -> Space {
        auto space = Space{std::forward<Args>(args)..., offset_};
        offset_ += space.dof_count();
        return space;
    }
};

void DG_stokes() {
    auto elems = 16;
    auto p = 4;
    auto c = -1;
    auto eta = 10.0 * (p + 1) * (p + 2);

    auto stokes = stokes_cavity{};
    // auto stokes = stokes_polynomial{};
    // auto stokes = stokes_nonpoly{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);

    auto mesh = ads::regular_mesh{xs, ys};
    auto quad = ads::quadrature{&mesh, std::max(p + 1, 2)};

    auto spaces = space_factory{};

    auto Vx = spaces.next<ads::space>(&mesh, bx, by);
    auto Vy = spaces.next<ads::space>(&mesh, bx, by);
    auto P = spaces.next<ads::space>(&mesh, bx, by);

    auto n = Vx.dof_count() + Vy.dof_count() + P.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto M = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    // clang-format off
    assemble(Vx,    quad, M, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy,    quad, M, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(P, Vx, quad, M, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P, Vy, quad, M, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(Vx, P, quad, M, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy, P, quad, M, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    // clang-format on
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.facets(), Vx, quad, M, [eta](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, quad, M, [eta](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), P, Vx, quad, M, [](auto p, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{jump(vx).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vy, quad, M, [](auto p, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{0, jump(vy).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, P, quad, M, [](auto ux, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{jump(ux).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, P, quad, M, [](auto uy, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{0, jump(uy).val};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.interior_facets(), P, quad, M, [](auto p, auto q, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return h * jump(p).val * jump(q).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Vx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Vy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), Vx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function(&Vx, F.data());
    auto sol_vy = ads::bspline_function(&Vy, F.data());
    auto sol_p = ads::bspline_function(&P, F.data());

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_p = error(mesh, quad, L2{}, sol_p, stokes.p(mean));
    auto err = sum_norms(err_vx, err_vy, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    const auto vx_J =
        integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            const auto d = jump(sol_vx(x, edge));
            return 1 / h * d * d;
        });
    const auto vy_J =
        integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            const auto d = jump(sol_vy(x, edge));
            return 1 / h * d * d;
        });
    const auto r_P =
        integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            auto d = jump(sol_p(x, edge));
            return h * d * d;
        });
    fmt::print("vx seminorm = {:.6}\n", std::sqrt(vx_J));
    fmt::print("vy seminorm = {:.6}\n", std::sqrt(vy_J));
    fmt::print("p  seminorm = {:.6}\n", std::sqrt(r_P));

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", sol_vx, sol_vy, sol_p);
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

void DGiGRM_stokes() {
    auto elems = 16;
    auto p = 4;
    auto eta = 1.0 * (p + 1) * (p + 2);

    // auto stokes = stokes_cavity{};
    // auto stokes = stokes_polynomial{};
    auto stokes = stokes_nonpoly{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto mesh = ads::regular_mesh{xs, ys};
    auto quad = ads::quadrature{&mesh, std::max(p + 1, 2)};

    // Test
    auto p_test = 4;
    auto c_test = -1;

    auto Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto By = ads::make_bspline_basis(xs, p_test, c_test);

    auto tests = space_factory{};

    auto Wx = tests.next<ads::space>(&mesh, Bx, By);
    auto Wy = tests.next<ads::space>(&mesh, Bx, By);
    auto Q = tests.next<ads::space>(&mesh, Bx, By);

    auto N = Wx.dof_count() + Wy.dof_count() + Q.dof_count();
    fmt::print("Test  DoFs: {:10L}\n", N);

    // Trial
    auto p_trial = 4;
    auto c_trial = 3;  // >= 0

    auto bx = ads::make_bspline_basis(xs, p_trial, c_trial);
    auto by = ads::make_bspline_basis(ys, p_trial, c_trial);

    auto trials = space_factory{};

    auto Vx = trials.next<ads::space>(&mesh, bx, by);
    auto Vy = trials.next<ads::space>(&mesh, bx, by);
    auto P = trials.next<ads::space>(&mesh, bx, by);

    auto n = Vx.dof_count() + Vy.dof_count() + P.dof_count();
    fmt::print("Trial DoFs: {:10L}\n", n);
    fmt::print("Total:      {:10L}\n", N + n);

    auto F = std::vector<double>(N + n);
    auto problem = ads::mumps::problem{F};
    auto solver = ads::mumps::solver{};

    auto G = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto B = [&problem, N](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, N + col + 1, val);
            problem.add(N + col + 1, row + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    // clang-format off
    assemble(Wx, quad, G, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Wy, quad, G, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Q,  quad, G, [](auto p,  auto q,  auto /*x*/) { return p.val * q.val;           });

    assemble(Vx, Wx, quad, B, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy, Wy, quad, B, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(P,  Wx, quad, B, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P,  Wy, quad, B, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(Vx,  Q, quad, B, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy,  Q, quad, B, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    // clang-format on
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.facets(), Wx, quad, G, [](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return 1/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Wy, quad, G, [](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return 1/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.interior_facets(), Q, quad, G, [](auto p, auto q, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return h * jump(p).val * jump(q).val;
    });

    assemble_facets(mesh.facets(), Vx, Wx, quad, B, [eta](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, Wy, quad, B, [eta](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), P, Wx, quad, B, [](auto p, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{jump(vx).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wy, quad, B, [](auto p, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{0, jump(vy).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.boundary_facets(), Vx, Q, quad, B, [](auto ux, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{jump(ux).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.boundary_facets(), Vy, Q, quad, B, [](auto uy, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{0, jump(uy).val};
        return - dot(u, n) * avg(q).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Wx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Wy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), Wx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function(&Vx, F.data() + N);
    auto sol_vy = ads::bspline_function(&Vy, F.data() + N);
    auto sol_p = ads::bspline_function(&P, F.data() + N);

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_p = error(mesh, quad, L2{}, sol_p, stokes.p(mean));
    auto err = sum_norms(err_vx, err_vy, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto r_vx = ads::bspline_function(&Wx, F.data());
    auto r_vy = ads::bspline_function(&Wy, F.data());
    auto r_p = ads::bspline_function(&Q, F.data());

    auto norm_r_vx = norm(mesh, quad, L2{}, r_vx);
    auto J_r_vx =
        std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            auto d = jump(r_vx(x, edge));
            return 1 / h * d * d;
        }));
    auto norm_r_vy = norm(mesh, quad, L2{}, r_vy);
    auto J_r_vy =
        std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            auto d = jump(r_vy(x, edge));
            return 1 / h * d * d;
        }));
    auto norm_r_p = norm(mesh, quad, L2{}, r_p);
    auto q_r_p = std::sqrt(
        integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
            const auto h = length(edge.span);
            auto d = jump(r_p(x, edge));
            return h * d * d;
        }));
    auto res_norm = sum_norms(norm_r_vx, J_r_vx, norm_r_vy, J_r_vy, norm_r_p, q_r_p);
    fmt::print("||r_vx|| = {:.6}\n", norm_r_vx);
    fmt::print("||r_vy|| = {:.6}\n", norm_r_vy);
    fmt::print("||r_p||  = {:.6}\n", norm_r_p);
    fmt::print("|r_vx|   = {:.6}\n", J_r_vx);
    fmt::print("|r_vy|   = {:.6}\n", J_r_vy);
    fmt::print("|r_p|    = {:.6}\n", q_r_p);
    fmt::print("res norm = {:.6}\n", res_norm);

    auto val = alt_norm(mesh, quad, L2_2{}, [](auto X) {
        auto [x, y] = X;
        return ads::value_type{x * x - y * y, 2 * x, -2 * y};
        // return 1;
        // return ads::facet_value<double>{0, 1};
    });
    fmt::print("H1 norm : {}\n", val);

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", sol_vx, sol_vy, sol_p, r_vx, r_vy, r_p);
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

class poisson3_type1 : private poisson_base<poisson3_type1> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto g([[maybe_unused]] point X) const noexcept -> double { return 0.0; }
};

class poisson3_type2 : private poisson_base<poisson3_type2> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x * x + 0.5 * y * y + 0.3 * z * z + sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return -3.6 + 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto g(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x * x + 0.5 * y * y + 0.3 * z * z;
    }
};

void poisson_3D() {
    auto elems = 32;
    auto p = 2;
    auto c = p - 1;
    auto eta = 10.0;

    auto poisson = poisson3_type2{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);
    auto bz = ads::make_bspline_basis(zs, p, c);

    auto mesh = ads::regular_mesh3{xs, ys, zs};
    auto quad = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto space = ads::space3{&mesh, bx, by, bz};

    auto n = space.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto out = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int J, double val) { F[J] += val; };
    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(space, quad, out, [](auto u, auto v, auto /*x*/) { return dot(grad(u), grad(v)); });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.boundary_facets(), space, quad, out, [eta](auto u, auto v, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(v)), n) * jump(u).val
               - dot(grad(avg(u)), n) * jump(v).val
               + eta/h * jump(u).val * jump(v).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(space, quad, rhs, [&poisson](auto v, auto x) { return v.val * poisson.f(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), space, quad, rhs, [eta,&poisson](auto v, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = poisson.g(x);
        return - dot(grad(v), n) * g
               + eta/h * g * v.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    fmt::print("Computing error\n");

    auto u = ads::bspline_function3(&space, F.data());

    auto t_before_err = std::chrono::steady_clock::now();
    auto err = error(mesh, quad, L2{}, u, poisson.u());
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("error = {:.6}\n", err);

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    // fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

void DG_poisson_3D() {
    auto elems = 8;
    auto p = 4;
    auto c = -1;
    auto eta = 10.0 * (p + 1) * (p + 2);

    auto poisson = poisson3_type2{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);
    auto bz = ads::make_bspline_basis(zs, p, c);

    auto mesh = ads::regular_mesh3{xs, ys, zs};
    auto quad = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto space = ads::space3{&mesh, bx, by, bz};

    auto n = space.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto out = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int J, double val) { F[J] += val; };
    using ads::dot;
    using ads::grad;

    // auto executor = ads::galois_executor{12};
    // auto executor = ads::sequential_executor{};

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(space, quad, out, [](auto u, auto v, auto /*x*/) { return dot(grad(u), grad(v)); });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.facets(), space, quad, out, [eta](auto u, auto v, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(v)), n) * jump(u).val
               - dot(grad(avg(u)), n) * jump(v).val
               + eta/h * jump(u).val * jump(v).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(space, quad, rhs, [&poisson](auto v, auto x) { return v.val * poisson.f(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), space, quad, rhs, [eta,&poisson](auto v, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = poisson.g(x);
        return - dot(grad(v), n) * g
               + eta/h * g * v.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    fmt::print("Computing error\n");

    auto u = ads::bspline_function3(&space, F.data());

    auto t_before_err = std::chrono::steady_clock::now();
    auto err = error(mesh, quad, L2{}, u, poisson.u());
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("error = {:.6}\n", err);

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    // fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

class stokes3_type1 : private stokes_base<stokes3_type1> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::p, base::vx, base::vy, base::vz, base::fx, base::fy, base::fz;

    auto p(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x * (1 - x) - 1. / 6;
    }

    auto vx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x * x * (1 - x) * (1 - x) * (2 * y - 6 * y * y + 4 * y * y * y);
    }

    auto vy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return -y * y * (1 - y) * (1 - y) * (2 * x - 6 * x * x + 4 * x * x * x);
    }

    auto vz(point /*X*/) const noexcept -> double { return 0.0; }

    auto fx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return (12 - 24 * y) * x * x * x * x                         //
             + (-24 + 48 * y) * x * x * x                            //
             + (-48 * y + 72 * y * y - 48 * y * y * y + 12) * x * x  //
             + (-2 + 24 * y - 72 * y * y + 48 * y * y * y) * x       //
             + 1 - 4 * y + 12 * y * y - 8 * y * y * y;
    }

    auto fy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return (8 - 48 * y + 48 * y * y) * x * x * x                                //
             + (-12 + 72 * y - 72 * y * y) * x * x                                  //
             + (4 - 24 * y + 48 * y * y - 48 * y * y * y + 24 * y * y * y * y) * x  //
             - 12 * y * y + 24 * y * y * y - 12 * y * y * y * y;
    }

    auto fz(point /*X*/) const noexcept -> double { return 0.0; }
};

class stokes3_type2 : private stokes_base<stokes3_type2> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::p, base::vx, base::vy, base::vz, base::fx, base::fy, base::fz;

    auto p(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        return x2 * (2 * y - 1) * (2 * y - 1) * (x + y - z) * (x + y - z) - 5. / 54;
    }

    auto vx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto xm12 = (1 - x) * (1 - x);
        const auto ym12 = (1 - y) * (1 - y);
        const auto zm12 = (1 - z) * (1 - z);
        return x2 * y2 * xm12 * (2 * y - 2) + 2 * x2 * y * xm12 * ym12
             - x2 * z2 * xm12 * (2 * z - 2) - 2 * x2 * z * xm12 * zm12;
    }

    auto vy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto xm12 = (1 - x) * (1 - x);
        const auto ym12 = (1 - y) * (1 - y);
        const auto zm12 = (1 - z) * (1 - z);
        return -x2 * y2 * ym12 * (2 * x - 2) - 2 * x * y2 * xm12 * ym12
             + y2 * z2 * ym12 * (2 * z - 2) + 2 * y2 * z * ym12 * zm12;
    }

    auto vz(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto xm12 = (1 - x) * (1 - x);
        const auto ym12 = (1 - y) * (1 - y);
        const auto zm12 = (1 - z) * (1 - z);
        return x2 * z2 * zm12 * (2 * x - 2) + 2 * x * z2 * xm12 * zm12
             - y2 * z2 * zm12 * (2 * y - 2) - 2 * y * z2 * ym12 * zm12;
    }

    auto fx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto x3 = x2 * x;
        const auto y3 = y2 * y;
        const auto z3 = z2 * z;
        const auto x4 = x3 * x;
        return -24 * x4 * y + 24 * x4 * z + 48 * x3 * y - 48 * x3 * z - 48 * x2 * y3 + 72 * x2 * y2
             - 48 * x2 * y + 48 * x2 * z3 - 72 * x2 * z2 + 48 * x2 * z
             + x2 * (2 * y - 1) * (2 * y - 1) * (2 * x + 2 * y - 2 * z) + 48 * x * y3 - 72 * x * y2
             + 24 * x * y - 48 * x * z3 + 72 * x * z2 - 24 * x * z
             + 2 * x * (2 * y - 1) * (2 * y - 1) * (x + y - z) * (x + y - z) - 8 * y3 + 12 * y2
             - 4 * y + 8 * z3 - 12 * z2 + 4 * z;
    }

    auto fy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto x3 = x2 * x;
        const auto y3 = y2 * y;
        const auto z3 = z2 * z;
        const auto y4 = y3 * y;
        return 48 * x3 * y2 - 48 * x3 * y + 8 * x3 - 72 * x2 * y2 + 72 * x2 * y
             + x2 * (2 * y - 1) * (2 * y - 1) * (2 * x + 2 * y - 2 * z)
             + x2 * (8 * y - 4) * (x + y - z) * (x + y - z) - 12 * x2 + 24 * x * y4 - 48 * x * y3
             + 48 * x * y2 - 24 * x * y + 4 * x - 24 * y4 * z + 48 * y3 * z - 48 * y2 * z3
             + 72 * y2 * z2 - 48 * y2 * z + 48 * y * z3 - 72 * y * z2 + 24 * y * z - 8 * z3
             + 12 * z2 - 4 * z;
    }

    auto fz(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto z2 = z * z;
        const auto x3 = x2 * x;
        const auto y3 = y2 * y;
        const auto z3 = z2 * z;
        const auto z4 = z3 * z;
        return -48 * x3 * z2 + 48 * x3 * z - 8 * x3 + 72 * x2 * z2 - 72 * x2 * z
             + x2 * (2 * y - 1) * (2 * y - 1) * (-2 * x - 2 * y + 2 * z) + 12 * x2 - 24 * x * z4
             + 48 * x * z3 - 48 * x * z2 + 24 * x * z - 4 * x + 48 * y3 * z2 - 48 * y3 * z + 8 * y3
             - 72 * y2 * z2 + 72 * y2 * z - 12 * y2 + 24 * y * z4 - 48 * y * z3 + 48 * y * z2
             - 24 * y * z + 4 * y;
    }
};

class stokes3_cavity : private stokes_base<stokes3_cavity> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::p, base::vx, base::vy, base::vz, base::fx, base::fy, base::fz;

    auto p(point /*X*/) const noexcept -> double { return 0.0; }

    auto vx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return std::abs(y - 1) < 1e-5 ? 1.0 : 0.0;
        // if (std::abs(y - 1) < 1e-5) {
        //     return sin(pi * x) * sin(pi * z);
        // } else {
        //     return 0.0;
        // }
    }

    auto vy(point /*X*/) const noexcept -> double { return 0.0; }

    auto vz(point /*X*/) const noexcept -> double { return 0.0; }

    auto fx(point /*X*/) const noexcept -> double { return 0.0; }

    auto fy(point /*X*/) const noexcept -> double { return 0.0; }

    auto fz(point /*X*/) const noexcept -> double { return 0.0; }
};

template <typename Vx, typename Vy, typename Vz, typename P>
auto save_to_file3(const std::string& path, Vx&& vx, Vy&& vy, Vz&& vz, P&& pressure) -> void {
    constexpr auto res = 50;
    auto extent = fmt::format("0 {0} 0 {0} 0 {0}", res);

    auto out = fmt::output_file(path);
    out.print("<?xml version=\"1.0\"?>\n");
    out.print("<VTKFile type=\"ImageData\" version=\"0.1\">\n");
    out.print("  <ImageData WholeExtent=\"{}\" origin=\"0 0 0\" spacing=\"1 1 1\">\n", extent);
    out.print("    <Piece Extent=\"{}\">\n", extent);
    out.print("      <PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n", extent);

    out.print("        <DataArray Name=\"Velocity\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"3\">\n");
    for (auto z : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
                const auto X = ads::point3_t{x, y, z};
                out.print("{:.7} {:.7} {:.7}\n", vx(X), vy(X), vz(X));
            }
        }
    }
    out.print("        </DataArray>\n");

    out.print("        <DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" "
              "NumberOfComponents=\"1\">\n");
    for (auto z : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
                const auto X = ads::point3_t{x, y, z};
                out.print("{:.7}\n", pressure(X));
            }
        }
    }
    out.print("        </DataArray>\n");

    out.print("      </PointData>\n");
    out.print("    </Piece>\n");
    out.print("  </ImageData>\n");
    out.print("</VTKFile>\n");
}

void DG_stokes_3D() {
    auto elems = 20;
    auto p = 1;
    auto c = -1;
    auto eta = 10.0 * (p + 1) * (p + 2);

    // auto stokes = stokes3_type1{};
    auto stokes = stokes3_cavity{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);
    auto bz = ads::make_bspline_basis(zs, p, c);

    auto mesh = ads::regular_mesh3{xs, ys, zs};
    auto quad = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto spaces = space_factory{};

    auto Vx = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto Vy = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto Vz = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto P = spaces.next<ads::space3>(&mesh, bx, by, bz);

    auto n = Vx.dof_count() + Vy.dof_count() + Vz.dof_count() + P.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F = std::vector<double>(n);
    auto problem = ads::mumps::problem{F.data(), n};
    auto solver = ads::mumps::solver{};

    auto M = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    // clang-format off
    assemble(Vx,    quad, M, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy,    quad, M, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Vz,    quad, M, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(P, Vx, quad, M, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P, Vy, quad, M, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(P, Vz, quad, M, [](auto p,  auto vz, auto /*x*/) { return - p.val * vz.dz;         });
    assemble(Vx, P, quad, M, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy, P, quad, M, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    assemble(Vz, P, quad, M, [](auto uz, auto  q, auto /*x*/) { return   uz.dz * q.val;         });
    // clang-format on
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.facets(), Vx, quad, M, [eta](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, quad, M, [eta](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Vz, quad, M, [eta](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vz)), n) * jump(uz).val
               - dot(grad(avg(uz)), n) * jump(vz).val
               + eta/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.facets(), P, Vx, quad, M, [](auto p, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{jump(vx).val, 0, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vy, quad, M, [](auto p, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, jump(vy).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vz, quad, M, [](auto p, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, 0, jump(vz).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, P, quad, M, [](auto ux, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{jump(ux).val, 0, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, P, quad, M, [](auto uy, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, jump(uy).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vz, P, quad, M, [](auto uz, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, 0, jump(uz).val};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.interior_facets(), P, quad, M, [](auto p, auto q, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return h * jump(p).val * jump(q).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Vx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Vy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    assemble_rhs(Vz, quad, rhs, [&stokes](auto vz, auto x) { return vz.val * stokes.fz(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), Vx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vz, quad, rhs, [eta,&stokes](auto vz, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vz(x);
        return - dot(grad(vz), n) * g
               + eta/h * g * vz.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function3(&Vx, F.data());
    auto sol_vy = ads::bspline_function3(&Vy, F.data());
    auto sol_vz = ads::bspline_function3(&Vz, F.data());
    auto sol_p = ads::bspline_function3(&P, F.data());

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_vz = error(mesh, quad, L2{}, sol_vz, stokes.vz());
    auto err_p = error(mesh, quad, L2{}, sol_p, stokes.p(mean));
    auto err = sum_norms(err_vx, err_vy, err_vz, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("vz error = {:.6}\n", err_vz);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto vf = integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(sol_p(x, face));
        return h * d * d;
    });
    fmt::print("Pressure jump seminorm = {:.6}\n", std::sqrt(vf));

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file3("result.vti", sol_vx, sol_vy, sol_vz, sol_p);
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}

void DGiGRM_stokes_3D() {
    // n=8, (2,-1) (2,1) works well
    auto elems = 4;
    auto p = 2;
    // auto eta = 10.0 * (p + 1) * (p + 2);
    auto eta = 1.0 * (p + 1) * (p + 2);

    auto stokes = stokes3_type2{};
    // auto stokes = stokes3_cavity{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto mesh = ads::regular_mesh3{xs, ys, zs};
    auto quad = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    // Test
    auto p_test = p;
    auto c_test = -1;

    auto Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto By = ads::make_bspline_basis(ys, p_test, c_test);
    auto Bz = ads::make_bspline_basis(zs, p_test, c_test);

    auto tests = space_factory{};

    auto Wx = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Wy = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Wz = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Q = tests.next<ads::space3>(&mesh, Bx, By, Bz);

    auto N = Wx.dof_count() + Wy.dof_count() + Wz.dof_count() + Q.dof_count();
    fmt::print("Test  DoFs: {:10L}\n", N);

    // Trial
    auto p_trial = p;
    auto c_trial = 1;  // >= 0

    auto bx = ads::make_bspline_basis(xs, p_trial, c_trial);
    auto by = ads::make_bspline_basis(ys, p_trial, c_trial);
    auto bz = ads::make_bspline_basis(zs, p_trial, c_trial);

    auto trials = space_factory{};

    auto Vx = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto Vy = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto Vz = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto P = trials.next<ads::space3>(&mesh, bx, by, bz);

    auto n = Vx.dof_count() + Vy.dof_count() + Vz.dof_count() + P.dof_count();
    fmt::print("Trial DoFs: {:10L}\n", n);
    fmt::print("Total:      {:10L}\n", N + n);

    auto F = std::vector<double>(N + n);
    auto problem = ads::mumps::problem{F.data(), N + n};
    auto solver = ads::mumps::solver{};

    auto G = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto B = [&problem, N](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, N + col + 1, val);
            problem.add(N + col + 1, row + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    // clang-format off
    assemble(Wx, quad, G, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Wy, quad, G, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Wz, quad, G, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(Q,  quad, G, [](auto p,  auto q,  auto /*x*/) { return p.val * q.val;           });

    assemble(Vx, Wx, quad, B, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy, Wy, quad, B, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Vz, Wz, quad, B, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(P,  Wx, quad, B, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P,  Wy, quad, B, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(P,  Wz, quad, B, [](auto p,  auto vz, auto /*x*/) { return - p.val * vz.dz;         });
    assemble(Vx,  Q, quad, B, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy,  Q, quad, B, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    assemble(Vz,  Q, quad, B, [](auto uz, auto  q, auto /*x*/) { return   uz.dz * q.val;         });
    // clang-format on
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    // clang-format off
    assemble_facets(mesh.facets(), Wx, quad, G, [](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Wy, quad, G, [](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Wz, quad, G, [](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.interior_facets(), Q, quad, G, [](auto p, auto q, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return h * jump(p).val * jump(q).val;
    });

    assemble_facets(mesh.facets(), Vx, Wx, quad, B, [eta](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, Wy, quad, B, [eta](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Vz, Wz, quad, B, [eta](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vz)), n) * jump(uz).val
               - dot(grad(avg(uz)), n) * jump(vz).val
               + eta/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.facets(), P, Wx, quad, B, [](auto p, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{jump(vx).val, 0, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wy, quad, B, [](auto p, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, jump(vy).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wz, quad, B, [](auto p, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, 0, jump(vz).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, Q, quad, B, [](auto ux, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{jump(ux).val, 0, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, Q, quad, B, [](auto uy, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, jump(uy).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vz, Q, quad, B, [](auto uz, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, 0, jump(uz).val};
        return - dot(u, n) * avg(q).val;
    });
    // clang-format on
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Wx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Wy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    assemble_rhs(Wz, quad, rhs, [&stokes](auto vz, auto x) { return vz.val * stokes.fz(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    // clang-format off
    assemble_rhs(mesh.boundary_facets(), Wx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wz, quad, rhs, [eta,&stokes](auto vz, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vz(x);
        return - dot(grad(vz), n) * g
               + eta/h * g * vz.val;
    });
    // clang-format on
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function3(&Vx, F.data() + N);
    auto sol_vy = ads::bspline_function3(&Vy, F.data() + N);
    auto sol_vz = ads::bspline_function3(&Vz, F.data() + N);
    auto sol_p = ads::bspline_function3(&P, F.data() + N);

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_vz = error(mesh, quad, L2{}, sol_vz, stokes.vz());
    auto err_p = error(mesh, quad, L2{}, sol_p, stokes.p(mean));
    auto err = sum_norms(err_vx, err_vy, err_vz, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("vz error = {:.6}\n", err_vz);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto r_vx = ads::bspline_function3(&Wx, F.data());
    auto r_vy = ads::bspline_function3(&Wy, F.data());
    auto r_vz = ads::bspline_function3(&Wz, F.data());
    auto r_p = ads::bspline_function3(&Q, F.data());

    auto norm_r_vx = norm(mesh, quad, L2{}, r_vx);
    // clang-format off
    auto J_r_vx = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vx(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_vy = norm(mesh, quad, L2{}, r_vy);
    auto J_r_vy = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vy(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_vz = norm(mesh, quad, L2{}, r_vz);
    auto J_r_vz = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vz(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_p = norm(mesh, quad, L2{}, r_p);
    auto q_r_p = std::sqrt(integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_p(x, face));
        return h * d * d;
    }));
    // clang-format on
    auto res_norm =
        sum_norms(norm_r_vx, J_r_vx, norm_r_vy, J_r_vy, norm_r_vz, J_r_vz, norm_r_p, q_r_p);

    fmt::print("||r_vx|| = {:.6}\n", norm_r_vx);
    fmt::print(" |r_vx|  = {:.6}\n", J_r_vx);
    fmt::print("||r_vy|| = {:.6}\n", norm_r_vy);
    fmt::print(" |r_vy|  = {:.6}\n", J_r_vy);
    fmt::print("||r_vz|| = {:.6}\n", norm_r_vz);
    fmt::print(" |r_vz|  = {:.6}\n", J_r_vz);
    fmt::print("||r_p||  = {:.6}\n", norm_r_p);
    fmt::print(" |r_p|   = {:.6}\n", q_r_p);
    fmt::print("res norm = {:.6}\n", res_norm);

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file3("result.vti", sol_vx, sol_vy, sol_vz, [&](auto x) { return sol_p(x) - mean; });
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output - t_before_output));
}
