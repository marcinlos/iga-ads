#include <iostream>

#include "ads/lin/band_matrix.hpp"
#include "ads/basis_data.hpp"
#include "ads/mass_matrix.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/bspline/bspline.hpp"
#include "ads/solver.hpp"
#include "ads/output_manager.hpp"
#include "ads/projection.hpp"
#include "ads/problems/tumor.hpp"

using namespace ads;



void tumor_problem() {
    int px = 2;
    int py = 2;

    int ex = 12;
    int ey = 12;

    double a = 0, b = 1;
    double dt = 1e-7;

    auto Bx = bspline::create_basis(a, b, px, ex);
    auto By = bspline::create_basis(a, b, py, ey);

    int nx = Bx.dofs();
    int ny = By.dofs();

    lin::band_matrix Mx{ px, px, nx };
    lin::band_matrix My{ py, py, ny };
    std::array<std::size_t, 2> ns{ nx, ny };

    lin::tensor<double, 2>
        u_prev{ ns }, n_prev{ ns }, mi_u_prev{ ns }, mi_n_prev{ ns },
        u{ ns }, n{ ns }, mi_u{ ns }, mi_n{ ns },
        buffer{ ns };

    int order = 2;
    basis_data basis_x(Bx, order);
    basis_data basis_y(By, order);

    auto g = [](double x, double y) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(8 * (dx * dx + dy * dy), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    };

    gram_matrix_1d(Mx, basis_x);
    gram_matrix_1d(My, basis_y);

    lin::solver_ctx ctx_x{ Mx }, ctx_y{ My };
    lin::factorize(Mx, ctx_x);
    lin::factorize(My, ctx_y);

    dim_data dim_data_x{ Mx, ctx_x }, dim_data_y{ My, ctx_y };

    problems::tumor::params params;
    params.delta = 0.01;
    params.epsilon = 0.005;
    params.P0 = 0.1;
    params.Gamma = 0.045;
    params.chi0 = 0.05;


    compute_projection(u, basis_x, basis_y, g);
    compute_projection(n, basis_x, basis_y, g);
    compute_projection(mi_u, basis_x, basis_y, g);
    compute_projection(mi_n, basis_x, basis_y, g);


    auto solve_all = [&]() {
        ads_solve(u, buffer, dim_data_x, dim_data_y);
        ads_solve(n, buffer, dim_data_x, dim_data_y);
        ads_solve(mi_u, buffer, dim_data_x, dim_data_y);
        ads_solve(mi_n, buffer, dim_data_x, dim_data_y);
    };

    solve_all();

    auto info = problems::tumor::info2d(u, n, basis_x, basis_y, params);
    std::cout << "Initial energy: " << info.energy << ", mass: " << info.mass << std::endl;

    int N = 100;
    ads::output_manager<2> mgr{ Bx, By, N };
    mgr.to_file(u, "init.vti");

    int iter_count = 100;

    for (int iter = 0; iter < iter_count; ++ iter) {

        using std::swap;
        swap(u_prev, u);
        swap(n_prev, n);
        swap(mi_u_prev, mi_u);
        swap(mi_n_prev, mi_n);

        problems::tumor::rhs2d(
                u_prev, n_prev, mi_u_prev, mi_n_prev,
                u, n, mi_u, mi_n, basis_x, basis_y, params, dt);
        solve_all();

        auto info = problems::tumor::info2d(u, n, basis_x, basis_y, params);
        std::cout << "Iter " << iter << ", energy: " << info.energy
                  << ", mass: " << info.mass << std::endl;

        if (iter % 100 == 0) {
            mgr.to_file(u, "out_%d.vti", iter);
        }
    }
}


int main() {
    tumor_problem();
}
