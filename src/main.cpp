#include <iostream>

#include "ads/lin/band_matrix.hpp"
#include "ads/basis_data.hpp"
#include "ads/mass_matrix.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/bspline/bspline.hpp"
#include "ads/solver.hpp"
#include "ads/output_manager.hpp"
#include "ads/projection.hpp"
#include "ads/problems/heat.hpp"

using namespace ads;


int main() {
    int p = 2;
    int elements = 20;
    double a = 0, b = 1;
    double dt = 1e-5;

    auto B = bspline::create_basis(a, b, p, elements);
    int n = B.dofs();

    lin::band_matrix M { p, p, n };
    lin::tensor<double, 3> u_prev {{n, n, n}}, u{{n, n, n}}, buffer({n, n, n});

    int order = 1;
    basis_data basis(B, order);

    auto g = [](double x, double y, double z) {
        double dx = x - 0.5;
        double dy = y - 0.5;
        double dz = z - 0.5;
        double r2 = std::min(8 * (dx * dx + dy * dy + dz * dz), 1.0);
        return (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    };

    gram_matrix_1d(M, basis);

    lin::solver_ctx ctx { M };
    lin::factorize(M, ctx);

    dim_data dim_data { M, ctx };

    compute_projection(u, basis, basis, basis, g);
    ads_solve(u, buffer, dim_data, dim_data, dim_data);

    int N = 30;
    ads::output_manager mgr { B, B, B, N };

    int iter_count = 100;

    for (int iter = 0; iter < iter_count; ++ iter) {
        using std::swap;
        swap(u_prev, u);
        std::cout << "Iter " << iter << std::endl;

        heat_rhs(u_prev, u, basis, basis, basis, dt);
        ads_solve(u, buffer, dim_data, dim_data, dim_data);

        if (iter % 10 == 0) {
            mgr.to_file(u, "out_%d.vti", iter);
        }

    }
}
