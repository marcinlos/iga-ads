#include <iostream>

#include "ads/simulation/config.hpp"
#include "ads/simulation/dimension.hpp"
#include "ads/projection.hpp"
#include "ads/solver.hpp"
#include "ads/lin/tensor/for_each.hpp"

using namespace ads;

double f(double x, double y, double z, double w) {
    return 1;
}

int main() {
    constexpr std::size_t N = 4;

    // number of derivatives required, for projection we don't need any
    int ders = 0;
    dim_config dim { 2, 4 };

    // we need a separate object for each dimensions since it includes buffers
    dimension x { dim, ders };
    dimension y { dim, ders };
    dimension z { dim, ders };
    dimension w { dim, ders };

    // factorize once, then just call dgbtrs instead of dgbtrf
    x.factorize_matrix();
    y.factorize_matrix();
    z.factorize_matrix();
    w.factorize_matrix();

    // dimensions of RHS arrays
    std::array<std::size_t, N> shape{ x.dofs(), y.dofs(), z.dofs(), w.dofs() };

    using rhs_type = lin::tensor<double, N>;

    // RHS vector and buffer for solver
    rhs_type u{ shape }, buffer{ shape };

    // compute projection of f on B-spline basis
    projector<N> project{{ x.basis, y.basis, z.basis, w.basis }};
    project(u, f);

    // solve the system
    ads_solve(u, buffer, x.data(), y.data(), z.data(), w.data());

    // check - all values should be 1.0
    lin::for_each_multiindex([&u](int ix, int iy, int iz, int iw) {
        double val = u(ix, iy, iz, iw);
        std::cout << "u(" << ix << ", " << iy << ", " << iz << ", " << iw << ") = " << val << std::endl;
    }, u);
}
