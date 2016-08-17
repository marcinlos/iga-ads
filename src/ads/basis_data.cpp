#include "ads/basis_data.hpp"
#include "ads/quad/gauss.hpp"
#include "ads/util.hpp"


namespace ads {

    basis_data::basis_data(const bspline::basis& basis, int derivatives)
    : first_dofs(bspline::first_nonzero_dofs(basis))
    , degree(basis.degree)
    , elements(basis.elements())
    , quad_order(basis.degree + 1)
    , knot(basis.knot)
    , basis(basis)
    , w(quad::gauss::Ws[quad_order])
    {
        int p = basis.degree;
        int q = quad_order;
        x = new double*[elements];
        J = new double[elements];
        b = new double***[elements];

        bspline::eval_ctx ctx(p);

        for (int e = 0; e < elements; ++ e) {
            x[e] = new double[q];
            double x1 = basis.knot[p + e];
            double x2 = basis.knot[p + e + 1];
            J[e] = 0.5 * (x2 - x1);

            for (int k = 0; k < q; ++ k) {
                double tx = ads::quad::gauss::Xs[q][k];
                double t = 0.5 * (tx + 1);
                x[e][k] = ads::lerp(t, x1, x2);
            }

            b[e] = new double**[quad_order];
            for (int k = 0; k < q; ++ k) {
                b[e][k] = new double*[derivatives + 1];
                for (int d = 0; d <= derivatives; ++ d) {
                    b[e][k][d] = new double[p + 1];
                }
                for (int i = 0; i < basis.dofs_per_element(); ++ i) {
                    eval_basis_with_derivatives(e + p, x[e][k], basis, b[e][k], derivatives, ctx);
                }
            }
        }
    }

}
