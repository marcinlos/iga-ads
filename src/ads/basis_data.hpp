#ifndef ADS_BASIS_DATA_HPP_
#define ADS_BASIS_DATA_HPP_

#include <boost/range/counting_range.hpp>

#include "ads/bspline/bspline.hpp"
#include "ads/util.hpp"
#include "ads/quad/gauss.hpp"

namespace ads {

typedef int element_id;

struct basis_data {
    using dof_range_type = decltype(boost::counting_range(0, 0));

    std::vector<int> first_dofs;
    int degree;
    int elements;
    int quad_order;
    bspline::knot_vector knot;
    const bspline::basis& basis;
    double**** b;
    double** x;
    const double* w;
    double* J;

    basis_data(const bspline::basis& basis, int derivatives)
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
            J[e] = 0.5 / elements;
            x[e] = new double[q];
            double x1 = basis.knot[p + e];
            double x2 = basis.knot[p + e + 1];

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

    int first_dof(element_id e) const {
        return first_dofs[e];
    }

    int last_dof(element_id e) const {
        return first_dof(e) + degree;
    }

    dof_range_type dof_range(element_id e) const {
        return boost::counting_range(first_dof(e), last_dof(e) + 1);
    }

};

}


#endif /* ADS_BASIS_DATA_HPP_ */
