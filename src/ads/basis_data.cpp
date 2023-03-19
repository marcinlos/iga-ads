// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/basis_data.hpp"

#include "ads/quad/gauss.hpp"
#include "ads/util.hpp"

namespace ads {

basis_data::~basis_data() {
    for (int e = 0; e < elements; ++e) {
        for (int k = 0; k < quad_order; ++k) {
            for (int d = 0; d <= derivatives; ++d) {
                delete[] b[e][k][d];
            }
            delete[] b[e][k];
        }
        delete[] b[e];
        delete[] x[e];
    }
    delete[] b;
    delete[] J;
    delete[] x;
}

basis_data::basis_data(const basis_data& other)
: first_dofs{other.first_dofs}
, element_ranges{other.element_ranges}
, degree{other.degree}
, elements{other.elements}
, dofs{other.dofs}
, derivatives{other.derivatives}
, quad_order{other.quad_order}
, elem_division{other.elem_division}
, points{other.points}
, basis{other.basis}
, b{new double***[elements]}
, x{new double*[elements]}
, w{other.w}
, J{new double[elements]} {
    for (int e = 0; e < elements; ++e) {
        x[e] = new double[quad_order];
        b[e] = new double**[quad_order];

        for (int k = 0; k < quad_order; ++k) {
            x[e][k] = other.x[e][k];
            b[e][k] = new double*[derivatives + 1];

            for (int d = 0; d <= derivatives; ++d) {
                b[e][k][d] = new double[degree + 1];

                for (int i = 0; i < degree + 1; ++i) {
                    b[e][k][d][i] = other.b[e][k][d][i];
                }
            }
        }

        J[e] = other.J[e];
    }
}

basis_data::basis_data(bspline::basis basis, int derivatives, int quad_order, int elem_division)
: first_dofs(bspline::first_nonzero_dofs(basis))
, element_ranges(bspline::elements_supporting_dofs(basis))
, degree(basis.degree)
, elements(basis.elements() * elem_division)
, dofs(basis.dofs())
, derivatives(derivatives)
, quad_order(quad_order)
, elem_division(elem_division)
, points(elements + 1)
, basis(std::move(basis))
, w(quad::gauss::Ws[quad_order]) {
    int p = basis.degree;
    int q = quad_order;
    x = new double*[elements];
    J = new double[elements];
    b = new double***[elements];

    bspline::eval_ctx ctx(p);

    // compute points of the subdivided elements
    for (int e = 0; e < this->basis.elements(); ++e) {
        double x1 = this->basis.points[e];
        double x2 = this->basis.points[e + 1];
        for (int k = 0; k <= elem_division; ++k) {
            points[e * elem_division + k] = ads::lerp(k, elem_division, x1, x2);
        }
    }

    for (int e = 0; e < elements; ++e) {
        x[e] = new double[q];
        double x1 = points[e];
        double x2 = points[e + 1];
        J[e] = 0.5 * (x2 - x1);

        for (int k = 0; k < q; ++k) {
            double tx = ads::quad::gauss::Xs[q][k];
            double t = 0.5 * (tx + 1);
            x[e][k] = ads::lerp(t, x1, x2);
        }

        b[e] = new double**[quad_order];
        for (int k = 0; k < q; ++k) {
            b[e][k] = new double*[derivatives + 1];
            for (int d = 0; d <= derivatives; ++d) {
                b[e][k][d] = new double[p + 1];
            }
            int span = find_span(x[e][k], this->basis);
            eval_basis_with_derivatives(span, x[e][k], this->basis, b[e][k], derivatives, ctx);
        }
    }
}

}  // namespace ads
