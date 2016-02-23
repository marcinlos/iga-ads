#ifndef ADS_BSPLINE_EVAL_HPP_
#define ADS_BSPLINE_EVAL_HPP_

#include "ads/bspline/bspline.hpp"

namespace ads {
namespace bspline {


template <typename U>
double eval(double x, const U& u, const basis& b, eval_ctx& ctx) {
    int span = find_span(x, b);
    double* bvals = ctx.basis_vals();
    eval_basis(span, x, b, bvals, ctx);
    int offset = span - b.degree; // first nonzero function on element

    double value = 0;
    for (int i = 0; i <= b.degree; ++ i) {
        value += u(i + offset) * bvals[i];
    }
    return value;
}


template <typename U>
double eval(double x, double y, const U& u, const basis& bx, const basis& by, eval_ctx& cx, eval_ctx& cy) {

    int spanx = find_span(x, bx);
    int spany = find_span(y, by);

    double* bvx = cx.basis_vals();
    double* bvy = cy.basis_vals();

    eval_basis(spanx, x, bx, bvx, cx);
    eval_basis(spany, y, by, bvy, cy);

    int offsetx = spanx - bx.degree;
    int offsety = spany - by.degree;

    double value = 0;
    for (int ix = 0; ix < bx.dofs_per_element(); ++ ix) {
    for (int iy = 0; iy < by.dofs_per_element(); ++ iy) {
        int ixx = ix + offsetx;
        int iyy = iy + offsety;
        value += u(ixx, iyy) * bvx[ix] * bvy[iy];
    }
    }
    return value;
}


template <typename U>
double eval(double x, double y, double z,
        const U& u,
        const basis& bx, const basis& by, const basis& bz,
        eval_ctx& cx, eval_ctx& cy, eval_ctx& cz) {

    int spanx = find_span(x, bx);
    int spany = find_span(y, by);
    int spanz = find_span(z, bz);

    double* bvx = cx.basis_vals();
    double* bvy = cy.basis_vals();
    double* bvz = cz.basis_vals();

    eval_basis(spanx, x, bx, bvx, cx);
    eval_basis(spany, y, by, bvy, cy);
    eval_basis(spanz, z, bz, bvz, cz);

    int offsetx = spanx - bx.degree;
    int offsety = spany - by.degree;
    int offsetz = spanz - bz.degree;

    double value = 0;
    for (int ix = 0; ix < bx.dofs_per_element(); ++ ix) {
        for (int iy = 0; iy < by.dofs_per_element(); ++ iy) {
            for (int iz = 0; iz < bz.dofs_per_element(); ++ iz) {
                int ixx = ix + offsetx;
                int iyy = iy + offsety;
                int izz = iz + offsetz;
                value += u(ixx, iyy, izz) * bvx[ix] * bvy[iy] * bvz[iz];
            }
        }
    }
    return value;
}



}
}


#endif /* ADS_BSPLINE_EVAL_HPP_ */
