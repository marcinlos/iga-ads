#ifndef ADS_BSPLINE_EVAL_HPP_
#define ADS_BSPLINE_EVAL_HPP_

#include "ads/bspline/bspline.hpp"
#include "ads/util/function_value.hpp"

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
function_value_1d eval_ders(double x, const U& u, const basis& b, eval_ders_ctx& ctx) {
    int span = find_span(x, b);
    double** bvals = ctx.basis_vals();
    eval_basis_with_derivatives(span, x, b, bvals, 1, ctx);
    int offset = span - b.degree; // first nonzero function on element

    double value = 0;
    double dx = 0;
    for (int i = 0; i <= b.degree; ++ i) {
        double uu = u(i + offset);
        value += uu * bvals[0][i];
        dx += uu * bvals[1][i];
    }
    return { value, dx };
}

template <typename U>
double eval(double x, double y, const U& u, const basis& bx, const basis& by, eval_ctx& cx,
        eval_ctx& cy) {

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
function_value_2d eval_ders(double x, double y, const U& u, const basis& bx, const basis& by,
        eval_ders_ctx& cx, eval_ders_ctx& cy) {

    int spanx = find_span(x, bx);
    int spany = find_span(y, by);

    double** bvx = cx.basis_vals();
    double** bvy = cy.basis_vals();

    eval_basis_with_derivatives(spanx, x, bx, bvx, 1, cx);
    eval_basis_with_derivatives(spany, y, by, bvy, 1, cy);

    int offsetx = spanx - bx.degree;
    int offsety = spany - by.degree;

    double value = 0;
    double dx = 0;
    double dy = 0;

    for (int ix = 0; ix < bx.dofs_per_element(); ++ ix) {
        for (int iy = 0; iy < by.dofs_per_element(); ++ iy) {
            int ixx = ix + offsetx;
            int iyy = iy + offsety;
            double uu = u(ixx, iyy);
            value += uu * bvx[0][ix] * bvy[0][iy];
            dx += uu * bvx[1][ix] * bvy[0][iy];
            dy += uu * bvx[0][ix] * bvy[1][iy];
        }
    }
    return { value, dx, dy };
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

template <typename U>
function_value_3d eval_ders(double x, double y, double z, const U& u, const basis& bx,
        const basis& by, const basis& bz, eval_ders_ctx& cx, eval_ders_ctx& cy, eval_ders_ctx& cz) {

    int spanx = find_span(x, bx);
    int spany = find_span(y, by);
    int spanz = find_span(z, bz);

    double** bvx = cx.basis_vals();
    double** bvy = cy.basis_vals();
    double** bvz = cz.basis_vals();

    eval_basis_with_derivatives(spanx, x, bx, bvx, 1, cx);
    eval_basis_with_derivatives(spany, y, by, bvy, 1, cy);
    eval_basis_with_derivatives(spanz, z, bz, bvz, 1, cz);

    int offsetx = spanx - bx.degree;
    int offsety = spany - by.degree;
    int offsetz = spanz - bz.degree;

    double value = 0;
    double dx = 0;
    double dy = 0;
    double dz = 0;

    for (int ix = 0; ix < bx.dofs_per_element(); ++ ix) {
        for (int iy = 0; iy < by.dofs_per_element(); ++ iy) {
            for (int iz = 0; iz < bz.dofs_per_element(); ++ iz) {
                int ixx = ix + offsetx;
                int iyy = iy + offsety;
                int izz = iz + offsetz;
                double uu = u(ixx, iyy, izz);
                value += uu * bvx[0][ix] * bvy[0][iy] * bvz[0][iz];
                dx += uu * bvx[1][ix] * bvy[0][iy] * bvz[0][iz];
                dy += uu * bvx[0][ix] * bvy[1][iy] * bvz[0][iz];
                dz += uu * bvx[0][ix] * bvy[0][iy] * bvz[1][iz];
            }
        }
    }
    return { value, dx, dy, dz };
}

}
}

#endif /* ADS_BSPLINE_EVAL_HPP_ */
