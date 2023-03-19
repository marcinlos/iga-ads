// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_FORM_MATRIX_HPP
#define ADS_FORM_MATRIX_HPP

#include "ads/basis_data.hpp"
#include "ads/lin/band_matrix.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/util/function_value/function_value_1d.hpp"

namespace ads {

void gram_matrix_1d(lin::band_matrix& M, const basis_data& d);

void stiffness_matrix_1d(lin::band_matrix& M, const basis_data& d);

void advection_matrix_1d(lin::band_matrix& M, const basis_data& d);

void gram_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

void stiffness_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

void advection_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

template <typename Form>
void form_matrix(lin::band_matrix& M, const basis_data& d, Form&& form) {
    for (element_id e = 0; e < d.elements; ++e) {
        for (int q = 0; q < d.quad_order; ++q) {
            int first = d.first_dof(e);
            int last = d.last_dof(e);
            for (int a = 0; a + first <= last; ++a) {
                for (int b = 0; b + first <= last; ++b) {
                    int ia = a + first;
                    int ib = b + first;
                    auto va = d.b[e][q][0][a];
                    auto vb = d.b[e][q][0][b];
                    auto da = d.b[e][q][1][a];
                    auto db = d.b[e][q][1][b];
                    auto fa = function_value_1d{va, da};
                    auto fb = function_value_1d{vb, db};
                    M(ia, ib) += form(fb, fa) * d.w[q] * d.J[e];
                }
            }
        }
    }
}

}  // namespace ads

#endif  // ADS_FORM_MATRIX_HPP
