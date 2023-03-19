// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/form_matrix.hpp"

namespace ads {

void gram_matrix_1d(lin::band_matrix& M, const basis_data& d) {
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
                    M(ia, ib) += va * vb * d.w[q] * d.J[e];
                }
            }
        }
    }
}

void stiffness_matrix_1d(lin::band_matrix& M, const basis_data& d) {
    for (element_id e = 0; e < d.elements; ++e) {
        for (int q = 0; q < d.quad_order; ++q) {
            int first = d.first_dof(e);
            int last = d.last_dof(e);
            for (int a = 0; a + first <= last; ++a) {
                for (int b = 0; b + first <= last; ++b) {
                    int ia = a + first;
                    int ib = b + first;
                    auto da = d.b[e][q][1][a];
                    auto db = d.b[e][q][1][b];
                    M(ia, ib) += da * db * d.w[q] * d.J[e];
                }
            }
        }
    }
}

void advection_matrix_1d(lin::band_matrix& M, const basis_data& d) {
    for (element_id e = 0; e < d.elements; ++e) {
        for (int q = 0; q < d.quad_order; ++q) {
            int first = d.first_dof(e);
            int last = d.last_dof(e);
            for (int a = 0; a + first <= last; ++a) {
                for (int b = 0; b + first <= last; ++b) {
                    int ia = a + first;
                    int ib = b + first;
                    auto va = d.b[e][q][0][a];
                    auto db = d.b[e][q][1][b];
                    M(ia, ib) += va * db * d.w[q] * d.J[e];
                }
            }
        }
    }
}

void gram_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V) {
    for (element_id e = 0; e < V.elements; ++e) {
        for (int q = 0; q < V.quad_order; ++q) {
            for (int a = 0; a + V.first_dof(e) <= V.last_dof(e); ++a) {
                for (int b = 0; b + U.first_dof(e) <= U.last_dof(e); ++b) {
                    int ia = a + V.first_dof(e);
                    int ib = b + U.first_dof(e);
                    auto va = V.b[e][q][0][a];
                    auto vb = U.b[e][q][0][b];
                    M(ia, ib) += va * vb * V.w[q] * V.J[e];
                }
            }
        }
    }
}

void stiffness_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V) {
    for (element_id e = 0; e < V.elements; ++e) {
        for (int q = 0; q < V.quad_order; ++q) {
            for (int a = 0; a + V.first_dof(e) <= V.last_dof(e); ++a) {
                for (int b = 0; b + U.first_dof(e) <= U.last_dof(e); ++b) {
                    int ia = a + V.first_dof(e);
                    int ib = b + U.first_dof(e);
                    auto da = V.b[e][q][1][a];
                    auto db = U.b[e][q][1][b];
                    M(ia, ib) += da * db * V.w[q] * V.J[e];
                }
            }
        }
    }
}

void advection_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V) {
    for (element_id e = 0; e < V.elements; ++e) {
        for (int q = 0; q < V.quad_order; ++q) {
            for (int a = 0; a + V.first_dof(e) <= V.last_dof(e); ++a) {
                for (int b = 0; b + U.first_dof(e) <= U.last_dof(e); ++b) {
                    int ia = a + V.first_dof(e);
                    int ib = b + U.first_dof(e);
                    auto va = V.b[e][q][0][a];
                    auto db = U.b[e][q][1][b];
                    M(ia, ib) += va * db * V.w[q] * V.J[e];
                }
            }
        }
    }
}

}  // namespace ads
