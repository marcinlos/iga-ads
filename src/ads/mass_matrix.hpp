#ifndef ADS_MASS_MATRIX_HPP_
#define ADS_MASS_MATRIX_HPP_

#include "ads/lin/band_matrix.hpp"
#include "ads/basis_data.hpp"


namespace ads {

void gram_matrix_1d(lin::band_matrix& M, const basis_data& d) {
    for (element_id e = 0; e < d.elements; ++ e) {
        for (int q = 0; q < d.quad_order; ++ q) {
            int first = d.first_dof(e);
            int last = d.last_dof(e);
            for (int a = 0; a + first <= last; ++ a) {
                for (int b = 0; b + first <= last; ++ b) {
                    int ia = a + first;
                    int ib = b + first;
                    M(ia, ib) += d.b[e][q][0][a] * d.b[e][q][0][b] * d.w[q] * d.J[e];
                }
            }
        }
    }
}


}


#endif /* ADS_MASS_MATRIX_HPP_ */
