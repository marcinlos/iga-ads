#ifndef ADS_PROJECTION_HPP_
#define ADS_PROJECTION_HPP_

#include "ads/basis_data.hpp"


namespace ads {


template <typename Rhs, typename Function>
void compute_projection(Rhs& u, const basis_data& d1, const basis_data& d2, const basis_data& d3, Function f) {
    for (element_id e1 = 0; e1 < d1.elements; ++ e1) {
    for (element_id e2 = 0; e2 < d2.elements; ++ e2) {
    for (element_id e3 = 0; e3 < d3.elements; ++ e3) {
        double J = d1.J[e1] * d2.J[e2] * d3.J[e3];
        int first1 = d1.first_dof(e1);
        int last1  = d1.last_dof(e1);
        int first2 = d2.first_dof(e2);
        int last2  = d2.last_dof(e2);
        int first3 = d3.first_dof(e3);
        int last3  = d3.last_dof(e3);

        for (int q1 = 0; q1 < d1.quad_order; ++ q1) {
        for (int q2 = 0; q2 < d2.quad_order; ++ q2) {
        for (int q3 = 0; q3 < d3.quad_order; ++ q3) {
            double w = d1.w[q1] * d2.w[q2] * d3.w[q3];
            double x1 = d1.x[e1][q1];
            double x2 = d2.x[e2][q2];
            double x3 = d3.x[e3][q3];

            for (int a1 = 0; a1 + first1 <= last1; ++ a1) {
            for (int a2 = 0; a2 + first2 <= last2; ++ a2) {
            for (int a3 = 0; a3 + first3 <= last3; ++ a3) {
                int i1 = a1 + first1;
                int i2 = a2 + first2;
                int i3 = a3 + first3;

                double B1 = d1.b[e1][q1][0][a1];
                double B2 = d2.b[e2][q2][0][a2];
                double B3 = d3.b[e3][q3][0][a3];
                double B = B1 * B2 * B3;
                u(i1, i2, i3) += f(x1, x2, x3) * B * w * J;
            }
            }
            }
        }
        }
        }
    }
    }
    }
}


}


#endif /* ADS_PROJECTION_HPP_ */
