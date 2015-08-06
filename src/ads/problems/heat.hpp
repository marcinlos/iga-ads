#ifndef ADS_PROBLEMS_HEAT_HPP_
#define ADS_PROBLEMS_HEAT_HPP_

#include "ads/basis_data.hpp"


namespace ads {
namespace problems {
namespace heat {


template <typename Rhs>
void rhs2d(const Rhs& U_prev, Rhs& U, const basis_data& d1, const basis_data& d2, double dt) {
    zero(U);

    for (element_id e1 = 0; e1 < d1.elements; ++ e1) {
    for (element_id e2 = 0; e2 < d2.elements; ++ e2) {
        double J = d1.J[e1] * d2.J[e2];
        int first1 = d1.first_dof(e1);
        int last1 = d1.last_dof(e1);
        int first2 = d2.first_dof(e2);
        int last2 = d2.last_dof(e2);

        for (int q1 = 0; q1 < d1.quad_order; ++ q1) {
        for (int q2 = 0; q2 < d2.quad_order; ++ q2) {
            double w = d1.w[q1] * d2.w[q2];

            for (int a1 = 0; a1 + first1 <= last1; ++ a1) {
            for (int a2 = 0; a2 + first2 <= last2; ++ a2) {
                int i1 = a1 + first1;
                int i2 = a2 + first2;

                double B1 = d1.b[e1][q1][0][a1];
                double B2 = d2.b[e2][q2][0][a2];
                double dB1 = d1.b[e1][q1][1][a1];
                double dB2 = d2.b[e2][q2][1][a2];

                double v = B1 * B2;
                double dxv = dB1 * B2;
                double dyv = B1 * dB2;

                double u = 0;
                double dxu = 0, dyu = 0;
                for (int b1 = 0; b1 + first1 <= last1; ++ b1) {
                for (int b2 = 0; b2 + first2 <= last2; ++ b2) {
                    int j1 = b1 + first1;
                    int j2 = b2 + first2;

                    double B1 = d1.b[e1][q1][0][b1];
                    double B2 = d2.b[e2][q2][0][b2];
                    double dB1 = d1.b[e1][q1][1][b1];
                    double dB2 = d2.b[e2][q2][1][b2];

                    double B = B1 * B2;
                    double dxB = dB1 * B2;
                    double dyB = B1 * dB2;

                    double c = U_prev(j1, j2);
                    u   += c * B;
                    dxu += c * dxB;
                    dyu += c * dyB;
                }
                }

                double gradient_prod = dxu * dxv + dyu * dyv;
                double val = u * v - dt * gradient_prod;
                U(i1, i2) += val * w * J;
            }
            }
        }
        }
    }
    }
}



template <typename Rhs>
double energy2d(Rhs& U, const basis_data& d1, const basis_data& d2) {

    double E = 0;
    for (element_id e1 = 0; e1 < d1.elements; ++ e1) {
    for (element_id e2 = 0; e2 < d2.elements; ++ e2) {
        double J = d1.J[e1] * d2.J[e2];
        int first1 = d1.first_dof(e1);
        int last1 = d1.last_dof(e1);
        int first2 = d2.first_dof(e2);
        int last2 = d2.last_dof(e2);

        for (int q1 = 0; q1 < d1.quad_order; ++ q1) {
        for (int q2 = 0; q2 < d2.quad_order; ++ q2) {
            double w = d1.w[q1] * d2.w[q2];

            double u = 0;
            for (int a1 = 0; a1 + first1 <= last1; ++ a1) {
            for (int a2 = 0; a2 + first2 <= last2; ++ a2) {
                int i1 = a1 + first1;
                int i2 = a2 + first2;

                double B1 = d1.b[e1][q1][0][a1];
                double B2 = d2.b[e2][q2][0][a2];

                double B = B1 * B2;
                double c = U(i1, i2);
                u += c * B;
            }
            }
            E += u * w * J;
        }
        }
    }
    }
    return E;
}




template <typename Rhs>
void rhs3d(const Rhs& U_prev, Rhs& U, const basis_data& d1, const basis_data& d2, const basis_data& d3, double dt) {
    zero(U);

    for (element_id e1 = 0; e1 < d1.elements; ++ e1) {
    for (element_id e2 = 0; e2 < d2.elements; ++ e2) {
    for (element_id e3 = 0; e3 < d3.elements; ++ e3) {
        double J = d1.J[e1] * d2.J[e2] * d3.J[e3];
        int first1 = d1.first_dof(e1);
        int last1 = d1.last_dof(e1);
        int first2 = d2.first_dof(e2);
        int last2 = d2.last_dof(e2);
        int first3 = d3.first_dof(e3);
        int last3 = d3.last_dof(e3);

        for (int q1 = 0; q1 < d1.quad_order; ++ q1) {
        for (int q2 = 0; q2 < d2.quad_order; ++ q2) {
        for (int q3 = 0; q3 < d3.quad_order; ++ q3) {
            double w = d1.w[q1] * d2.w[q2] * d3.w[q3];

            for (int a1 = 0; a1 + first1 <= last1; ++ a1) {
            for (int a2 = 0; a2 + first2 <= last2; ++ a2) {
            for (int a3 = 0; a3 + first3 <= last3; ++ a3) {
                int i1 = a1 + first1;
                int i2 = a2 + first2;
                int i3 = a3 + first3;

                double B1 = d1.b[e1][q1][0][a1];
                double B2 = d2.b[e2][q2][0][a2];
                double B3 = d3.b[e3][q3][0][a3];
                double dB1 = d1.b[e1][q1][1][a1];
                double dB2 = d2.b[e2][q2][1][a2];
                double dB3 = d3.b[e3][q3][1][a3];
                double v = B1 * B2 * B3;
                double dxv = dB1 * B2 * B3;
                double dyv = B1 * dB2 * B3;
                double dzv = B1 * B2 * dB3;

                double u = 0;
                double dxu = 0, dyu = 0, dzu = 0;
                for (int b1 = 0; b1 + first1 <= last1; ++ b1) {
                for (int b2 = 0; b2 + first2 <= last2; ++ b2) {
                for (int b3 = 0; b3 + first3 <= last3; ++ b3) {
                    int j1 = b1 + first1;
                    int j2 = b2 + first2;
                    int j3 = b3 + first3;

                    double B1 = d1.b[e1][q1][0][b1];
                    double B2 = d2.b[e2][q2][0][b2];
                    double B3 = d3.b[e3][q3][0][b3];
                    double dB1 = d1.b[e1][q1][1][b1];
                    double dB2 = d2.b[e2][q2][1][b2];
                    double dB3 = d3.b[e3][q3][1][b3];
                    double B = B1 * B2 * B3;
                    double dxB = dB1 * B2 * B3;
                    double dyB = B1 * dB2 * B3;
                    double dzB = B1 * B2 * dB3;

                    double c = U_prev(j1, j2, j3);
                    u += c * B;
                    dxu += c * dxB;
                    dyu += c * dyB;
                    dzu += c * dzB;
                }
                }
                }

                double gradient_prod = dxu * dxv + dyu * dyv + dzu * dzv;
                double val = u * v - dt * gradient_prod;
                U(i1, i2, i3) += val * w * J;
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
}
}



#endif /* ADS_PROBLEMS_HEAT_HPP_ */
