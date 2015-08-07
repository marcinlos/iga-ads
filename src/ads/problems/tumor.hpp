#ifndef ADS_PROBLEMS_TUMOR_HPP_
#define ADS_PROBLEMS_TUMOR_HPP_

#include "ads/basis_data.hpp"


namespace ads {
namespace problems {
namespace tumor {

struct params {

    double delta;
    double epsilon;
    double P0;
    double Gamma;
    double chi0;

    double P(double u) {
        return u >= 0 ? delta * P0 * u : 0;
    }

    double f(double u) {
        double v = 1 - u;
        return Gamma * u * u * v * v;
    }

    double f_der(double u) {
        double v = 1 - u;
        return 2 * Gamma * u * v * (v + u);
    }

    double chi(double u, double n) {
        return chi0 * u * n;
    }

    double Du_chi(double u, double n) {
        return chi0 * n;
    }

    double Dn_chi(double u, double n) {
        return chi0 * u;
    }

    double gamma_u(double u, double mi_u, double mi_n) {
        return P(u) * (mi_n - mi_u);
    }

    double gamma_n(double u, double mi_u, double mi_n) {
        return - gamma_u(u, mi_u, mi_n);
    }

    double energy(double u, double n, double dxu, double dyu) {
        double grad2 = dxu * dxu + dyu * dyu;
        return f(u) + epsilon * epsilon / 2 * grad2 + chi(u, n) + n * n / (2 * delta);
    }
};


template <typename Rhs>
void rhs2d(
    const Rhs& U_prev,
    const Rhs& N_prev,
    const Rhs& mi_U_prev,
    const Rhs& mi_N_prev,
    Rhs& U,
    Rhs& N,
    Rhs& mi_U,
    Rhs& mi_N,
    const basis_data& d1,
    const basis_data& d2,
    params params,
    double dt)
{
    zero(U);
    zero(N);
    zero(mi_U);
    zero(mi_N);

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

                double Bx = d1.b[e1][q1][0][a1];
                double By = d2.b[e2][q2][0][a2];
                double dBx = d1.b[e1][q1][1][a1];
                double dBy = d2.b[e2][q2][1][a2];

                double v = Bx * By;
                double dxv = dBx * By;
                double dyv = Bx * dBy;

                double u = 0, n = 0, mi_u = 0, mi_n = 0;
                double dxu = 0, dyu = 0;
                double dxn = 0, dyn = 0;
                double dx_mi_u = 0, dy_mi_u = 0;
                double dx_mi_n = 0, dy_mi_n = 0;

                for (int b1 = 0; b1 + first1 <= last1; ++ b1) {
                for (int b2 = 0; b2 + first2 <= last2; ++ b2) {
                    int j1 = b1 + first1;
                    int j2 = b2 + first2;

                    double Bx  = d1.b[e1][q1][0][b1];
                    double By  = d2.b[e2][q2][0][b2];
                    double dBx = d1.b[e1][q1][1][b1];
                    double dBy = d2.b[e2][q2][1][b2];

                    double B   = Bx * By;
                    double dxB = dBx * By;
                    double dyB = Bx * dBy;

                    double cu = U_prev(j1, j2);
                    u   += cu * B;
                    dxu += cu * dxB;
                    dyu += cu * dyB;

                    double cn = N_prev(j1, j2);
                    n   += cn * B;
                    dxn += cn * dxB;
                    dyn += cn * dyB;

                    double c_mi_u = mi_U_prev(j1, j2);
                    mi_u    += c_mi_u * B;
                    dx_mi_u += c_mi_u * dxB;
                    dy_mi_u += c_mi_u * dyB;

                    double c_mi_n = mi_N_prev(j1, j2);
                    mi_n    += c_mi_n * B;
                    dx_mi_n += c_mi_n * dxB;
                    dy_mi_n += c_mi_n * dyB;
                }
                }

                double grad_u = dxu * dxv + dyu * dyv;
                double grad_mi_u = dx_mi_u * dxv + dy_mi_u * dyv;
                double grad_mi_n = dx_mi_n * dxv + dy_mi_n * dyv;

                double gamma_u = params.gamma_u(u, mi_u, mi_n);
                double gamma_n = params.gamma_n(u, mi_u, mi_n);

                double Du_chi = params.Du_chi(u, n);
                double Dn_chi = params.Dn_chi(u, n);

                double e2 = params.epsilon * params.epsilon;
                double df = params.f_der(u);

                double u_delta = gamma_u - grad_mi_u;
                U(i1, i2) += (u * v + dt * u_delta) * w * J;

                double n_delta = gamma_n - grad_mi_n;
                N(i1, i2) += (n * v + dt * n_delta) * w * J;


                double mi_u_val = df + e2 * grad_u + Du_chi;
                mi_U(i1, i2) += mi_u_val * w * J;

                double mi_n_val = Dn_chi + n / params.delta;
                mi_N(i1, i2) += mi_n_val * w * J;
            }
            }
        }
        }
    }
    }
}


struct solution_info {
    double energy;
    double mass;
};



template <typename Rhs>
solution_info info2d(Rhs& U, Rhs& N, const basis_data& d1, const basis_data& d2, params params) {

    double E = 0;
    double mass = 0;

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

            double u = 0, n = 0;
            double dxu = 0, dyu = 0;
            for (int a1 = 0; a1 + first1 <= last1; ++ a1) {
            for (int a2 = 0; a2 + first2 <= last2; ++ a2) {
                int i1 = a1 + first1;
                int i2 = a2 + first2;

                double Bx  = d1.b[e1][q1][0][a1];
                double By  = d2.b[e2][q2][0][a2];
                double dBx = d1.b[e1][q1][1][a1];
                double dBy = d2.b[e2][q2][1][a2];

                double B = Bx * By;
                double dxB = dBx * By;
                double dyB = Bx * dBy;

                double c = U(i1, i2);
                u   += c * B;
                dxu += c * dxB;
                dyu += c * dyB;
                n += N(i1, i2) * B;
            }
            }
            E += params.energy(u, n, dxu, dyu) * w * J;
            mass += (u + n) * w * J;
        }
        }
    }
    }
    return { E, mass };
}


}
}
}


#endif /* ADS_PROBLEMS_TUMOR_HPP_ */
