#ifndef ADS_PROBLEMS_MULTISTEP_BASE_HPP_
#define ADS_PROBLEMS_MULTISTEP_BASE_HPP_


#include "ads/simulation.hpp"
#include "problems/multistep/scheme.hpp"


namespace ads {
namespace problems {


class multistep_base {
protected:
    using coeffs = std::vector<double>;

    int s;
    coeffs as, bs;
    int order;

    std::vector<int> fibo;

public:
    multistep_base(scheme scm, int order)
    : s{ scm.s }
    , as{ std::move(scm.as) }
    , bs{ std::move(scm.bs) }
    , order{ order }
    , fibo(order)
    {
        compute_fib(fibo, order);
    }

protected:

    // Compute M + eta K
    void fill_matrix(lin::band_matrix& M, const basis_data& d, double eta) {
        for (element_id e = 0; e < d.elements; ++ e) {
            for (int q = 0; q < d.quad_order; ++ q) {
                int first = d.first_dof(e);
                int last = d.last_dof(e);
                for (int a = 0; a + first <= last; ++ a) {
                    for (int b = 0; b + first <= last; ++ b) {
                        int ia = a + first;
                        int ib = b + first;
                        auto va = d.b[e][q][0][a];
                        auto vb = d.b[e][q][0][b];
                        auto da = d.b[e][q][1][a];
                        auto db = d.b[e][q][1][b];
                        auto val = va * vb + eta * da * db;
                        M(ia, ib) += val * d.w[q] * d.J[e];
                    }
                }
            }
        }
    }

    void fix_dof(int k, const dimension& dim, lin::band_matrix& K) {
        int last = dim.dofs() - 1;
        for (int i = clamp(k - dim.p, 0, last); i <= clamp(k + dim.p, 0, last); ++ i) {
            K(k, i) = 0;
        }
        K(k, k) = 1;
    }

    void compute_fib(std::vector<int>& f, int d) const {
        lin::tensor<int, 2> buf{{d, d}};

        buf(0, 0) = 1;
        for (int i = 1; i < d; ++ i) {
            buf(i, 0) = buf(i - 1, 0);
            for (int j = 1; j < i; ++ j) {
                buf(i, j) = - buf(i - 1, j - 1) + buf(i - 1, j);
            }
            buf(i, i) = - buf(i - 1, i - 1);
        }

        for (int i = 0; i < d; ++ i) {
            f[i] = buf(d - 1, i);
        }
    }

};

}
}


#endif /* ADS_PROBLEMS_MULTISTEP_BASE_HPP_ */
