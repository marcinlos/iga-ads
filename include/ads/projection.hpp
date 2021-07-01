// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_PROJECTION_HPP_
#define ADS_PROJECTION_HPP_

#include "ads/basis_data.hpp"


namespace ads {

template <std::size_t N>
struct projector {
    std::array<basis_data, N> bases;

    using index_type = std::array<int, N>;

    using point_type = std::array<double, N>;

    double jacobian(index_type e) const {
        double j = 1;
        for (std::size_t i = 0u; i < N; ++ i) {
            j *= bases[i].J[e[i]];
        }
        return j;
    }

    double weight(index_type q) const {
        double w = 1;
        for (std::size_t i = 0u; i < N; ++ i) {
            w *= bases[i].w[q[i]];
        }
        return w;
    }

    point_type point(index_type e, index_type q) const {
        point_type p;
        for (std::size_t i = 0u; i < N; ++ i) {
            p[i] = bases[i].x[e[i]][q[i]];
        }
        return p;
    }

    template <typename Rhs, typename F>
    void operator () (Rhs& u, F&& f) {
        auto ev = evaluator<Rhs, F>{ *this, u, std::forward<F>(f) };
        ev.eval();
    }

    template <typename Rhs, typename F>
    struct evaluator {
        projector& proj;
        Rhs& u;
        F&& f;

        void eval() {
            index_type e;
            element_loop(e, 0);
        }

        void element_loop(index_type& e, std::size_t level) {
            if (level < N) {
                auto& basis = proj.bases[level];
                for (int i = 0; i < basis.elements; ++ i) {
                    e[level] = i;
                    element_loop(e, level + 1);
                }
            } else {
                double J = proj.jacobian(e);
                index_type q;
                quad_loop(e, q, J, 0);
            }
        }

        void quad_loop(const index_type& e, index_type& q, double J, std::size_t level) {
            if (level < N) {
                auto& basis = proj.bases[level];
                for (int i = 0; i < basis.quad_order; ++ i) {
                    q[level] = i;
                    quad_loop(e, q, J, level + 1);
                }
            } else {
                double w = proj.weight(q);
                auto x = proj.point(e, q);
                index_type a;
                dof_loop(e, q, a, J, w, x, 0);
            }
        }

        void dof_loop(const index_type& e, const index_type& q, index_type& a, double J, double w, const point_type& x, std::size_t level) {
            if (level < N) {
                auto& basis = proj.bases[level];
                for (int i = 0; i <= basis.degree; ++ i) {
                    a[level] = i;
                    dof_loop(e, q, a, J, w, x, level + 1);
                }
            } else {
                double B = 1;
                index_type idx = a;
                for (std::size_t i = 0u; i < N; ++ i) {
                    B *= proj.bases[i].b[e[i]][q[i]][0][a[i]];
                    idx[i] += proj.bases[i].first_dof(e[i]);
                }

                rhs(idx) += call_f(x) * B * w * J;
            }
        }

        double& rhs(index_type a) {
            return indexer(a);
        }

        double call_f(point_type x) {
            return call(x);
        }

        template <typename... Idx>
        std::enable_if_t<sizeof...(Idx) == N, double&>
        indexer(const index_type&, Idx... indices) {
            return u(indices...);
        }

        template <typename... Idx>
        std::enable_if_t<sizeof...(Idx) < N, double&>
        indexer(const index_type& idx, Idx... indices) {
            return indexer(idx, idx[N - sizeof...(Idx) - 1], indices...);
        }

        template <typename... Xs>
        std::enable_if_t<sizeof...(Xs) == N, double>
        call(const point_type&, Xs... xs) {
            return f(xs...);
        }

        template <typename... Xs>
        std::enable_if_t<sizeof...(Xs) < N, double>
        call(const point_type& x, Xs... xs) {
            return call(x, x[N - sizeof...(Xs) - 1], xs...);
        }
    };

};


template <typename Rhs, typename Function>
void compute_projection(Rhs& u, const basis_data& d1, Function f) {
    projector<1> project {{d1}};
    project(u, f);
}


template <typename Rhs, typename Function>
void compute_projection(Rhs& u, const basis_data& d1, const basis_data& d2, Function&& f) {
    projector<2> project {{d1, d2}};
    project(u, f);
}


template <typename Rhs, typename Function>
void compute_projection(Rhs& u, const basis_data& d1, const basis_data& d2, const basis_data& d3, Function&& f) {
    projector<3> project {{d1, d2, d3}};
    project(u, f);
}

}

#endif /* ADS_PROJECTION_HPP_ */
