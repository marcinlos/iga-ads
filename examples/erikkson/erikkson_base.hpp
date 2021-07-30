// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ERIKKSON_ERIKKSON_BASE_HPP
#define ERIKKSON_ERIKKSON_BASE_HPP

#include <fstream>

#include "ads/bspline/eval.hpp"
#include "ads/executor/galois.hpp"
#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "solution.hpp"

namespace ads {

class erikkson_base : public simulation_2d {
protected:
    using Base = simulation_2d;
    using Base::Base;

    template <typename RHS>
    void stationary_bc(RHS& u, dimension& Ux, dimension& Uy) {
        dirichlet_bc(u, boundary::left, Ux, Uy, [](double t) { return std::sin(M_PI * t); });
        dirichlet_bc(u, boundary::right, Ux, Uy, 0);
        dirichlet_bc(u, boundary::bottom, Ux, Uy, 0);
        dirichlet_bc(u, boundary::top, Ux, Uy, 0);
    }

    template <typename RHS>
    void skew_bc(RHS& u, dimension& Ux, dimension& Uy) {
        dirichlet_bc(u, boundary::left, Ux, Uy, [](double t) { return t < 0.5 ? 1 : 0; });
        dirichlet_bc(u, boundary::right, Ux, Uy, 0);
        dirichlet_bc(u, boundary::bottom, Ux, Uy, 1);
        dirichlet_bc(u, boundary::top, Ux, Uy, 0);
    }

    template <typename RHS>
    void zero_bc(RHS& u, dimension& Ux, dimension& Uy) {
        for_boundary_dofs(Ux, Uy, [&](index_type i) { u(i[0], i[1]) = 0; });
    }

    void plot_horizontal(const char* filename, double y0, const vector_type& u, const dimension& Ux,
                         const dimension& Uy) const {
        std::ofstream out{filename};
        auto ctx_x = bspline::eval_ctx{Ux.B.degree};
        auto ctx_y = bspline::eval_ctx{Uy.B.degree};

        auto print = [&](double xx) {
            auto val = bspline::eval(xx, y0, u, Ux.B, Uy.B, ctx_x, ctx_y);
            out << std::setprecision(16) << xx << " " << val << std::endl;
        };

        print(Ux.a);
        auto N = Ux.basis.quad_order;
        for (auto e : Ux.element_indices()) {
            std::vector<double> qs(Ux.basis.x[e], Ux.basis.x[e] + N);
            std::sort(begin(qs), end(qs));
            for (auto xx : qs) {
                print(xx);
            }
        }
        print(Ux.b);
    }

    void plot_vertical(const char* filename, double x0, const vector_type& u, const dimension& Ux,
                       const dimension& Uy) const {
        std::ofstream out{filename};
        auto ctx_x = bspline::eval_ctx{Ux.B.degree};
        auto ctx_y = bspline::eval_ctx{Uy.B.degree};

        auto print = [&](double yy) {
            auto val = bspline::eval(x0, yy, u, Ux.B, Uy.B, ctx_x, ctx_y);
            out << std::setprecision(16) << yy << " " << val << std::endl;
        };

        print(Uy.a);
        auto N = Uy.basis.quad_order;
        for (auto e : Uy.element_indices()) {
            std::vector<double> qs(Uy.basis.x[e], Uy.basis.x[e] + N);
            std::sort(begin(qs), end(qs));
            for (auto yy : qs) {
                print(yy);
            }
        }
        print(Uy.b);
    }

    void plot_middle(const char* filename, const vector_type& u, const dimension& Ux,
                     const dimension& Uy) const {
        plot_horizontal(filename, 0.5, u, Ux, Uy);
    }

    void print_solution(const char* filename, const vector_type& u, const dimension& Ux,
                        const dimension& Uy) const {
        std::ofstream out{filename};
        for (auto dof : dofs(Ux, Uy)) {
            out << dof[0] << " " << dof[1] << " " << u(dof[0], dof[1]) << std::endl;
        }
    }

    auto exact(double eps) const {
        return [eps](point_type x) { return erikkson_exact(x[0], x[1], eps); };
    }

    auto exact_nonstationary(double t) const {
        return [t](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };
    }
};

}  // namespace ads

#endif  // ERIKKSON_ERIKKSON_BASE_HPP
