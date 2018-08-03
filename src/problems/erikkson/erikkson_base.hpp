#ifndef PROBLEMS_ERIKKSON_ERIKKSON_BASE_HPP_
#define PROBLEMS_ERIKKSON_ERIKKSON_BASE_HPP_

#include <fstream>
#include "ads/simulation.hpp"
#include "ads/simulation/utils.hpp"
#include "problems/erikkson/solution.hpp"
#include "ads/executor/galois.hpp"
#include "ads/bspline/eval.hpp"


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

    void plot_middle(const char* filename, const vector_type& u, const dimension& Ux, const dimension& Uy) const {
        std::ofstream out{filename};
        bspline::eval_ctx ctx_x{ Ux.B.degree }, ctx_y{ Uy.B.degree };

        auto print = [&](double xx) {
            auto val = bspline::eval(xx, 0.5, u, Ux.B, Uy.B, ctx_x, ctx_y);
            out << std::setprecision(16) << xx << " " << val << std::endl;
        };

        print(0);
        auto N = Ux.basis.quad_order;
        for (auto e : Ux.element_indices()) {
            std::vector<double> qs(Ux.basis.x[e], Ux.basis.x[e] + N);
            std::sort(begin(qs), end(qs));
            for (auto xx : qs) {
                print(xx);
            }
        }
        print(1);
    }

    void print_solution(const char* filename, const vector_type& u, const dimension& Ux, const dimension& Uy) const {
        std::ofstream out{filename};
        for (auto dof : dofs(Ux, Uy)) {
            out << dof[0] << " " << dof[1] << " " << u(dof[0], dof[1]) << std::endl;
        }
    }

    auto exact(double eps) const {
        return [&](point_type x) { return erikkson_exact(x[0], x[1], eps); };
    }

    auto exact_nonstationary(double t) const {
        return [&](point_type x) { return erikkson_nonstationary_exact(x[0], x[1], t); };
    }
};

}

#endif // PROBLEMS_ERIKKSON_ERIKKSON_BASE_HPP_
