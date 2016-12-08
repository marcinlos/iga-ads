#ifndef ADS_PROBLEMS_VALIDATION_HPP_
#define ADS_PROBLEMS_VALIDATION_HPP_


#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"
#include "ads/executor/galois.hpp"


namespace ads {
namespace problems {


class validation : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    output_manager<2> output;
    galois_executor executor{4};

public:
    validation(const config_2d& config)
    : Base{ config }
    , u{ shape() }
    , u_prev{ shape() }
    , output{ x.B, y.B, 200 }
    { }

    double init_state(double x, double y) {
        return fi(x, y) * sc(0);
    };

private:
    void solve(vector_type& v) {
        for (int i = 0; i < y.dofs(); ++ i) {
            v(0, i) = 0;
            v(x.dofs() - 1, i) = 0;
        }
        for (int i = 0; i < x.dofs(); ++ i) {
            v(i, 0) = 0;
            v(i, y.dofs() - 1) = 0;
        }
        Base::solve(v);
    }

    static constexpr double k = 2 * M_PI * M_PI;

    double solution(double x, double y, double t) const {
        return sc(t) * fi(x, y);
    }

    double fi(double x, double y) const {
        return std::sin(x * M_PI) * std::sin(y * M_PI);
    }

    double sc(double t) const {
        return std::exp(-k*t);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();
        Base::prepare_matrices();
    }

    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double /*t*/) override {
        compute_rhs();
        solve(u);
    }

    void after_step(int iter, double t) override {
        if (iter % 1000 == 0) {
            output.to_file(u, "out_%d.data", iter);
            validate(t);
        }
    }

    void compute_rhs() {
        auto& rhs = u;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weigth(q);
                value_type u = eval_fun(u_prev, e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double gradient_prod = grad_dot(u, v);
                    double val = u.val * v.val - steps.dt * gradient_prod ;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() {
                update_global_rhs(rhs, U, e);
            });
        });
    }

    void validate(double t) {
        double L2 = 0, Max = 0;

        executor.for_each(elements(), [&](index_type e) {
            double localL2 = 0, localMax = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                auto x = point(e, q);
                double w = weigth(q);

                auto uval = eval_fun(u, e, q).val;
                auto sol = solution(x[0], x[1], t);
                auto e = std::abs(uval - sol);

                localMax = std::max(localMax, e);
                localL2 += w * J * e * e;
            }
            executor.synchronized([&]() {
                L2 += localL2;
                Max = std::max(Max, localMax);
            });
        });
        L2 = std::sqrt(L2);
        std::cout << t << "  " << L2 << "  " << Max << std::endl;
    }
};

}
}

#endif /* ADS_PROBLEMS_VALIDATION_HPP_ */
