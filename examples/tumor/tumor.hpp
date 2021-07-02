// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef TUMOR_TUMOR_HPP
#define TUMOR_TUMOR_HPP

#include <cmath>

#include <boost/format.hpp>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "params.hpp"
#include "skin.hpp"
#include "state.hpp"
#include "vasculature.hpp"


namespace bsp = ads::bspline;

namespace tumor {

class tumor_2d : public ads::simulation_2d {
private:
    static constexpr std::size_t Dim = 2;
    using Base = ads::simulation_2d;

    state<Dim> now, prev;
    params p;

    vasc::vasculature vasculature;

    ads::output_manager<2> output;

    int save_every;

    int vasc_update_every = 10;

    ads::galois_executor executor{4};

    ads::bspline::eval_ctx xctx;
    ads::bspline::eval_ctx yctx;
    ads::bspline::eval_ders_ctx xdctx;
    ads::bspline::eval_ders_ctx ydctx;

public:
    tumor_2d(const ads::config_2d& config, const params& params, int save_every, vasc::vasculature vasculature);

private:

    struct constant {
        double value;

        constant(double v)
        : value{ v }
        { }

        double operator ()(double, double) const {
            return value;
        }
    };

    double init_tumor(double x, double y) {
        double dx = (x - 1500) / 600;
        double dy = (y - 1500) / 600;
        double r2 = std::min(dx * dx + dy * dy, 1.0);
        return 0.8 * (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

    double init_M(double x, double y) {
        auto lay = p.skin.layer_at(x, y, x);
        return lay == skin_model::layer::dermis ? 1 : p.init_M;
    }

    void dirichlet() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();
    }

    void apply_boundary_conditions(vector_type& v) {
        for (int i = 0; i < y.dofs(); ++ i) {
            v(0, i) = 0;
            v(x.dofs() - 1, i) = 0;
        }
        for (int i = 0; i < x.dofs(); ++ i) {
            v(i, 0) = 0;
            v(i, y.dofs() - 1) = 0;
        }
    }

    point_type normalize(point_type p) const {
        return { (p[0] - x.a) / (x.b - x.a), (p[1] - y.a) / (y.b - y.a) };
    }

    void prepare_matrices() {
        dirichlet();
        Base::prepare_matrices();
    }

    void plot_vasculature(int iter) {
        using boost::format;
        vasculature.plot_veins(str(format("vasculature_%d.data") % iter));
        vasculature.plot_oxygen(str(format("oxygen_%d.data") % iter));
    }

    void before() override;

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(now, prev);

        now.clear();
    }

    void step(int iter, double /*t*/) override {
        compute_rhs();
        solve_all();
        update_vasculature(iter);
    }

    void save_to_file(int iter);

    void after_step(int iter, double /*t*/) override {
        int next_iter = iter + 1;
        if (next_iter % save_every == 0) {
            save_to_file(next_iter);
        }
    }

    void solve_all();

    value_type ensure_positive(value_type v) const {
        if (v.val < 0) {
            v.val = 0;
        }
        return v;
    }

    void compute_rhs();

    void update_vasculature(int iter);

    double& val(vector_type& v, index_type idx) {
        return v(idx[0], idx[1]);
    }
};

}

#endif // TUMOR_TUMOR_HPP
