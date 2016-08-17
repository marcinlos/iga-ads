#ifndef ADS_PROBLEMS_TUMOR_TUMOR_HPP_
#define ADS_PROBLEMS_TUMOR_TUMOR_HPP_


#include "problems/tumor/skin.hpp"
#include "problems/tumor/state.hpp"
#include "problems/tumor/vasculature.hpp"
#include "problems/tumor/params.hpp"

#include "ads/simulation.hpp"
#include "ads/output_manager.hpp"

#include "ads/executor/galois.hpp"


#include <cmath>
#include <boost/format.hpp>


namespace ads {
namespace tumor {

class tumor : public simulation_2d {
private:
    using Base = simulation_2d;

    state now, prev;
    params p;

    vasc::vasculature vasculature;

    output_manager<2> output;

    int save_every = 100;

    int vasc_update_every = 10;

    galois_executor executor{8};

    bspline::eval_ctx xctx;
    bspline::eval_ctx yctx;
    bspline::eval_ders_ctx xdctx;
    bspline::eval_ders_ctx ydctx;

public:
    tumor(const config_2d& config, const params& params, vasc::vasculature vasculature);

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
        double dx = x - 0.5;
        double dy = y - 0.5;
        double r2 = std::min(12 * (dx * dx + dy * dy), 1.0);
        return 0.8 * (r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1);
    }

    double init_fibronectin(double x, double y) {
        auto lay = p.skin.layer_at(x, y, x);
        return lay == skin_model::layer::dermis || lay == skin_model::layer::hypodermis ? 0.8 : 0;
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
}



#endif /* ADS_PROBLEMS_TUMOR_TUMOR_HPP_ */
