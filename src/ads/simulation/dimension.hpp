#ifndef ADS_SIMULATION_DIMENSION_HPP_
#define ADS_SIMULATION_DIMENSION_HPP_

#include <boost/range/counting_range.hpp>

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/band_solve.hpp"
#include "ads/solver.hpp"
#include "ads/mass_matrix.hpp"
#include "ads/bspline/bspline.hpp"
#include "ads/basis_data.hpp"
#include "ads/util.hpp"
#include "ads/simulation/config.hpp"


namespace ads {

class dimension {
public:
    using element_range_type = decltype(boost::counting_range(0, 0));

    int p;
    int elements;
    double a;
    double b;
    bspline::basis B;
    lin::band_matrix M;
    basis_data basis;
    lin::solver_ctx ctx;

    dimension(const dim_config& config, int derivatives);

    int dofs() const {
        return B.dofs();
    }

    element_range_type element_indices() const {
        return boost::counting_range(0, elements);
    }

    dim_data data() {
        return {M, ctx};
    }

    void fix_dof(int k);

    void fix_left() {
        fix_dof(0);
    }

    void fix_right() {
        int last = dofs() - 1;
        fix_dof(last);
    }

    void factorize_matrix() {
        lin::factorize(M, ctx);
    }

private:

    static bspline::basis bspline_basis(const dim_config& config) {
        return bspline::create_basis(config.a, config.b, config.p, config.elements);
    }
};


}


#endif /* ADS_SIMULATION_DIMENSION_HPP_ */
