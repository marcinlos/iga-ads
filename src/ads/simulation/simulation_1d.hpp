#ifndef ADS_SIMULATION_SIMULATION_1D_HPP_
#define ADS_SIMULATION_SIMULATION_1D_HPP_

#include <array>
#include <cstddef>
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"


namespace ads {


class simulation_1d : public simulation_base {
protected:
    dimension x;

    using vector_type = lin::vector;

    void solve(vector_type& rhs) {
        ads_solve(rhs, x.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, f);
    }

    std::array<std::size_t, 1> shape() const {
        return {x.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
    }

public:
    simulation_1d(const config_1d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    { }
};


}


#endif /* ADS_SIMULATION_SIMULATION_1D_HPP_ */
