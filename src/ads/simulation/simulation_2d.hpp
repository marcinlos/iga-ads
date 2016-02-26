#ifndef ADS_SIMULATION_SIMULATION_2D_HPP_
#define ADS_SIMULATION_SIMULATION_2D_HPP_

#include <array>
#include <cstddef>
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"


namespace ads {


class simulation_2d : public simulation_base {
protected:
    dimension x, y;

    using vector_type = lin::tensor<double, 2>;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, f);
    }

    std::array<std::size_t, 2> shape() const {
        return {x.dofs(), y.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
    }

public:
    simulation_2d(const config_2d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , buffer{shape()}
    { }
};


#endif /* ADS_SIMULATION_SIMULATION_2D_HPP_ */
