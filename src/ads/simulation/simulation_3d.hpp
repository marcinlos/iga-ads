#ifndef ADS_SIMULATION_SIMULATION_3D_HPP_
#define ADS_SIMULATION_SIMULATION_3D_HPP_

#include <array>
#include <cstddef>
#include "ads/simulation/dimension.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/solver.hpp"
#include "ads/projection.hpp"
#include "ads/simulation/simulation_base.hpp"


namespace ads {


class simulation_3d : public simulation_base {
protected:
    dimension x, y, z;

    using vector_type = lin::tensor<double, 3>;
    vector_type buffer;

    void solve(vector_type& rhs) {
        ads_solve(rhs, buffer, x.data(), y.data(), z.data());
    }

    template <typename Function>
    void projection(vector_type& v, Function f) {
        compute_projection(v, x.basis, y.basis, z.basis, f);
    }

    std::array<std::size_t, 3> shape() const {
        return {x.dofs(), y.dofs(), z.dofs()};
    }

    void prepare_matrices() {
        x.factorize_matrix();
        y.factorize_matrix();
        z.factorize_matrix();
    }

public:
    simulation_3d(const config_3d& config)
    : simulation_base{config.steps}
    , x{config.x, config.derivatives}
    , y{config.y, config.derivatives}
    , z{config.z, config.derivatives}
    , buffer{shape()}
    { }
};


}


#endif /* ADS_SIMULATION_SIMULATION_3D_HPP_ */
