#ifndef PROBLEMS_TUMOR_SKIN_HPP_
#define PROBLEMS_TUMOR_SKIN_HPP_

#include <cstddef>

namespace tumor {

class skin_model {
public:
    enum layer {
        stratum_corneum,
        stratum_spinosum,
        basement_membrame,
        dermis,
        hypodermis,
        count
    };

    static constexpr std::size_t layer_count = static_cast<std::size_t>(layer::count);

    double diffusion_coefficient[layer_count] = {
        20 * 1e-6, // stratum corneum
        83 * 1e-6, // stratum spinosum
        0.83 * 1e-6,  // basement membrame
        41.5 * 1e-6,  // dermis
        20 * 1e-6,  // hypodermis
    };

    double stratum_corneum_top = 3000;
    double stratum_spinosum_top = 2820;
    double basement_membrame_top = 2440;
    double dermis_top = 2400;
    double hypodermis_top = 600;

    double top[layer_count] = {
        stratum_corneum_top,
        stratum_spinosum_top,
        basement_membrame_top,
        dermis_top,
        hypodermis_top
    };

    layer layer_at(double /*x*/, double /*y*/, double z) const {
        for (std::size_t i = 1; i < layer_count; ++ i) {
            if (z > top[i]) {
                return static_cast<layer>(i - 1);
            }
        }
        return layer::hypodermis;
    }

    double diffusion(double x, double y, double z) const {
        return diffusion_coefficient[layer_at(x, y, z)];
    }

    double init_M(double x, double y, double z) const {
        return layer_at(x, y, z) == dermis ? 1.0 : 0.1;
    }
};


}

#endif /* PROBLEMS_TUMOR_SKIN_HPP_ */
