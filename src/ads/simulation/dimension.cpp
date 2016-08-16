#include "ads/simulation/dimension.hpp"

namespace ads {

    dimension::dimension(const dim_config& config, int derivatives)
    : p{config.p}
    , elements{config.elements}
    , a{config.a}
    , b{config.b}
    , B(bspline_basis(config))
    , M{p, p, B.dofs()}
    , basis(B, derivatives)
    , ctx{M}
    {
        gram_matrix_1d(M, basis);
    }

    void dimension::fix_dof(int k) {
        int last = dofs() - 1;
        for (int i = clamp(k - p, 0, last); i <= clamp(k + p, 0, last); ++ i) {
            M(k, i) = 0;
        }
        M(k, k) = 1;
    }

}
