#include "ads/simulation/dimension.hpp"

namespace ads {

    dimension::dimension(bspline::basis b, int quad_order, int derivatives, int elem_division)
    : p{b.degree}
    , elements{b.elements() * elem_division}
    , a{b.begin()}
    , b{b.end()}
    , B{ std::move(b) }
    , M{p, p, B.dofs()}
    , basis(B, derivatives, quad_order, elem_division)
    , ctx{M}
    {
        gram_matrix_1d(M, basis);
    }

    dimension::dimension(const dim_config& config, int derivatives)
    : dimension{ bspline_basis(config), config.quad_order, derivatives }
    { }

    void dimension::fix_dof(int k) {
        int last = dofs() - 1;
        for (int i = clamp(k - p, 0, last); i <= clamp(k + p, 0, last); ++ i) {
            M(k, i) = 0;
        }
        M(k, k) = 1;
    }

}
