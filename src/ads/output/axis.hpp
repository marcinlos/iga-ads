#ifndef ADS_OUTPUT_AXIS_HPP_
#define ADS_OUTPUT_AXIS_HPP_

#include <vector>
#include "ads/util.hpp"
#include "ads/bspline/bspline.hpp"
#include "ads/output/range.hpp"

namespace ads {
namespace output {

inline std::vector<double> linspace(const bspline::basis& basis, std::size_t n) {
    std::vector<double> xs(n + 1);
    for (std::size_t i = 0; i <= n; ++ i) {
        xs[i] = lerp(i, n, basis.begin(), basis.end());
    }
    return xs;
}

struct axis {
    const bspline::basis& basis;
    bspline::eval_ctx ctx;
    std::vector<double> points;

    using range_type = decltype(output::from_container(points));

    axis(const bspline::basis& basis, std::size_t intervals)
    : basis{ basis }
    , ctx{ basis.degree }
    , points{ linspace(basis, intervals) }
    { }

    std::size_t size() const {
        return points.size();
    }

    std::size_t intervals() const {
        return size() - 1;
    }

    range_type range() const {
        return output::from_container(points);
    }

    double operator [](int i) const {
        return points[i];
    }
};

}
}


#endif /* ADS_OUTPUT_AXIS_HPP_ */
