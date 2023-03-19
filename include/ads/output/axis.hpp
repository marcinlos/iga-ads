// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_AXIS_HPP
#define ADS_OUTPUT_AXIS_HPP

#include <vector>

#include "ads/bspline/bspline.hpp"
#include "ads/output/range.hpp"
#include "ads/util.hpp"

namespace ads::output {

struct axis {
    const bspline::basis& basis;
    bspline::eval_ctx ctx;
    std::vector<double> points;

    using range_type = decltype(output::from_container(points));

    axis(const bspline::basis& basis, std::size_t intervals)
    : basis{basis}
    , ctx{basis.degree}
    , points{linspace(basis.begin(), basis.end(), intervals)} { }

    int size() const { return narrow_cast<int>(points.size()); }

    int intervals() const { return size() - 1; }

    range_type range() const { return output::from_container(points); }

    double operator[](int i) const { return points[i]; }
};

}  // namespace ads::output

#endif  // ADS_OUTPUT_AXIS_HPP
