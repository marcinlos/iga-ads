// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_MT_SOURCE_HPP
#define MAXWELL_MT_SOURCE_HPP

#include <array>
#include <cmath>
#include <vector>

#include "ads/bspline/bspline.hpp"
#include "ads/quad/gauss.hpp"
#include "ads/util/function_value.hpp"
#include "spaces.hpp"
#include "state.hpp"

class mt_source {
private:

    using point_type = std::array<double, 3>;
    using index_type = std::array<int, 3>;
    using index_1d_iter_type = boost::counting_iterator<int>;
    using index_iter_type = ads::util::iter_product3<index_1d_iter_type, index_type>;
    using index_range = boost::iterator_range<index_iter_type>;
    using value_type = ads::function_value_3d;

public:
    auto apply_forcing(double /*t*/, state& /*rhs*/, space const& /*V*/) -> void { }
};

#endif  // MAXWELL_MT_SOURCE_HPP
