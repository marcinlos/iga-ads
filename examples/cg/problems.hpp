// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CG_PROBLEMS_HPP
#define CG_PROBLEMS_HPP

#include <array>

#include "../erikkson/solution.hpp"
#include "ads/util/function_value/function_value_2d.hpp"


using point = std::array<double, 2>;
using value = ads::function_value_2d;

struct erikkson {

    double epsilon;

    double diffusion(point /*x*/) const {
        return epsilon;
    }

    point beta(point /*x*/) const {
        return {1, 0};
    }

    double forcing(point /*x*/) const {
        return 0;
    }

    double boundary(point x) const {
        return x[0] == 0 ? std::sin(M_PI * x[1]) : 0;
    }

    value solution(point x) const {
        return ads::erikkson_exact(x[0], x[1], epsilon);
    }
};

struct problem1 {

    double epsilon;

    double diffusion(point /*x*/) const {
        return epsilon;
    }

    point beta(point /*x*/) const {
        return {1, 1};
    }

    double forcing(point x) const {
        return ads::erikkson2_forcing(x[0], x[1], epsilon);
    }

    double boundary(point /*x*/) const {
        return 0;
    }

    value solution(point x) const {
        return ads::erikkson2_exact(x[0], x[1], epsilon);
    }
};


#endif // CG_PROBLEMS_HPP
