// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_PROBLEMS_FLOW_PUMPS_HPP_
#define ADS_PROBLEMS_FLOW_PUMPS_HPP_

#include <vector>

#include "ads/util.hpp"
#include "geometry.hpp"


namespace ads::problems {

struct pumps {
    std::vector<ads::vec3d> sources;
    std::vector<ads::vec3d> sinks;

    static constexpr double radius = 0.15;
    static constexpr double pumping_strength = 1;
    static constexpr double draining_strength = 1e5;

    double pumping(double x, double y, double z) const {
        ads::vec3d v{x, y, z};
        double p = 0;
        for (const auto& pos : sources) {
            double dist = len(v - pos);
            p += pumping_strength * ads::falloff(0, radius, dist);
        }
        return p;
    }

    double draining(double u, double x, double y, double z) const {
        ads::vec3d v{x, y, z};
        double p = 0;
        for (const auto& pos : sinks) {
            double dist = len(v - pos);
            double s = draining_strength * ads::falloff(0, radius, dist);
            p += u * s;
        }
        return p;
    }

    auto pumping_fun() const {
        return [this](double x, double y, double z) {
            return pumping(x, y, z);
        };
    }
};

}

#endif /* ADS_PROBLEMS_FLOW_PUMPS_HPP_ */
