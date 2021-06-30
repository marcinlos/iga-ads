#ifndef ADS_PROBLEMS_FLOW_ENVIRONMENT_HPP_
#define ADS_PROBLEMS_FLOW_ENVIRONMENT_HPP_

#include <random>

#include "flow.hpp"
#include "geometry.hpp"


namespace ads::problems {

using ads::vec3d;

struct path {
    std::vector<vec3d> points;

    double dist(const vec3d& p) const {
        double dist = std::numeric_limits<double>::infinity();
        for (auto i = 0u; i < points.size() - 1; ++ i) {
            dist = std::min(dist, dist_from_segment(p, points[i], points[i + 1]));
        }
        return dist;
    }
};


class environment {
    constexpr static double GROUND = 0.2;

    mutable std::mt19937 rng;

    std::vector<path> paths;

    const double MIN = 1.0;
    const double MAX = 1000.0;

    const std::size_t PATH_COUNT = 20;
    const std::size_t MIN_PATH_LEN = 10;
    const std::size_t MAX_PATH_LEN = 20;

    vec3d random_vector(double a, double b) const {
        std::uniform_real_distribution<> dist{a, b};
        return {dist(rng), dist(rng), dist(rng)};
    }

    std::size_t random_path_length() const {
        std::uniform_int_distribution<> dist{MIN_PATH_LEN, MAX_PATH_LEN};
        return dist(rng);
    }

public:

    environment(std::mt19937::result_type seed)
    : rng{ seed }
    {
        double step = 0.05;

        for (auto i = 0u; i < PATH_COUNT; ++ i) {
            auto length = random_path_length();
            std::vector<vec3d> path;
            auto p = random_vector(0.15, 0.85);
            auto dp = random_vector(-1, 1);
            path.push_back(p);

            for (auto j = 1u; j < length; ++ j) {
                auto ddp = random_vector(-1, 1);
                double cos = dot(dp, ddp) / (len(dp) * len(ddp));
                ddp = ddp - 0.2 * cos * dp;
                dp = dp + 0.4 * ddp;
                p = p + step * dp;
                path.push_back(p);
            }
            paths.push_back({path});
        }
    }

    double permeability(double x, double y, double z) const {
        if (z < GROUND) {
            return 0.2;
        } else {
            vec3d v{x, y, z};
            double dist = std::numeric_limits<double>::infinity();
            for (auto& path : paths) {
                dist = std::min(dist, path.dist(v));
            }
            return lerp(MIN, MAX, falloff(0, 0.06, dist));
        }
    }

    double init_state(double x, double y, double z) const {
        vec3d v{x, y, z};
        double dist = std::numeric_limits<double>::infinity();
        for (auto& path : paths) {
            dist = std::min(dist, path.dist(v));
        }
        double network = lerp(0.0, 1.0, falloff(0, 0.1, dist));
        return 0.1 * network * bump(0.3, 0.5, x, y, z);
    }


    auto permeability_fun() const {
        return [this](double x, double y, double z) {
            return permeability(x, y, z);
        };
    }

    struct helper_init_state {
        const environment* owner;
        double operator ()(double x, double y, double z) const {
            return owner->init_state(x, y, z);
        }
    };

    helper_init_state init_state_fun() const {
        return helper_init_state{this};
    }

};

}

#endif /* ADS_PROBLEMS_FLOW_ENVIRONMENT_HPP_ */
