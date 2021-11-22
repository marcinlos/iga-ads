// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef MAXWELL_PROBLEMS_HPP
#define MAXWELL_PROBLEMS_HPP

#include <array>
#include <cmath>

#include "ads/util/function_value.hpp"
#include "ads/util/math/vec.hpp"

using value_type = ads::function_value_3d;
using point_type = std::array<double, 3>;

inline auto near(double a, double b) -> bool {
    return std::abs(a - b) < 1e-8;
}

inline auto normal(point_type x) -> ads::math::vec<3> {
    if (near(x[0], 0))
        return {-1, 0, 0};
    else if (near(x[0], 1))
        return {1, 0, 0};
    else if (near(x[1], 0))
        return {0, -1, 0};
    else if (near(x[1], 1))
        return {0, 1, 0};
    else if (near(x[2], 0))
        return {0, 0, -1};
    else if (near(x[2], 1))
        return {0, 0, 1};
    else
        return {0, 0, 0};
}

inline auto dot(ads::math::vec<3> a, point_type b) -> double {
    return a.x * b[0] + a.y * b[1] + a.z * b[2];
}

template <typename Self>
class maxwell_problem {
private:
    auto self() const -> Self const& { return static_cast<Self const&>(*this); };

public:
    auto E1_at(double t) const {
        return [this, t](point_type x) { return self().E1(x, t); };
    }

    auto E2_at(double t) const {
        return [this, t](point_type x) { return self().E2(x, t); };
    }

    auto E3_at(double t) const {
        return [this, t](point_type x) { return self().E3(x, t); };
    }

    auto H1_at(double t) const {
        return [this, t](point_type x) { return self().H1(x, t); };
    }

    auto H2_at(double t) const {
        return [this, t](point_type x) { return self().H2(x, t); };
    }

    auto H3_at(double t) const {
        return [this, t](point_type x) { return self().H3(x, t); };
    }

    auto E1_val_at(double t) const {
        return [this, t](point_type x) { return self().E1(x, t).val; };
    }

    auto E2_val_at(double t) const {
        return [this, t](point_type x) { return self().E2(x, t).val; };
    }

    auto E3_val_at(double t) const {
        return [this, t](point_type x) { return self().E3(x, t).val; };
    }

    auto H1_val_at(double t) const {
        return [this, t](point_type x) { return self().H1(x, t).val; };
    }

    auto H2_val_at(double t) const {
        return [this, t](point_type x) { return self().H2(x, t).val; };
    }

    auto H3_val_at(double t) const {
        return [this, t](point_type x) { return self().H3(x, t).val; };
    }

    auto init_E1() const { return E1_val_at(0); }
    auto init_E2() const { return E2_val_at(0); }
    auto init_E3() const { return E3_val_at(0); }

    auto init_H1() const { return H1_val_at(0); }
    auto init_H2() const { return H2_val_at(0); }
    auto init_H3() const { return H3_val_at(0); }
};

struct maxwell_manufactured1 : maxwell_problem<maxwell_manufactured1> {
    double k;
    double l;
    double d;

    maxwell_manufactured1(double k, double l)
    : k{k}
    , l{l}
    , d{std::hypot(k, l)} { }

    auto eps(point_type) const noexcept -> double { return 1.0; }

    auto mu(point_type) const noexcept -> double { return 1.0; }

    value_type E1(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = cos(d * M_PI * t);
        auto val = sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);

        auto dx1 = 0.0;
        auto dx2 = k * M_PI * cos(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
        auto dx3 = l * M_PI * sin(k * M_PI * x[1]) * cos(l * M_PI * x[2]);

        return value_type{val, dx1, dx2, dx3} * time;
    }

    value_type E2(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = cos(d * M_PI * t);
        auto val = sin(k * M_PI * x[0]) * sin(l * M_PI * x[2]);

        auto dx1 = k * M_PI * cos(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
        auto dx2 = 0.0;
        auto dx3 = l * M_PI * sin(k * M_PI * x[0]) * cos(l * M_PI * x[2]);

        return 2 * value_type{val, dx1, dx2, dx3} * time;
    }

    value_type E3(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = cos(d * M_PI * t);
        auto val = sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);

        auto dx1 = k * M_PI * cos(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
        auto dx2 = l * M_PI * sin(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
        auto dx3 = 0.0;

        return 3 * value_type{val, dx1, dx2, dx3} * time;
    }

    value_type H1(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = sin(d * M_PI * t);

        auto val2 = sin(k * M_PI * x[0]) * cos(l * M_PI * x[2]);
        auto dv2x1 = k * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[2]);
        auto dv2x2 = 0.0;
        auto dv2x3 = -l * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
        auto v2 = value_type{val2, dv2x1, dv2x2, dv2x3};

        auto val3 = sin(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
        auto dv3x1 = k * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
        auto dv3x2 = -l * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
        auto dv3x3 = 0.0;
        auto v3 = value_type{val3, dv3x1, dv3x2, dv3x3};

        return (2 * (l / d) * v2 + 3 * (-l / d) * v3) * time;
    }

    value_type H2(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = sin(d * M_PI * t);

        auto val1 = sin(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
        auto dv1x1 = 0.0;
        auto dv1x2 = k * M_PI * cos(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
        auto dv1x3 = -l * M_PI * sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
        auto v1 = value_type{val1, dv1x1, dv1x2, dv1x3};

        auto val3 = cos(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
        auto dv3x1 = -k * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
        auto dv3x2 = l * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
        auto dv3x3 = 0.0;
        auto v3 = value_type{val3, dv3x1, dv3x2, dv3x3};

        return ((-l / d) * v1 + 3 * (k / d) * v3) * time;
    }

    value_type H3(point_type x, double t) const {
        using std::cos;
        using std::sin;

        auto time = sin(d * M_PI * t);

        auto val1 = cos(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
        auto dv1x1 = 0.0;
        auto dv1x2 = -k * M_PI * sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
        auto dv1x3 = l * M_PI * cos(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
        auto v1 = value_type{val1, dv1x1, dv1x2, dv1x3};

        auto val2 = cos(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
        auto dv2x1 = -k * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
        auto dv2x2 = 0.0;
        auto dv2x3 = l * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[2]);
        auto v2 = value_type{val2, dv2x1, dv2x2, dv2x3};

        return ((k / d) * v1 + 2 * (-k / d) * v2) * time;
    }
};

#endif  // MAXWELL_PROBLEMS_HPP
