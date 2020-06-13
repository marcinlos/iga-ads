#include <cmath>
#include <array>

#include "ads/util/function_value.hpp"


namespace ads {

    using value_type = function_value_3d;
    using point_type = std::array<double, 3>;
    using value_vec = std::array<value_type, 3>;


    struct maxwell_manufactured1 {
        double k;
        double l;
        double d;

        double eps = 1;
        double mu = 1;

        maxwell_manufactured1(double k, double l)
        : k{k}, l{l}, d{std::hypot(k, l)}
        { }

        value_type E1(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = cos(d * M_PI * t);
            auto val = sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);

            auto dx1 = 0.0;
            auto dx2 = k * M_PI * cos(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
            auto dx3 = l * M_PI * sin(k * M_PI * x[1]) * cos(l * M_PI * x[2]);

            return value_type{val, dx1, dx2, dx3} * time;
        }

        value_type E2(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = cos(d * M_PI * t);
            auto val = sin(k * M_PI * x[0]) * sin(l * M_PI * x[2]);

            auto dx1 = k * M_PI * cos(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
            auto dx2 = 0.0;
            auto dx3 = l * M_PI * sin(k * M_PI * x[0]) * cos(l * M_PI * x[2]);

            return 2 * value_type{val, dx1, dx2, dx3} * time;
        }

        value_type E3(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = cos(d * M_PI * t);
            auto val = sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);

            auto dx1 = k * M_PI * cos(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
            auto dx2 = l * M_PI * sin(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
            auto dx3 = 0.0;

            return 3 * value_type{val, dx1, dx2, dx3} * time;
        }

        value_type H1(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = sin(d * M_PI * t);

            auto val2 = sin(k * M_PI * x[0]) * cos(l * M_PI * x[2]);
            auto dv2x1 =   k * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[2]);
            auto dv2x2 =   0.0;
            auto dv2x3 = - l * M_PI * sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
            auto v2 = value_type{val2, dv2x1, dv2x2, dv2x3};

            auto val3 = sin(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
            auto dv3x1 =   k * M_PI * cos(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
            auto dv3x2 = - l * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
            auto dv3x3 =   0.0;
            auto v3 = value_type{val3, dv3x1, dv3x2, dv3x3};

            return (2 * (-l/d) * v2 + 3 * (-l/d) * v3) * time;
        }

        value_type H2(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = sin(d * M_PI * t);

            auto val1 = sin(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
            auto dv1x1 =   0.0;
            auto dv1x2 =   k * M_PI * cos(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
            auto dv1x3 = - l * M_PI * sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
            auto v1 = value_type{val1, dv1x1, dv1x2, dv1x3};

            auto val3 = cos(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
            auto dv3x1 = - k * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[1]);
            auto dv3x2 =   l * M_PI * sin(k * M_PI * x[0]) * cos(l * M_PI * x[1]);
            auto dv3x3 = 0.0;
            auto v3 = value_type{val3, dv3x1, dv3x2, dv3x3};

            return ((-l/d) * v1 + 3 * (k/d) * v3) * time;
        }

        value_type H3(point_type x, double t) const {
            using std::sin;
            using std::cos;

            auto time = sin(d * M_PI * t);

            auto val1 = cos(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
            auto dv1x1 =   0.0;
            auto dv1x2 = - k * M_PI * sin(k * M_PI * x[1]) * sin(l * M_PI * x[2]);
            auto dv1x3 =   l * M_PI * cos(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
            auto v1 = value_type{val1, dv1x1, dv1x2, dv1x3};

            auto val2 = cos(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
            auto dv2x1 = - k * M_PI * sin(k * M_PI * x[0]) * sin(l * M_PI * x[2]);
            auto dv2x2 =   0.0;
            auto dv2x3 =   l * M_PI * cos(k * M_PI * x[1]) * cos(l * M_PI * x[2]);
            auto v2 = value_type{val2, dv2x1, dv2x2, dv2x3};

            return ((k/d) * v1 + 2 * (k/d) * v2) * time;
        }

        auto E1_at(double t) const {
            return [this,t](point_type x) { return E1(x, t); };
        }

        auto E2_at(double t) const {
            return [this,t](point_type x) { return E2(x, t); };
        }

        auto E3_at(double t) const {
            return [this,t](point_type x) { return E3(x, t); };
        }

        auto H1_at(double t) const {
            return [this,t](point_type x) { return H1(x, t); };
        }

        auto H2_at(double t) const {
            return [this,t](point_type x) { return H2(x, t); };
        }

        auto H3_at(double t) const {
            return [this,t](point_type x) { return H3(x, t); };
        }

        auto E1_val_at(double t) const {
            return [this,t](point_type x) { return E1(x, t).val; };
        }

        auto E2_val_at(double t) const {
            return [this,t](point_type x) { return E2(x, t).val; };
        }

        auto E3_val_at(double t) const {
            return [this,t](point_type x) { return E3(x, t).val; };
        }

        auto H1_val_at(double t) const {
            return [this,t](point_type x) { return H1(x, t).val; };
        }

        auto H2_val_at(double t) const {
            return [this,t](point_type x) { return H2(x, t).val; };
        }

        auto H3_val_at(double t) const {
            return [this,t](point_type x) { return H3(x, t).val; };
        }
    };
}
