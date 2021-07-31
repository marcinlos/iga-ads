// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include <exception>
#include <fstream>
#include <string>

#include "ads/executor/galois.hpp"
#include "ads/simulation.hpp"

namespace ads {

class diff_computer2d : public simulation_2d {
    galois_executor executor{4};

    vector_type read_solution(const std::string& path) {
        vector_type u{{{x.dofs(), y.dofs()}}};
        std::ifstream in{path};
        if (!in) {
            throw std::invalid_argument{path};
        }
        while (in) {
            int i;
            int j;
            double val;
            in >> i >> j >> val;
            u(i, j) = val;
        }
        return u;
    }

public:
    explicit diff_computer2d(const ads::config_2d& cfg)
    : simulation_2d(cfg) { }

    double error_L2(const std::string& path1, const std::string& path2) {
        auto u = read_solution(path1);
        auto v = read_solution(path2);

        double L2 = 0;

        executor.for_each(elements(), [&](index_type e) {
            double localL2 = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);

                auto uval = eval_fun(u, e, q);
                auto vval = eval_fun(v, e, q);
                auto e = uval - vval;
                localL2 += w * J * e.val * e.val;
            }
            executor.synchronized([&]() { L2 += localL2; });
        });
        return std::sqrt(L2);
    }

    double error_H1(const std::string& path1, const std::string& path2) {
        auto u = read_solution(path1);
        auto v = read_solution(path2);

        double H1 = 0;

        executor.for_each(elements(), [&](index_type e) {
            double localH1 = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);

                auto uval = eval_fun(u, e, q);
                auto vval = eval_fun(v, e, q);
                auto e = uval - vval;
                localH1 += w * J * (e.val * e.val + e.dx * e.dx + e.dy * e.dy);
            }
            executor.synchronized([&]() { H1 += localH1; });
        });
        return std::sqrt(H1);
    }
};

class diff_computer3d : public simulation_3d {
    galois_executor executor{4};

    struct state {
        vector_type ux, uy, uz;
    };

    state read_solution(const std::string& path) {
        std::array<int, 3> dims{{x.dofs(), y.dofs(), z.dofs()}};
        vector_type ux{dims};
        vector_type uy{dims};
        vector_type uz{dims};
        std::ifstream in{path};
        if (!in) {
            throw std::invalid_argument{path};
        }
        while (in) {
            int i;
            int j;
            int k;
            in >> i >> j >> k;
            in >> ux(i, j, k) >> uy(i, j, k) >> uz(i, j, k);
        }
        return {std::move(ux), std::move(uy), std::move(uz)};
    }

public:
    explicit diff_computer3d(const ads::config_3d& cfg)
    : simulation_3d(cfg) { }

    double error_L2(const std::string& path1, const std::string& path2) {
        auto a = read_solution(path1);
        auto b = read_solution(path2);

        double L2 = 0;

        executor.for_each(elements(), [&](index_type e) {
            double localL2 = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);

                auto ex = eval_fun(a.ux, e, q) - eval_fun(b.ux, e, q);
                auto ey = eval_fun(a.uy, e, q) - eval_fun(b.uy, e, q);
                auto ez = eval_fun(a.uz, e, q) - eval_fun(b.uz, e, q);

                localL2 += w * J * (ex.val * ex.val + ey.val * ey.val + ez.val * ez.val);
            }
            executor.synchronized([&]() { L2 += localL2; });
        });
        return std::sqrt(L2);
    }

    double error_H1(const std::string& path1, const std::string& path2) {
        auto a = read_solution(path1);
        auto b = read_solution(path2);

        double H1 = 0;

        executor.for_each(elements(), [&](index_type e) {
            double localH1 = 0;
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);

                auto ex = eval_fun(a.ux, e, q) - eval_fun(b.ux, e, q);
                auto ey = eval_fun(a.uy, e, q) - eval_fun(b.uy, e, q);
                auto ez = eval_fun(a.uz, e, q) - eval_fun(b.uz, e, q);
                auto norm2 = [](auto e) {
                    return e.val * e.val + e.dx * e.dx + e.dy * e.dy + e.dz * e.dz;
                };
                localH1 += w * J * (norm2(ex) + norm2(ey) + norm2(ez));
            }
            executor.synchronized([&]() { H1 += localH1; });
        });
        return std::sqrt(H1);
    }
};

}  // namespace ads

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: ./error <dim> <p> <n> <file1> <file2>" << std::endl;
        std::exit(-1);
    }
    auto d = std::atoi(argv[1]);
    auto p = std::atoi(argv[2]);
    auto n = std::atoi(argv[3]);
    auto* file1 = argv[4];
    auto* file2 = argv[5];

    ads::dim_config dim{p, n};
    int ders = 1;

    ads::timesteps_config steps{4000, 2.7e-2};

    if (d == 2) {
        ads::config_2d c{dim, dim, steps, ders};
        ads::diff_computer2d diff{c};
        std::cout << diff.error_L2(file1, file2) << ' ' << diff.error_H1(file1, file2) << std::endl;
    } else if (d == 3) {
        ads::config_3d c{dim, dim, dim, steps, ders};
        ads::diff_computer3d diff{c};
        std::cout << diff.error_L2(file1, file2) << ' ' << diff.error_H1(file1, file2) << std::endl;
    }
}
