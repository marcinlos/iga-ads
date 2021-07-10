// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_MANAGER_HPP
#define ADS_OUTPUT_MANAGER_HPP

#include <cstddef>
#include <ostream>
#include <vector>

#include "ads/bspline/bspline.hpp"
#include "ads/bspline/eval.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/output/axis.hpp"
#include "ads/output/gnuplot.hpp"
#include "ads/output/grid.hpp"
#include "ads/output/output_format.hpp"
#include "ads/output/output_manager_base.hpp"
#include "ads/output/vtk.hpp"

namespace ads {

const output::output_format DEFAULT_FMT = output::fixed_format(10, 18);

template <std::size_t Dim>
struct output_manager;

template <>
struct output_manager<1> : output_manager_base<output_manager<1>> {
private:
    output::axis x;
    lin::tensor<double, 1> vals;
    output::gnuplot_printer<1> output{DEFAULT_FMT};

public:
    output_manager(const bspline::basis& bx, std::size_t n)
    : x{bx, n}
    , vals{{x.size()}} { }

    using output_manager_base::to_file;

    template <typename Solution>
    void write(const Solution& sol, std::ostream& os) {
        for (int i = 0; i < x.size(); ++i) {
            vals(i) = bspline::eval(x[i], sol, x.basis, x.ctx);
        }
        auto grid = make_grid(x.range());
        output.print(os, grid, vals);
    }
};

template <>
struct output_manager<2> : output_manager_base<output_manager<2>> {
private:
    output::axis x, y;
    lin::tensor<double, 2> vals;
    output::gnuplot_printer<2> output{DEFAULT_FMT};

public:
    output_manager(const bspline::basis& bx, const bspline::basis& by, std::size_t n)
    : x{bx, n}
    , y{by, n}
    , vals{{x.size(), y.size()}} { }

    using output_manager_base::to_file;

    template <typename Solution>
    void write(const Solution& sol, std::ostream& os) {
        for (int i = 0; i < x.size(); ++i) {
            for (int j = 0; j < y.size(); ++j) {
                vals(i, j) = bspline::eval(x[i], y[j], sol, x.basis, y.basis, x.ctx, y.ctx);
            }
        }
        auto grid = make_grid(x.range(), y.range());
        output.print(os, grid, vals);
    }
};

template <>
struct output_manager<3> : output_manager_base<output_manager<3>> {
private:
    using value_array = lin::tensor<double, 3>;
    output::axis x, y, z;
    value_array vals;
    output::vtk output{DEFAULT_FMT};

public:
    output_manager(const bspline::basis& bx, const bspline::basis& by, const bspline::basis& bz,
                   std::size_t n)
    : x{bx, n}
    , y{by, n}
    , z{bz, n}
    , vals{{x.size(), y.size(), z.size()}} { }

    using output_manager_base::to_file;

    template <typename Solution>
    void evaluate(const Solution& sol, value_array& out) {
        for (int i = 0; i < x.size(); ++i) {
            for (int j = 0; j < y.size(); ++j) {
                for (int k = 0; k < z.size(); ++k) {
                    out(i, j, k) = bspline::eval(x[i], y[j], z[k], sol, x.basis, y.basis, z.basis,
                                                 x.ctx, y.ctx, z.ctx);
                }
            }
        }
    }

    template <typename Solution>
    value_array evaluate(const Solution& sol) {
        value_array out{{x.size(), y.size(), z.size()}};
        evaluate(sol, out);
        return out;
    }

    template <typename Solution>
    void write(const Solution& sol, std::ostream& os) {
        evaluate(sol, vals);
        auto grid = make_grid(x.range(), y.range(), z.range());
        output.print(os, grid, vals);
    }

    template <typename... Values>
    void to_file(const std::string& file_pattern, int iter, const Values&... values) {
        auto name = str(boost::format(file_pattern) % iter);
        to_file(name, values...);
    }

    template <typename... Values>
    void to_file(const std::string& output_file, const Values&... values) {
        std::ofstream os{output_file};
        auto grid = make_grid(x.range(), y.range(), z.range());
        output.print(os, grid, values...);
    }

private:
    template <typename... Values>
    void print(std::ostream& os, const Values&... values) {
        auto grid = make_grid(x.range(), y.range(), z.range());
        output.print(os, grid, values...);
    }
};

}  // namespace ads

#endif  // ADS_OUTPUT_MANAGER_HPP
