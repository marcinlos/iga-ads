#ifndef ADS_OUTPUT_MANAGER_HPP_
#define ADS_OUTPUT_MANAGER_HPP_

#include <cstddef>
#include <vector>
#include <fstream>
#include <boost/format.hpp>
#include "ads/bspline/bspline.hpp"
#include "ads/bspline/multidimensional.hpp"
#include "ads/output/output_format.hpp"
#include "ads/output/range.hpp"
#include "ads/output/grid.hpp"
#include "ads/output/vtk.hpp"


namespace ads {

struct output_manager {

    std::size_t n;
    const bspline::basis& bx;
    const bspline::basis& by;
    const bspline::basis& bz;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    lin::tensor<double, 3> vals;
    bspline::eval_ctx cx, cy, cz;

    output_manager(const bspline::basis& bx, const bspline::basis& by, const bspline::basis& bz, std::size_t n)
    : n { n }
    , bx { bx }
    , by { by }
    , bz { bz }
    , xs(n + 1)
    , ys(n + 1)
    , zs(n + 1)
    , vals {{ n + 1, n + 1, n + 1 }}
    , cx { bx.degree }
    , cy { by.degree }
    , cz { bz.degree }
    {
        for (std::size_t i = 0; i <= n; ++ i) {
            xs[i] = ads::lerp(i, n, bx.begin(), bx.end());
            ys[i] = ads::lerp(i, n, by.begin(), by.end());
            zs[i] = ads::lerp(i, n, bz.begin(), bz.end());
        }
    }

    template <typename Solution>
    void to_file(const Solution& sol, const std::string& file_pattern, int iter) {
        for (std::size_t i = 0; i <= n; ++ i) {
        for (std::size_t j = 0; j <= n; ++ j) {
        for (std::size_t k = 0; k <= n; ++ k) {
            vals(i, j, k) = bspline::eval(xs[i], ys[j], zs[k], sol, bx, by, bz, cx, cy, cz);
        }
        }
        }

        auto xrange = output::from_container(xs);
        auto yrange = output::from_container(ys);
        auto zrange = output::from_container(zs);
        auto grid = output::make_grid(xrange, yrange, zrange);

        output::vtk output { output::fixed_format(10, 18) };
        auto name = boost::format(file_pattern) % iter;
        std::ofstream os { name.str() };
        output.print(os, grid, vals);
    }
};


}


#endif /* ADS_OUTPUT_MANAGER_HPP_ */
