#include <fstream>

#include "problems/tumor/vasculature/plot.hpp"
#include "ads/util.hpp"
#include "ads/output_manager.hpp"
#include "ads/output/gnuplot.hpp"
#include "ads/output/range.hpp"
#include "ads/output/grid.hpp"

namespace ads {
namespace tumor {
namespace vasc {

    void plot(std::ostream& os, const val_array& v) {
        auto s = v.sizes();
        auto xs = linspace(0.0, 1.0, s[0] - 1);
        auto ys = linspace(0.0, 1.0, s[1] - 1);
        auto rx = output::from_container(xs);
        auto ry = output::from_container(xs);
        auto grid = output::make_grid(rx, ry);

        output::gnuplot_printer<2> printer { DEFAULT_FMT };
        printer.print(os, grid, v);
    }

    void plot(const std::string& file, const val_array& v) {
        std::ofstream os { file };
        plot(os, v);
    }

}
}
}
