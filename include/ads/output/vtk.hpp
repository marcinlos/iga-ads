// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_OUTPUT_VTK_HPP_
#define ADS_OUTPUT_VTK_HPP_

#include <sstream>

#include <boost/format.hpp>

#include "ads/lin/tensor/base.hpp"
#include "ads/lin/tensor/for_each.hpp"
#include "ads/output/grid.hpp"
#include "ads/output/output_base.hpp"
#include "ads/util/io.hpp"


namespace ads::output {

namespace impl {

template <typename... Iters>
struct vtk_print_helper : output_base {

    explicit vtk_print_helper(const output_format& format)
    : output_base { format }
    { }

    template <typename Value, typename... Values>
    void print(std::ostream& os, const grid<Iters...>&, const Value& value, const Values&... values) {
        util::stream_state_saver guard(os);
        prepare_stream(os);
        for_each_multiindex([this,&os,&value,&values...](auto... is) {
            this->print_row(os, value(is...), values(is...)...);
        }, value);
    }
};

}

struct vtk : output_base {

    explicit vtk(const output_format& format)
    : output_base { format }
    { }

    template <typename... Iters>
    void print_header(std::ostream& os, const grid<Iters...>& grid, std::size_t value_count) {
        using boost::format;
        constexpr auto paraview_dim = 3;
        auto origin = repeat(0, paraview_dim);
        auto extent = make_extent(dims(grid));
        auto spacing = repeat(1, paraview_dim);

        os << "<?xml version=\"1.0\"?>" << std::endl;
        os << "<VTKFile type=\"ImageData\" version=\"0.1\">" << std::endl;
        os << format("  <ImageData WholeExtent=\"%s\" origin=\"%s\" spacing=\"%s\">")
                % extent % origin % spacing<< std::endl;
        os << format("    <Piece Extent=\"%s\">") % extent << std::endl;
        os <<        "      <PointData Scalars=\"Result\">" << std::endl;
        os << format("        <DataArray Name=\"Result\"  type=\"Float32\" "
                "format=\"ascii\" NumberOfComponents=\"%d\">") % value_count << std::endl;
    }

    void print_end(std::ostream& os) {
        os << "        </DataArray>" << std::endl;
        os << "      </PointData>" << std::endl;
        os << "    </Piece>" << std::endl;
        os << "  </ImageData>" << std::endl;
        os << "</VTKFile>" << std::endl;
    }

    template <typename... Iters, typename... Values>
    void print(std::ostream& os, const grid<Iters...>& grid, const Values&... values) {
        auto value_count = sizeof...(Values);
        impl::vtk_print_helper<Iters...> printer { output_base::format };

        print_header(os, grid, value_count);
        printer.print(os, grid, values...);
        print_end(os);
    }

private:

    template <typename T>
    std::string repeat(T value, std::size_t n) {
        return join(' ', std::vector<T>(n, value));
    }

    template <typename T, std::size_t N>
    std::string make_extent(const std::array<T, N>& extents) {
        constexpr std::size_t dim = 3; // paraview needs 3D
        std::vector<T> exts(2 * dim);
        for (std::size_t i = 0; i < N; ++ i) {
            exts[2 * i + 1] = extents[i] - 1;
        }
        return join(' ', exts);
    }

    template <typename Sep, typename Cont>
    std::string join(Sep sep, const Cont& cont) {
        using std::begin;
        using std::endl;
        return join(sep, begin(cont), end(cont));
    }

    template <typename Sep, typename Iter>
    std::string join(Sep sep, Iter begin, Iter end) {
        std::stringstream ss;
        while (true) {
            ss << *begin++;
            if (begin != end) {
                ss << sep;
            } else {
                break;
            }
        }
        return ss.str();
    }

};

}

#endif /* ADS_OUTPUT_VTK_HPP_ */
