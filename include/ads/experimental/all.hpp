// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_EXPERIMENTAL_ALL_HPP
#define ADS_EXPERIMENTAL_ALL_HPP

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <optional>
#include <tuple>
#include <type_traits>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/counting_range.hpp>
#include <fmt/chrono.h>
#include <fmt/os.h>

#include "ads/bspline/bspline.hpp"
#include "ads/bspline/eval.hpp"
#include "ads/executor/galois.hpp"
#include "ads/executor/sequential.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/quad/gauss.hpp"
#include "ads/solver/mumps.hpp"
#include "ads/util.hpp"
#include "ads/util/function_value.hpp"
#include "ads/util/iter/product.hpp"

namespace ads {

using partition = std::vector<double>;

using simple_index = int;
using simple_index_iterator = boost::counting_iterator<simple_index>;
using simple_index_range = boost::iterator_range<simple_index_iterator>;

auto range(int start, int past_end) noexcept -> simple_index_range {
    return boost::counting_range(start, past_end);
}

auto empty_range() noexcept -> simple_index_range {
    return range(0, 0);
}

struct index_types {
    using index = std::tuple<int, int>;
    using index_iterator = util::iter_product2<simple_index_iterator, index>;
    using index_range = boost::iterator_range<index_iterator>;
};

struct interval {
    double left;
    double right;

    constexpr interval(double left, double right) noexcept
    : left{left}
    , right{right} { }
};

constexpr auto length(interval s) noexcept -> double {
    return std::abs(s.right - s.left);
}

auto subinterval(const partition& points, int i) noexcept -> interval {
    assert(i >= 0 && i < as_signed(points.size() - 1) && "Subinterval index out of range");
    const auto a = points[i];
    const auto b = points[i + 1];
    return {a, b};
}

constexpr auto lerp(double t, interval s) noexcept -> double {
    return (1 - t) * s.left + t * s.right;
}

auto operator<<(std::ostream& os, interval s) -> std::ostream& {
    return os << "[" << s.left << ", " << s.right << "]";
}

enum class facet_type : bool {
    interior = true,
    boundary = false,
};

class interval_mesh {
private:
    partition points_;

public:
    using point = double;
    using element_index = simple_index;
    using element_iterator = simple_index_iterator;
    using element_range = simple_index_range;
    using facet_index = simple_index;
    using facet_range = simple_index_range;

    explicit interval_mesh(partition points) noexcept
    : points_{std::move(points)} { }

    auto elements() const noexcept -> element_range { return range(0, element_count()); }

    auto element_count() const noexcept -> int { return narrow_cast<int>(points_.size()) - 1; }

    auto subinterval(element_index e) const noexcept -> interval {
        return ads::subinterval(points_, e);
    }

    struct point_data {
        point position;
        facet_type type;
        double normal;
    };

    auto facets() const noexcept -> facet_range {
        const auto facet_count = narrow_cast<int>(points_.size());
        return range(0, facet_count);
    }

    auto boundary_facets() const noexcept -> std::array<facet_index, 2> {
        const auto last = narrow_cast<int>(points_.size()) - 1;
        return {0, last};
    }

    auto interior_facets() const noexcept -> facet_range {
        const auto facet_count = narrow_cast<int>(points_.size());
        return range(1, facet_count - 1);
    }

    auto facet(facet_index i) const noexcept -> point_data {
        assert(i < as_signed(points_.size()) && "Point index out of range");
        // all points are positive except for the first one
        const auto normal = i > 0 ? 1.0 : -1.0;
        const auto type = (i > 0 && i < as_signed(points_.size()) - 1) ? facet_type::interior
                                                                       : facet_type::boundary;
        return {points_[i], type, normal};
    }
};

enum class orientation {
    horizontal,
    vertical,
};

class regular_mesh {
private:
    interval_mesh mesh_x_;
    interval_mesh mesh_y_;

public:
    using point = std::tuple<double, double>;
    using element_index = index_types::index;
    using element_iterator = index_types::index_iterator;
    using element_range = index_types::index_range;

    struct edge_index {
        simple_index ix;
        simple_index iy;
        orientation dir;
    };

    using facet_index = edge_index;
    // using facet_iterator   = void;
    // using facet_range      = void;

    regular_mesh(partition xs, partition ys) noexcept
    : mesh_x_{std::move(xs)}
    , mesh_y_{std::move(ys)} { }

    auto elements() const noexcept -> element_range {
        const auto rx = mesh_x_.elements();
        const auto ry = mesh_y_.elements();
        return util::product_range<element_index>(rx, ry);
    }

    struct element_data {
        interval span_x;
        interval span_y;
    };

    struct edge_data {
        interval span;
        double position;
        orientation direction;
        facet_type type;
        point normal;
    };

    using facet_data = edge_data;

    auto element(element_index e) const noexcept -> element_data {
        const auto [ix, iy] = e;
        const auto sx = mesh_x_.subinterval(ix);
        const auto sy = mesh_y_.subinterval(iy);
        return {sx, sy};
    }

    auto facets() const noexcept -> std::vector<edge_index> {
        auto indices = std::vector<edge_index>{};

        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.facets()) {
                indices.push_back({ix, iy, orientation::horizontal});
            }
        }
        for (auto ix : mesh_x_.facets()) {
            for (auto iy : mesh_y_.elements()) {
                indices.push_back({ix, iy, orientation::vertical});
            }
        }
        return indices;
    }

    auto boundary_facets() const noexcept -> std::vector<edge_index> {
        auto indices = std::vector<edge_index>{};

        for (auto iy : mesh_y_.boundary_facets()) {
            for (auto ix : mesh_x_.elements()) {
                indices.push_back({ix, iy, orientation::horizontal});
            }
        }
        for (auto ix : mesh_x_.boundary_facets()) {
            for (auto iy : mesh_y_.elements()) {
                indices.push_back({ix, iy, orientation::vertical});
            }
        }
        return indices;
    }

    auto interior_facets() const noexcept -> std::vector<edge_index> {
        auto indices = std::vector<edge_index>{};

        for (auto iy : mesh_y_.interior_facets()) {
            for (auto ix : mesh_x_.elements()) {
                indices.push_back({ix, iy, orientation::horizontal});
            }
        }
        for (auto ix : mesh_x_.interior_facets()) {
            for (auto iy : mesh_y_.elements()) {
                indices.push_back({ix, iy, orientation::vertical});
            }
        }
        return indices;
    }

    auto facet(edge_index e) const noexcept -> edge_data {
        const auto [ix, iy, dir] = e;

        if (dir == orientation::horizontal) {
            const auto [y, type, ny] = mesh_y_.facet(iy);
            const auto span = mesh_x_.subinterval(ix);
            const auto normal = point{0, ny};
            return {span, y, dir, type, normal};
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            const auto [x, type, nx] = mesh_x_.facet(ix);
            const auto span = mesh_y_.subinterval(iy);
            const auto normal = point{nx, 0};
            return {span, x, dir, type, normal};
        }
    }
};

auto operator<<(std::ostream& os, const regular_mesh::edge_index& edge) -> std::ostream& {
    const auto [ix, iy, dir] = edge;
    const char* sign = dir == orientation::horizontal ? "-" : "|";
    return os << "(" << ix << ", " << iy << ")[" << sign << "]";
}

using local_dof = simple_index;
using global_dof = simple_index;

auto spans_for_elements(const bspline::basis& b) -> std::vector<int> {
    auto spans = std::vector<int>{};
    spans.reserve(b.elements());

    for (int i = 0; i + 1 < as_signed(b.knot_size()); ++i) {
        if (b.knot[i] != b.knot[i + 1]) {
            spans.push_back(i);
        }
    }
    assert(as_signed(spans.size()) == b.elements());
    return spans;
}

class bspline_space {
private:
    bspline::basis basis_;
    std::vector<int> first_dofs_;
    std::vector<int> spans_;

public:
    using point = double;
    using element_index = simple_index;
    using facet_index = simple_index;
    using dof_index = simple_index;
    using dof_iterator = simple_index_iterator;
    using dof_range = simple_index_range;

    explicit bspline_space(bspline::basis basis)
    : basis_{std::move(basis)}
    , first_dofs_{first_nonzero_dofs(basis_)}
    , spans_{spans_for_elements(basis_)} { }

    auto basis() const noexcept -> const bspline::basis& { return basis_; }

    auto degree() const noexcept -> int { return basis_.degree; }

    auto dofs_per_element() const noexcept -> int { return degree() + 1; }

    auto dof_count() const noexcept -> int { return basis_.dofs(); }

    auto dof_count(element_index) const noexcept -> int { return dofs_per_element(); }

    auto facet_dof_count(facet_index f) const noexcept -> int {
        auto dofs = dofs_on_facet(f);
        using std::begin;
        using std::end;
        return narrow_cast<int>(std::distance(begin(dofs), end(dofs)));
    }

    auto dofs() const noexcept -> dof_range { return range(0, dof_count()); }

    auto dofs(element_index e) const noexcept -> dof_range {
        const auto first = first_dofs_[e];
        return range(first, first + dofs_per_element());
    }

    auto local_index(dof_index dof, element_index e) const noexcept -> local_dof {
        const auto first = first_dofs_[e];
        return dof - first;
    }

    auto first_dof(element_index e) const noexcept -> global_dof { return first_dofs_[e]; }

    auto last_dof(element_index e) const noexcept -> global_dof {
        return first_dof(e) + dofs_per_element() - 1;
    }

    auto dofs_on_facet(facet_index f) const noexcept -> dof_range {
        const auto last_element = basis_.elements() - 1;
        const auto elem_left = std::max(f - 1, 0);
        const auto elem_right = std::min(f, last_element);
        const auto first = first_dofs_[elem_left];
        const auto one_past_last = first_dofs_[elem_right] + dofs_per_element();
        return range(first, one_past_last);
    }

    auto facet_local_index(dof_index dof, facet_index f) const noexcept -> local_dof {
        const auto elem_left = std::max(f - 1, 0);
        const auto first = first_dofs_[elem_left];
        return dof - first;
    }

    auto span(element_index e) const noexcept -> int { return spans_[e]; }
};

class bspline_basis_values {
private:
    std::vector<double> buffer_;
    std::vector<double**> point_data_;  // indexed by point
    std::vector<double*> dof_data_;     // indexed by point and derivative

public:
    bspline_basis_values(int points, int dofs, int ders)
    : buffer_(points * dofs * (ders + 1))
    , point_data_(points)
    , dof_data_(points * (ders + 1)) {
        for (int i = 0; i < points * (ders + 1); ++i) {
            dof_data_[i] = &buffer_[i * dofs];
        }
        for (int i = 0; i < points; ++i) {
            point_data_[i] = &dof_data_[i * (ders + 1)];
        }
    }

    auto operator()(int point, local_dof i, int der) const -> double {
        return point_data_[point][der][i];
    }

    auto point_buffer(int point) noexcept -> double** { return point_data_[point]; }
};

auto evaluate_basis(const std::vector<double>& points, const bspline_space& space, int ders)
    -> bspline_basis_values {
    const auto point_count = static_cast<int>(points.size());
    const auto dof_count = space.dofs_per_element();

    auto values = bspline_basis_values{point_count, dof_count, ders};
    auto context = bspline::eval_ctx{space.degree()};

    for (int q = 0; q < as_signed(points.size()); ++q) {
        const auto x = points[q];
        auto* const buffer = values.point_buffer(q);
        const auto span = find_span(x, space.basis());
        eval_basis_with_derivatives(span, x, space.basis(), buffer, ders, context);
    }
    return values;
}

class bspline_basis_values_on_vertex {
private:
    std::optional<bspline_basis_values> left_;
    std::optional<bspline_basis_values> right_;
    local_dof left_last_;
    local_dof right_first_;

public:
    bspline_basis_values_on_vertex(std::optional<bspline_basis_values> left,
                                   std::optional<bspline_basis_values> right, local_dof left_last,
                                   local_dof right_first) noexcept
    : left_{std::move(left)}
    , right_{std::move(right)}
    , left_last_{left_last}
    , right_first_{right_first} {
        assert((left_ || right_) && "Neither left nor right adjacent element data specified");
    }

    auto operator()(local_dof i, int der) const noexcept -> double {
        if (left_ && i <= left_last_) {
            return left(i, der);
        } else {  // right_ has value
            return right(i, der);
        }
    }

    auto left(local_dof i, int der) const noexcept -> double {
        if (left_ && i <= left_last_) {
            return (*left_)(0, i, der);
        } else {
            return 0.0;
        }
    }

    auto right(local_dof i, int der) const noexcept -> double {
        if (right_ && i >= right_first_) {
            const auto idx = i - right_first_;
            return (*right_)(0, idx, der);
        } else {
            return 0.0;
        }
    }

    auto jump(local_dof i, int der, double normal) const noexcept -> double {
        const auto left_val = left(i, der);
        const auto right_val = right(i, der);
        return normal * (left_val - right_val);
    }

    auto average(local_dof i, int der) const noexcept -> double {
        const auto left_val = left(i, der);
        const auto right_val = right(i, der);
        const auto sum = left_val + right_val;
        if (left_ && right_) {
            return sum / 2;
        } else {
            // one of these is zero
            return sum;
        }
    }
};

auto evaluate_basis_at_point(double x, const bspline_space& space, int ders, int span)
    -> bspline_basis_values {
    const auto dof_count = space.dofs_per_element();

    auto values = bspline_basis_values{1, dof_count, ders};
    auto context = bspline::eval_ctx{space.degree()};

    auto* const buffer = values.point_buffer(0);

    eval_basis_with_derivatives(span, x, space.basis(), buffer, ders, context);

    return values;
}

auto element_left(bspline_space::facet_index f, const bspline::basis&) -> std::optional<int> {
    return f > 0 ? f - 1 : std::optional<int>{};
}

auto element_right(bspline_space::facet_index f, const bspline::basis& b) -> std::optional<int> {
    return f < b.elements() ? f : std::optional<int>{};
}

auto evaluate_basis(bspline_space::facet_index f, const bspline_space& space, int ders)
    -> bspline_basis_values_on_vertex {
    const auto& basis = space.basis();
    const auto x = basis.points[f];

    const auto maybe_elem_left = element_left(f, space.basis());
    const auto maybe_elem_right = element_right(f, space.basis());

    if (maybe_elem_left && maybe_elem_right) {
        const auto elem_left = maybe_elem_left.value();
        const auto elem_right = maybe_elem_right.value();
        const auto span_left = space.span(elem_left);
        const auto span_right = space.span(elem_right);
        const auto left_last = space.last_dof(elem_left);
        const auto right_first = space.first_dof(elem_right);
        const auto left_last_loc = space.facet_local_index(left_last, f);
        const auto right_first_loc = space.facet_local_index(right_first, f);

        auto vals_left = evaluate_basis_at_point(x, space, ders, span_left);
        auto vals_right = evaluate_basis_at_point(x, space, ders, span_right);

        return {std::move(vals_left), std::move(vals_right), left_last_loc, right_first_loc};

    } else if (maybe_elem_left) {
        const auto elem_left = maybe_elem_left.value();
        const auto span_left = space.span(elem_left);
        const auto left_last = space.last_dof(elem_left);
        const auto left_last_loc = space.facet_local_index(left_last, f);

        auto vals_left = evaluate_basis_at_point(x, space, ders, span_left);

        return {std::move(vals_left), {}, left_last_loc, {}};

    } else {  // maybe_elem_right
        assert(maybe_elem_right && "No elements adjacent to specified face");
        const auto elem_right = maybe_elem_right.value();
        const auto span_right = space.span(elem_right);
        const auto right_first = space.first_dof(elem_right);
        const auto right_first_loc = space.facet_local_index(right_first, f);

        auto vals_right = evaluate_basis_at_point(x, space, ders, span_right);

        return {{}, std::move(vals_right), {}, right_first_loc};
    }
}

class interval_quadrature_points {
private:
    std::vector<double> points_;
    const double* weights_;
    double scale_;

public:
    using point = double;
    using point_index = simple_index;
    using point_iterator = simple_index_iterator;
    using point_range = simple_index_range;

    interval_quadrature_points(std::vector<double> points, const double* weights, double scale)
    : points_{std::move(points)}
    , weights_{weights}
    , scale_{scale} { }

    auto points() const noexcept -> const std::vector<double>& { return points_; }

    auto indices() const noexcept -> point_range {
        const auto count = narrow_cast<int>(points_.size());
        return range(0, count);
    }

    auto coords(point_index q) const noexcept -> point {
        assert(q >= 0 && q < as_signed(points_.size()) && "Quadrature point index out of bounds");
        return points_[q];
    }

    auto weight(point_index q) const noexcept -> double {
        assert(q >= 0 && q < as_signed(points_.size()) && "Quadrature point index out of bounds");
        return weights_[q] * scale_;
    }

    struct point_data {
        point x;
        double weight;
    };

    auto data(point_index q) const noexcept -> point_data {
        const auto x = coords(q);
        const auto w = weight(q);
        return {x, w};
    }
};

class tensor_quadrature_points {
private:
    interval_quadrature_points ptx_;
    interval_quadrature_points pty_;

public:
    using point = std::tuple<double, double>;
    using point_index = index_types::index;
    using point_iterator = index_types::index_iterator;
    using point_range = index_types::index_range;

    tensor_quadrature_points(interval_quadrature_points ptx,
                             interval_quadrature_points pty) noexcept
    : ptx_{std::move(ptx)}
    , pty_{std::move(pty)} { }

    auto xs() const noexcept -> const std::vector<double>& { return ptx_.points(); }

    auto ys() const noexcept -> const std::vector<double>& { return pty_.points(); }

    auto indices() const noexcept -> point_range {
        const auto rx = ptx_.indices();
        const auto ry = pty_.indices();
        return util::product_range<point_index>(rx, ry);
    }

    auto coords(point_index q) const noexcept -> point {
        const auto [ix, iy] = q;
        const auto x = ptx_.coords(ix);
        const auto y = pty_.coords(iy);
        return {x, y};
    }

    auto weight(point_index q) const noexcept -> double {
        const auto [ix, iy] = q;
        const auto wx = ptx_.weight(ix);
        const auto wy = pty_.weight(iy);
        return wx * wy;
    }

    struct point_data {
        point x;
        double weight;
    };

    auto data(point_index q) const noexcept -> point_data {
        const auto x = coords(q);
        const auto w = weight(q);
        return {x, w};
    }
};

class edge_quadrature_points {
private:
    interval_quadrature_points points_;
    double position_;
    orientation direction_;

public:
    using point = std::tuple<double, double>;
    using point_index = interval_quadrature_points::point_index;
    using point_iterator = interval_quadrature_points::point_iterator;
    using point_range = interval_quadrature_points::point_range;

    edge_quadrature_points(interval_quadrature_points points, double position,
                           orientation direction) noexcept
    : points_{std::move(points)}
    , position_{position}
    , direction_{direction} { }

    auto points() const noexcept -> const std::vector<double>& { return points_.points(); }

    auto position() const noexcept -> double { return position_; }

    auto indices() const noexcept -> point_range { return points_.indices(); }

    auto coords(point_index q) const noexcept -> point {
        const auto s = points_.coords(q);
        if (direction_ == orientation::horizontal) {
            return {s, position_};
        } else {
            return {position_, s};
        }
    }

    auto weight(point_index q) const noexcept -> double { return points_.weight(q); }

    struct point_data {
        point x;
        double weight;
    };

    auto data(point_index q) const noexcept -> point_data {
        const auto x = coords(q);
        const auto w = weight(q);
        return {x, w};
    }
};

class quadrature {
private:
    regular_mesh* mesh_;
    int point_count_;

public:
    using point = regular_mesh::point;
    using element_index = regular_mesh::element_index;
    using facet_index = regular_mesh::facet_index;
    using point_set = tensor_quadrature_points;
    using edge_point_set = edge_quadrature_points;

    quadrature(regular_mesh* mesh, int point_count)
    : mesh_{mesh}
    , point_count_{point_count} {
        assert(point_count_ >= 2 && "Too few quadrature points");
    }

    auto coordinates(element_index e) const -> point_set {
        const auto element = mesh_->element(e);

        auto ptx = data_for_interval(element.span_x);
        auto pty = data_for_interval(element.span_y);

        return {std::move(ptx), std::move(pty)};
    }

    auto coordinates(facet_index f) const -> edge_point_set {
        const auto edge = mesh_->facet(f);
        auto pts = data_for_interval(edge.span);

        return {std::move(pts), edge.position, edge.direction};
    }

private:
    auto data_for_interval(interval target) const -> interval_quadrature_points {
        const auto size = length(target);
        const auto scale = size / 2;  // Gauss quadrature is defined for [-1, 1]
        const auto* weights = quad::gauss::Ws[point_count_];

        return {transform_points(target), weights, scale};
    }

    auto transform_points(interval target) const -> std::vector<double> {
        auto points = std::vector<double>(point_count_);

        for (int i = 0; i < point_count_; ++i) {
            const auto t = quad::gauss::Xs[point_count_][i];  // [-1, 1]
            const auto s = (t + 1) / 2;                       // [ 0, 1]
            points[i] = lerp(s, target);
        }
        return points;
    }
};

using value_type = ads::function_value_2d;

template <typename Value>
struct facet_value {
    Value avg;
    Value jump;
};

template <typename Value>
auto avg(const facet_value<Value>& val) noexcept -> Value {
    return val.avg;
}

template <typename Value>
auto jump(const facet_value<Value>& val) noexcept -> Value {
    return val.jump;
}

auto grad_dot(const value_type& u, const value_type& v) noexcept -> double {
    return u.dx * v.dx + u.dy * v.dy;
}

using point_t = regular_mesh::point;
auto dot(const point_t& a, const point_t& b) noexcept -> double {
    return std::get<0>(a) * std::get<0>(b) + std::get<1>(a) * std::get<1>(b);
}

auto grad(const value_type& v) noexcept -> point_t {
    return {v.dx, v.dy};
}

template <typename ValsX, typename ValsY>
inline auto eval_tensor_basis(const ValsX& vals_x, const ValsY& vals_y) noexcept -> value_type {
    const auto Bx = vals_x(0);
    const auto By = vals_y(0);
    const auto dBx = vals_x(1);
    const auto dBy = vals_y(1);

    const auto v = Bx * By;
    const auto dx = dBx * By;
    const auto dy = Bx * dBy;

    return {v, dx, dy};
}

class space {
private:
    regular_mesh* mesh_;
    bspline_space space_x_;
    bspline_space space_y_;
    global_dof dof_offset_;

    class evaluator;
    class edge_evaluator;

public:
    using point = regular_mesh::point;
    using element_index = regular_mesh::element_index;
    using facet_index = regular_mesh::facet_index;
    using dof_index = index_types::index;
    using dof_iterator = index_types::index_iterator;
    using dof_range = index_types::index_range;

    space(regular_mesh* mesh, bspline::basis bx, bspline::basis by,
          global_dof dof_offset = 0) noexcept
    : mesh_{mesh}
    , space_x_{std::move(bx)}
    , space_y_{std::move(by)}
    , dof_offset_{dof_offset} { }

    auto mesh() const noexcept -> const regular_mesh& { return *mesh_; }

    auto space_x() const noexcept -> const bspline_space& { return space_x_; }

    auto space_y() const noexcept -> const bspline_space& { return space_y_; }

    auto dof_count() const noexcept -> int {
        const auto nx = space_x_.dof_count();
        const auto ny = space_y_.dof_count();

        return nx * ny;
    }

    auto dof_count(element_index e) const noexcept -> int {
        const auto [ex, ey] = e;
        const auto nx = space_x_.dof_count(ex);
        const auto ny = space_y_.dof_count(ey);

        return nx * ny;
    }

    auto facet_dof_count(facet_index f) const noexcept -> int {
        const auto [fx, fy, dir] = f;

        if (dir == orientation::horizontal) {
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.facet_dof_count(fy);
            return ndofs_x * ndofs_y;
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            const auto ndofs_x = space_x_.facet_dof_count(fx);
            const auto ndofs_y = space_y_.dofs_per_element();
            return ndofs_x * ndofs_y;
        }
    }

    auto dofs() const noexcept -> dof_range {
        const auto dofs_x = space_x_.dofs();
        const auto dofs_y = space_y_.dofs();

        return util::product_range<dof_index>(dofs_x, dofs_y);
    }

    auto dofs(element_index e) const noexcept -> dof_range {
        const auto [ex, ey] = e;
        const auto dofs_x = space_x_.dofs(ex);
        const auto dofs_y = space_y_.dofs(ey);

        return util::product_range<dof_index>(dofs_x, dofs_y);
    }

    auto dofs_on_facet(facet_index f) const noexcept -> dof_range {
        const auto [ix, iy, dir] = f;

        if (dir == orientation::horizontal) {
            const auto dofs_x = space_x_.dofs(ix);
            const auto dofs_y = space_y_.dofs_on_facet(iy);
            return util::product_range<dof_index>(dofs_x, dofs_y);
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            const auto dofs_x = space_x_.dofs_on_facet(ix);
            const auto dofs_y = space_y_.dofs(iy);
            return util::product_range<dof_index>(dofs_x, dofs_y);
        }
    }

    auto local_index(dof_index dof, element_index e) const noexcept -> local_dof {
        const auto idx = index_on_element(dof, e);

        const auto ndofs_x = space_x_.dofs_per_element();
        const auto ndofs_y = space_y_.dofs_per_element();

        return linearized(idx, {ndofs_x, ndofs_y});
    }

    auto facet_local_index(dof_index dof, facet_index f) const noexcept -> local_dof {
        const auto idx = index_on_facet(dof, f);

        const auto [fx, fy, dir] = f;

        if (dir == orientation::horizontal) {
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.facet_dof_count(fy);
            return linearized(idx, {ndofs_x, ndofs_y});
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            const auto ndofs_x = space_x_.facet_dof_count(fx);
            const auto ndofs_y = space_y_.dofs_per_element();
            return linearized(idx, {ndofs_x, ndofs_y});
        }
    }

    auto global_index(dof_index dof) const noexcept -> global_dof {
        const auto ndofs_x = space_x_.dof_count();
        const auto ndofs_y = space_y_.dof_count();

        return dof_offset_ + linearized(dof, {ndofs_x, ndofs_y});
    }

    auto dof_evaluator(element_index e, const tensor_quadrature_points& points, int ders) const
        -> evaluator {
        auto data_x = evaluate_basis(points.xs(), space_x_, ders);
        auto data_y = evaluate_basis(points.ys(), space_y_, ders);

        return evaluator{this, e, std::move(data_x), std::move(data_y)};
    }

    auto dof_evaluator(facet_index f, const edge_quadrature_points& points, int ders) const
        -> edge_evaluator {
        const auto [fx, fy, dir] = f;
        if (dir == orientation::horizontal) {
            auto data_x = evaluate_basis(points.points(), space_x_, ders);
            auto data_y = evaluate_basis(fy, space_y_, ders);
            return edge_evaluator{this, f, std::move(data_x), std::move(data_y)};
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            auto data_x = evaluate_basis(fx, space_x_, ders);
            auto data_y = evaluate_basis(points.points(), space_y_, ders);
            return edge_evaluator{this, f, std::move(data_y), std::move(data_x)};
        }
    }

private:
    class evaluator {
    private:
        const space* space_;
        element_index element_;
        bspline_basis_values vals_x_;
        bspline_basis_values vals_y_;

    public:
        using point_index = tensor_quadrature_points::point_index;

        evaluator(const space* space, element_index element, bspline_basis_values vals_x,
                  bspline_basis_values vals_y) noexcept
        : space_{space}
        , element_{element}
        , vals_x_{std::move(vals_x)}
        , vals_y_{std::move(vals_y)} { }

        auto operator()(dof_index dof, point_index q) const noexcept -> value_type {
            const auto [qx, qy] = q;
            const auto [ix, iy] = space_->index_on_element(dof, element_);

            return eval_tensor_basis(
                [&, qx = qx, ix = ix](int der) { return vals_x_(qx, ix, der); },
                [&, qy = qy, iy = iy](int der) { return vals_y_(qy, iy, der); });
        }
    };

    class edge_evaluator {
    private:
        const space* space_;
        facet_index facet_;
        bspline_basis_values vals_interval_;
        bspline_basis_values_on_vertex vals_point_;

    public:
        using point_index = edge_quadrature_points::point_index;

        edge_evaluator(const space* space, facet_index facet, bspline_basis_values vals_interval,
                       bspline_basis_values_on_vertex vals_point) noexcept
        : space_{space}
        , facet_{facet}
        , vals_interval_{std::move(vals_interval)}
        , vals_point_{std::move(vals_point)} { }

        auto operator()(dof_index dof, point_index q, point normal) const noexcept
            -> facet_value<value_type> {
            const auto [ix, iy] = space_->index_on_facet(dof, facet_);
            const auto [nx, ny] = normal;

            if (facet_.dir == orientation::horizontal) {
                const auto avg = eval_tensor_basis(
                    [&, ix = ix](int der) { return vals_interval_(q, ix, der); },
                    [&, iy = iy](int der) { return vals_point_.average(iy, der); });

                const auto jump = eval_tensor_basis(
                    [&, ix = ix](int der) { return vals_interval_(q, ix, der); },
                    [&, iy = iy, ny = ny](int der) { return vals_point_.jump(iy, der, ny); });

                return {avg, jump};
            } else {
                const auto avg = eval_tensor_basis(
                    [&, ix = ix](int der) { return vals_point_.average(ix, der); },
                    [&, iy = iy](int der) { return vals_interval_(q, iy, der); });

                const auto jump = eval_tensor_basis(
                    [&, ix = ix, nx = nx](int der) { return vals_point_.jump(ix, der, nx); },
                    [&, iy = iy](int der) { return vals_interval_(q, iy, der); });

                return {avg, jump};
            }
        }
    };

    auto index_on_element(dof_index dof, element_index e) const noexcept -> dof_index {
        const auto [ex, ey] = e;
        const auto [dx, dy] = dof;

        const auto ix = space_x_.local_index(dx, ex);
        const auto iy = space_y_.local_index(dy, ey);

        return {ix, iy};
    }

    auto index_on_facet(dof_index dof, facet_index f) const noexcept -> dof_index {
        const auto [fx, fy, dir] = f;
        const auto [dx, dy] = dof;

        if (dir == orientation::horizontal) {
            const auto ix = space_x_.local_index(dx, fx);
            const auto iy = space_y_.facet_local_index(dy, fy);
            return {ix, iy};
        } else {
            assert(dir == orientation::vertical && "Invalid edge orientation");
            const auto ix = space_x_.facet_local_index(dx, fx);
            const auto iy = space_y_.local_index(dy, fy);
            return {ix, iy};
        }
    }

    auto linearized(dof_index dof, std::array<int, 2> bounds) const noexcept -> simple_index {
        const auto [ix, iy] = dof;
        const auto [nx, ny] = bounds;
        const auto order = reverse_ordering<2>{{nx, ny}};
        return order.linear_index(ix, iy);
    }
};

class bspline_function {
private:
    const space* space_;
    const double* coefficients_;
    mutable bspline::eval_ctx ctx_x_;
    mutable bspline::eval_ctx ctx_y_;
    mutable std::mutex ctx_lock_;

public:
    using point = space::point;

    bspline_function(const space* space, const double* coefficients)
    : space_{space}
    , coefficients_{coefficients}
    , ctx_x_{space->space_x().degree()}
    , ctx_y_{space->space_y().degree()} { }

    auto operator()(point p) const noexcept -> double { return eval_(p); }

    auto operator()(point p, const regular_mesh::edge_data& edge) const noexcept
        -> facet_value<double> {
        // TODO: implement this correctly
        if (edge.type == facet_type::interior) {
            const auto [px, py] = p;
            const auto [nx, ny] = edge.normal;
            const auto eps = 1e-16;
            const auto before = point{px - eps * nx, py - eps * ny};
            const auto after = point{px + eps * nx, py + eps * ny};

            const auto v0 = eval_(before);
            const auto v1 = eval_(after);
            return {0.5 * (v0 + v1), v0 - v1};
        } else {
            const auto v = eval_(p);
            return {v, v};
        }
    }

private:
    auto eval_(point p) const noexcept -> double {
        const auto [x, y] = p;

        auto coeffs = [this](int i, int j) {
            const auto idx = space_->global_index({i, j});
            return coefficients_[idx];
        };

        const auto& bx = space_->space_x().basis();
        const auto& by = space_->space_y().basis();

        std::scoped_lock guard{ctx_lock_};
        return bspline::eval(x, y, coeffs, bx, by, ctx_x_, ctx_y_);
    }
};

auto evenly_spaced(double a, double b, int elems) -> partition {
    assert(elems > 0 && "Invalid number of partition elements");
    auto xs = partition(elems + 1);

    for (int i = 0; i <= elems; ++i) {
        xs[i] = lerp(i, elems, a, b);
    }
    return xs;
}

auto make_bspline_basis(const partition& points, int p, int c) -> bspline::basis {
    const auto n = as_signed(points.size());
    const auto r = p - c;
    auto knot = bspline::knot_vector{};
    knot.reserve((p + 1) + (n - 1) * r + (p + 1));

    auto append = [&knot](int k, double x) { std::fill_n(back_inserter(knot), k, x); };

    append(p + 1, points[0]);
    for (int i = 1; i < n - 1; ++i) {
        append(r, points[i]);
    }
    append(p + 1, points[n - 1]);

    return bspline::basis{std::move(knot), p};
}

struct index_types3 {
    using index = std::tuple<int, int, int>;
    using index_iterator = util::iter_product3<simple_index_iterator, index>;
    using index_range = boost::iterator_range<index_iterator>;
};

enum class orientation3 {
    dir_x,
    dir_y,
    dir_z,
};

class regular_mesh3 {
private:
    interval_mesh mesh_x_;
    interval_mesh mesh_y_;
    interval_mesh mesh_z_;

public:
    using point = std::tuple<double, double, double>;
    using element_index = index_types3::index;
    using element_iterator = index_types3::index_iterator;
    using element_range = index_types3::index_range;

    struct face_index {
        simple_index ix;
        simple_index iy;
        simple_index iz;
        orientation3 dir;
    };

    using facet_index = face_index;

    regular_mesh3(partition xs, partition ys, partition zs) noexcept
    : mesh_x_{std::move(xs)}
    , mesh_y_{std::move(ys)}
    , mesh_z_{std::move(zs)} { }

    auto elements() const noexcept -> element_range {
        const auto rx = mesh_x_.elements();
        const auto ry = mesh_y_.elements();
        const auto rz = mesh_z_.elements();
        return util::product_range<element_index>(rx, ry, rz);
    }

    struct element_data {
        interval span_x;
        interval span_y;
        interval span_z;
    };

    struct face_data {
        // dir x: 1 - y, 2 - z
        // dir y: 1 - x, 2 - z
        // dir z: 1 - x, 2 - y
        interval span1;
        interval span2;
        double position;
        double diameter;
        orientation3 direction;
        facet_type type;
        point normal;
    };

    using facet_data = face_data;

    auto element(element_index e) const noexcept -> element_data {
        const auto [ix, iy, iz] = e;
        const auto sx = mesh_x_.subinterval(ix);
        const auto sy = mesh_y_.subinterval(iy);
        const auto sz = mesh_z_.subinterval(iz);
        return {sx, sy, sz};
    }

    auto facets() const noexcept -> std::vector<facet_index> {
        auto indices = std::vector<facet_index>{};

        for (auto ix : mesh_x_.facets()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_x});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.facets()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_y});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.facets()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_z});
                }
            }
        }
        return indices;
    }

    auto boundary_facets() const noexcept -> std::vector<facet_index> {
        auto indices = std::vector<facet_index>{};

        for (auto ix : mesh_x_.boundary_facets()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_x});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.boundary_facets()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_y});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.boundary_facets()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_z});
                }
            }
        }
        return indices;
    }

    auto interior_facets() const noexcept -> std::vector<facet_index> {
        auto indices = std::vector<facet_index>{};

        for (auto ix : mesh_x_.interior_facets()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_x});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.interior_facets()) {
                for (auto iz : mesh_z_.elements()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_y});
                }
            }
        }
        for (auto ix : mesh_x_.elements()) {
            for (auto iy : mesh_y_.elements()) {
                for (auto iz : mesh_z_.interior_facets()) {
                    indices.push_back({ix, iy, iz, orientation3::dir_z});
                }
            }
        }
        return indices;
    }

    auto facet(face_index f) const noexcept -> face_data {
        const auto [ix, iy, iz, dir] = f;

        if (dir == orientation3::dir_x) {
            // 1 - y, 2 - z
            const auto [x, type, nx] = mesh_x_.facet(ix);
            const auto span_y = mesh_y_.subinterval(iy);
            const auto span_z = mesh_z_.subinterval(iz);
            const auto normal = point{nx, 0, 0};
            const auto diam = diameter(span_y, span_z);
            return {span_y, span_z, x, diam, dir, type, normal};
        } else if (dir == orientation3::dir_y) {
            // 1 - x, 2 - z
            const auto span_x = mesh_x_.subinterval(ix);
            const auto [y, type, ny] = mesh_y_.facet(iy);
            const auto span_z = mesh_z_.subinterval(iz);
            const auto normal = point{0, ny, 0};
            const auto diam = diameter(span_x, span_z);
            return {span_x, span_z, y, diam, dir, type, normal};
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            // 1 - x, 2 - y
            const auto span_x = mesh_x_.subinterval(ix);
            const auto span_y = mesh_y_.subinterval(iy);
            const auto [z, type, nz] = mesh_z_.facet(iz);
            const auto normal = point{0, 0, nz};
            const auto diam = diameter(span_x, span_y);
            return {span_x, span_y, z, diam, dir, type, normal};
        }
    }

private:
    auto diameter(interval a, interval b) const noexcept -> double {
        return std::hypot(length(a), length(b));
    }
};

class tensor_quadrature_points3 {
private:
    interval_quadrature_points ptx_;
    interval_quadrature_points pty_;
    interval_quadrature_points ptz_;

public:
    using point = std::tuple<double, double, double>;
    using point_index = index_types3::index;
    using point_iterator = index_types3::index_iterator;
    using point_range = index_types3::index_range;

    tensor_quadrature_points3(interval_quadrature_points ptx, interval_quadrature_points pty,
                              interval_quadrature_points ptz) noexcept
    : ptx_{std::move(ptx)}
    , pty_{std::move(pty)}
    , ptz_{std::move(ptz)} { }

    auto xs() const noexcept -> const std::vector<double>& { return ptx_.points(); }

    auto ys() const noexcept -> const std::vector<double>& { return pty_.points(); }

    auto zs() const noexcept -> const std::vector<double>& { return ptz_.points(); }

    auto indices() const noexcept -> point_range {
        const auto rx = ptx_.indices();
        const auto ry = pty_.indices();
        const auto rz = ptz_.indices();
        return util::product_range<point_index>(rx, ry, rz);
    }

    auto coords(point_index q) const noexcept -> point {
        const auto [ix, iy, iz] = q;
        const auto x = ptx_.coords(ix);
        const auto y = pty_.coords(iy);
        const auto z = ptz_.coords(iz);
        return {x, y, z};
    }

    auto weight(point_index q) const noexcept -> double {
        const auto [ix, iy, iz] = q;
        const auto wx = ptx_.weight(ix);
        const auto wy = pty_.weight(iy);
        const auto wz = ptz_.weight(iz);
        return wx * wy * wz;
    }

    struct point_data {
        point x;
        double weight;
    };

    auto data(point_index q) const noexcept -> point_data {
        const auto x = coords(q);
        const auto w = weight(q);
        return {x, w};
    }
};

class face_quadrature_points {
private:
    // dir x: 1 - y, 2 - z
    // dir y: 1 - x, 2 - z
    // dir z: 1 - x, 2 - y
    interval_quadrature_points points1_;
    interval_quadrature_points points2_;
    double position_;
    orientation3 direction_;

public:
    using point = std::tuple<double, double, double>;
    using point_index = index_types::index;
    using point_iterator = index_types::index_iterator;
    using point_range = index_types::index_range;

    face_quadrature_points(interval_quadrature_points points1, interval_quadrature_points points2,
                           double position, orientation3 direction) noexcept
    : points1_{std::move(points1)}
    , points2_{std::move(points2)}
    , position_{position}
    , direction_{direction} { }

    auto points1() const noexcept -> const std::vector<double>& { return points1_.points(); }

    auto points2() const noexcept -> const std::vector<double>& { return points2_.points(); }

    auto position() const noexcept -> double { return position_; }

    auto indices() const noexcept -> point_range {
        const auto r1 = points1_.indices();
        const auto r2 = points2_.indices();
        return util::product_range<point_index>(r1, r2);
    }

    auto coords(point_index q) const noexcept -> point {
        const auto [q1, q2] = q;
        const auto s1 = points1_.coords(q1);
        const auto s2 = points2_.coords(q2);

        if (direction_ == orientation3::dir_x) {
            return {position_, s1, s2};
        } else if (direction_ == orientation3::dir_y) {
            return {s1, position_, s2};
        } else {
            return {s1, s2, position_};
        }
    }

    auto weight(point_index q) const noexcept -> double {
        const auto [q1, q2] = q;
        const auto w1 = points1_.weight(q1);
        const auto w2 = points2_.weight(q2);
        return w1 * w2;
    }

    struct point_data {
        point x;
        double weight;
    };

    auto data(point_index q) const noexcept -> point_data {
        const auto x = coords(q);
        const auto w = weight(q);
        return {x, w};
    }
};

class quadrature3 {
private:
    regular_mesh3* mesh_;
    int point_count_;

public:
    using point = regular_mesh3::point;
    using element_index = regular_mesh3::element_index;
    using facet_index = regular_mesh3::facet_index;
    using point_set = tensor_quadrature_points3;
    using face_point_set = face_quadrature_points;

    quadrature3(regular_mesh3* mesh, int point_count) noexcept
    : mesh_{mesh}
    , point_count_{point_count} {
        assert(point_count_ >= 2 && "Too few quadrature points");
    }

    auto coordinates(element_index e) const -> point_set {
        const auto element = mesh_->element(e);

        auto ptx = data_for_interval(element.span_x);
        auto pty = data_for_interval(element.span_y);
        auto ptz = data_for_interval(element.span_z);

        return {std::move(ptx), std::move(pty), std::move(ptz)};
    }

    auto coordinates(facet_index f) const -> face_point_set {
        const auto face = mesh_->facet(f);
        auto pts1 = data_for_interval(face.span1);
        auto pts2 = data_for_interval(face.span2);

        return {std::move(pts1), std::move(pts2), face.position, face.direction};
    }

private:
    // TODO: remove duplication
    auto data_for_interval(interval target) const -> interval_quadrature_points {
        const auto size = length(target);
        const auto scale = size / 2;  // Gauss quadrature is defined for [-1, 1]
        const auto* weights = quad::gauss::Ws[point_count_];

        return {transform_points(target), weights, scale};
    }

    auto transform_points(interval target) const -> std::vector<double> {
        auto points = std::vector<double>(point_count_);

        for (int i = 0; i < point_count_; ++i) {
            const auto t = quad::gauss::Xs[point_count_][i];  // [-1, 1]
            const auto s = (t + 1) / 2;                       // [ 0, 1]
            points[i] = lerp(s, target);
        }
        return points;
    }
};

using value_type3 = ads::function_value_3d;

auto grad_dot(const value_type3& u, const value_type3& v) noexcept -> double {
    return u.dx * v.dx + u.dy * v.dy + u.dz * v.dz;
}

using point3_t = regular_mesh3::point;

auto dot(const point3_t& a, const point3_t& b) noexcept -> double {
    return std::get<0>(a) * std::get<0>(b) + std::get<1>(a) * std::get<1>(b)
         + std::get<2>(a) * std::get<2>(b);
}

auto grad(const value_type3& v) noexcept -> point3_t {
    return {v.dx, v.dy, v.dz};
}

template <typename ValsX, typename ValsY, typename ValsZ>
inline auto eval_tensor_basis(const ValsX& vals_x, const ValsY& vals_y,
                              const ValsZ& vals_z) noexcept -> value_type3 {
    const auto Bx = vals_x(0);
    const auto dBx = vals_x(1);
    const auto By = vals_y(0);
    const auto dBy = vals_y(1);
    const auto Bz = vals_z(0);
    const auto dBz = vals_z(1);

    const auto v = Bx * By * Bz;
    const auto dx = dBx * By * Bz;
    const auto dy = Bx * dBy * Bz;
    const auto dz = Bx * By * dBz;

    return {v, dx, dy, dz};
}

class space3 {
private:
    regular_mesh3* mesh_;
    bspline_space space_x_;
    bspline_space space_y_;
    bspline_space space_z_;
    global_dof dof_offset_;

    class evaluator;
    class face_evaluator;

public:
    using point = regular_mesh3::point;
    using element_index = regular_mesh3::element_index;
    using facet_index = regular_mesh3::facet_index;
    using dof_index = index_types3::index;
    using dof_iterator = index_types3::index_iterator;
    using dof_range = index_types3::index_range;

    space3(regular_mesh3* mesh, bspline::basis bx, bspline::basis by, bspline::basis bz,
           global_dof dof_offset = 0) noexcept
    : mesh_{mesh}
    , space_x_{std::move(bx)}
    , space_y_{std::move(by)}
    , space_z_{std::move(bz)}
    , dof_offset_{dof_offset} { }

    auto mesh() const noexcept -> const regular_mesh3& { return *mesh_; }

    auto space_x() const noexcept -> const bspline_space& { return space_x_; }

    auto space_y() const noexcept -> const bspline_space& { return space_y_; }

    auto space_z() const noexcept -> const bspline_space& { return space_z_; }

    auto dof_count() const noexcept -> int {
        const auto nx = space_x_.dof_count();
        const auto ny = space_y_.dof_count();
        const auto nz = space_z_.dof_count();

        return nx * ny * nz;
    }

    auto dof_count(element_index e) const noexcept -> int {
        const auto [ex, ey, ez] = e;
        const auto nx = space_x_.dof_count(ex);
        const auto ny = space_y_.dof_count(ey);
        const auto nz = space_z_.dof_count(ez);

        return nx * ny * nz;
    }

    auto facet_dof_count(facet_index f) const noexcept -> int {
        const auto [fx, fy, fz, dir] = f;

        if (dir == orientation3::dir_x) {
            const auto ndofs_x = space_x_.facet_dof_count(fx);
            const auto ndofs_y = space_y_.dofs_per_element();
            const auto ndofs_z = space_z_.dofs_per_element();
            return ndofs_x * ndofs_y * ndofs_z;
        } else if (dir == orientation3::dir_y) {
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.facet_dof_count(fy);
            const auto ndofs_z = space_z_.dofs_per_element();
            return ndofs_x * ndofs_y * ndofs_z;
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.dofs_per_element();
            const auto ndofs_z = space_z_.facet_dof_count(fz);
            return ndofs_x * ndofs_y * ndofs_z;
        }
    }

    auto dofs() const noexcept -> dof_range {
        const auto dofs_x = space_x_.dofs();
        const auto dofs_y = space_y_.dofs();
        const auto dofs_z = space_z_.dofs();

        return util::product_range<dof_index>(dofs_x, dofs_y, dofs_z);
    }

    auto dofs(element_index e) const noexcept -> dof_range {
        const auto [ex, ey, ez] = e;
        const auto dofs_x = space_x_.dofs(ex);
        const auto dofs_y = space_y_.dofs(ey);
        const auto dofs_z = space_z_.dofs(ez);

        return util::product_range<dof_index>(dofs_x, dofs_y, dofs_z);
    }

    auto dofs_on_facet(facet_index f) const noexcept -> dof_range {
        const auto [ix, iy, iz, dir] = f;

        if (dir == orientation3::dir_x) {
            const auto dofs_x = space_x_.dofs_on_facet(ix);
            const auto dofs_y = space_y_.dofs(iy);
            const auto dofs_z = space_z_.dofs(iz);
            return util::product_range<dof_index>(dofs_x, dofs_y, dofs_z);
        } else if (dir == orientation3::dir_y) {
            const auto dofs_x = space_x_.dofs(ix);
            const auto dofs_y = space_y_.dofs_on_facet(iy);
            const auto dofs_z = space_z_.dofs(iz);
            return util::product_range<dof_index>(dofs_x, dofs_y, dofs_z);
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            const auto dofs_x = space_x_.dofs(ix);
            const auto dofs_y = space_y_.dofs(iy);
            const auto dofs_z = space_z_.dofs_on_facet(iz);
            return util::product_range<dof_index>(dofs_x, dofs_y, dofs_z);
        }
    }

    auto local_index(dof_index dof, element_index e) const noexcept -> local_dof {
        const auto idx = index_on_element(dof, e);

        const auto ndofs_x = space_x_.dofs_per_element();
        const auto ndofs_y = space_y_.dofs_per_element();
        const auto ndofs_z = space_z_.dofs_per_element();

        return linearized(idx, {ndofs_x, ndofs_y, ndofs_z});
    }

    auto facet_local_index(dof_index dof, facet_index f) const noexcept -> local_dof {
        const auto idx = index_on_facet(dof, f);

        const auto [fx, fy, fz, dir] = f;

        if (dir == orientation3::dir_x) {
            const auto ndofs_x = space_x_.facet_dof_count(fx);
            const auto ndofs_y = space_y_.dofs_per_element();
            const auto ndofs_z = space_z_.dofs_per_element();
            return linearized(idx, {ndofs_x, ndofs_y, ndofs_z});
        } else if (dir == orientation3::dir_y) {
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.facet_dof_count(fy);
            const auto ndofs_z = space_z_.dofs_per_element();
            return linearized(idx, {ndofs_x, ndofs_y, ndofs_z});
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.dofs_per_element();
            const auto ndofs_z = space_z_.facet_dof_count(fz);
            return linearized(idx, {ndofs_x, ndofs_y, ndofs_z});
        }
    }

    auto global_index(dof_index dof) const noexcept -> global_dof {
        const auto ndofs_x = space_x_.dof_count();
        const auto ndofs_y = space_y_.dof_count();
        const auto ndofs_z = space_z_.dof_count();

        return dof_offset_ + linearized(dof, {ndofs_x, ndofs_y, ndofs_z});
    }

    auto dof_evaluator(element_index e, const tensor_quadrature_points3& points, int ders) const
        -> evaluator {
        auto data_x = evaluate_basis(points.xs(), space_x_, ders);
        auto data_y = evaluate_basis(points.ys(), space_y_, ders);
        auto data_z = evaluate_basis(points.zs(), space_z_, ders);

        return evaluator{this, e, std::move(data_x), std::move(data_y), std::move(data_z)};
    }

    auto dof_evaluator(facet_index f, const face_quadrature_points& points, int ders) const
        -> face_evaluator {
        const auto [fx, fy, fz, dir] = f;

        if (dir == orientation3::dir_x) {
            // 1 - y, 2 - z
            auto data_x = evaluate_basis(fx, space_x_, ders);
            auto data_y = evaluate_basis(points.points1(), space_y_, ders);
            auto data_z = evaluate_basis(points.points2(), space_z_, ders);
            return face_evaluator{
                this, f, std::move(data_y), std::move(data_z), std::move(data_x),
            };
        } else if (dir == orientation3::dir_y) {
            // 1 - x, 2 - z
            auto data_x = evaluate_basis(points.points1(), space_x_, ders);
            auto data_y = evaluate_basis(fy, space_y_, ders);
            auto data_z = evaluate_basis(points.points2(), space_z_, ders);
            return face_evaluator{
                this, f, std::move(data_x), std::move(data_z), std::move(data_y),
            };
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            // 1 - x, 2 - y
            auto data_x = evaluate_basis(points.points1(), space_x_, ders);
            auto data_y = evaluate_basis(points.points2(), space_y_, ders);
            auto data_z = evaluate_basis(fz, space_z_, ders);
            return face_evaluator{
                this, f, std::move(data_x), std::move(data_y), std::move(data_z),
            };
        }
    }

private:
    class evaluator {
    private:
        const space3* space_;
        element_index element_;
        bspline_basis_values vals_x_;
        bspline_basis_values vals_y_;
        bspline_basis_values vals_z_;

    public:
        using point_index = tensor_quadrature_points3::point_index;

        evaluator(const space3* space, element_index element, bspline_basis_values vals_x,
                  bspline_basis_values vals_y, bspline_basis_values vals_z) noexcept
        : space_{space}
        , element_{element}
        , vals_x_{std::move(vals_x)}
        , vals_y_{std::move(vals_y)}
        , vals_z_{std::move(vals_z)} { }

        auto operator()(dof_index dof, point_index q) const noexcept -> value_type3 {
            const auto [qx, qy, qz] = q;
            const auto [ix, iy, iz] = space_->index_on_element(dof, element_);

            return eval_tensor_basis(
                [&, qx = qx, ix = ix](int der) { return vals_x_(qx, ix, der); },
                [&, qy = qy, iy = iy](int der) { return vals_y_(qy, iy, der); },
                [&, qz = qz, iz = iz](int der) { return vals_z_(qz, iz, der); });
        }
    };

    class face_evaluator {
    private:
        const space3* space_;
        facet_index facet_;

        // dir x: 1 - y, 2 - z
        // dir y: 1 - x, 2 - z
        // dir z: 1 - x, 2 - y
        bspline_basis_values vals_interval1_;
        bspline_basis_values vals_interval2_;

        bspline_basis_values_on_vertex vals_point_;

    public:
        // using point_index = face_quadrature_points::point_index;
        using point_index = std::tuple<int, int>;

        face_evaluator(const space3* space, facet_index facet, bspline_basis_values vals_interval1,
                       bspline_basis_values vals_interval2,
                       bspline_basis_values_on_vertex vals_point) noexcept
        : space_{space}
        , facet_{facet}
        , vals_interval1_{std::move(vals_interval1)}
        , vals_interval2_{std::move(vals_interval2)}
        , vals_point_{std::move(vals_point)} { }

        auto operator()(dof_index dof, point_index q, point normal) const noexcept
            -> facet_value<value_type3> {
            const auto [ix, iy, iz] = space_->index_on_facet(dof, facet_);
            const auto [nx, ny, nz] = normal;
            const auto [q1, q2] = q;

            if (facet_.dir == orientation3::dir_x) {
                // 1 - y, 2 - z
                const auto val_y = [&, iy = iy, q1 = q1](int der) {
                    return vals_interval1_(q1, iy, der);
                };
                const auto val_z = [&, iz = iz, q2 = q2](int der) {
                    return vals_interval2_(q2, iz, der);
                };
                const auto avg_x = [&, ix = ix](int der) { return vals_point_.average(ix, der); };
                const auto jump_x = [&, ix = ix, nx = nx](int der) {
                    return vals_point_.jump(ix, der, nx);
                };

                const auto avg = eval_tensor_basis(avg_x, val_y, val_z);
                const auto jump = eval_tensor_basis(jump_x, val_y, val_z);

                return {avg, jump};
            } else if (facet_.dir == orientation3::dir_y) {
                // dir y: 1 - x, 2 - z
                const auto val_x = [&, ix = ix, q1 = q1](int der) {
                    return vals_interval1_(q1, ix, der);
                };
                const auto val_z = [&, iz = iz, q2 = q2](int der) {
                    return vals_interval2_(q2, iz, der);
                };
                const auto avg_y = [&, iy = iy](int der) { return vals_point_.average(iy, der); };
                const auto jump_y = [&, iy = iy, ny = ny](int der) {
                    return vals_point_.jump(iy, der, ny);
                };

                const auto avg = eval_tensor_basis(val_x, avg_y, val_z);
                const auto jump = eval_tensor_basis(val_x, jump_y, val_z);

                return {avg, jump};
            } else {
                // dir z: 1 - x, 2 - y
                const auto val_x = [&, ix = ix, q1 = q1](int der) {
                    return vals_interval1_(q1, ix, der);
                };
                const auto val_y = [&, iy = iy, q2 = q2](int der) {
                    return vals_interval2_(q2, iy, der);
                };
                const auto avg_z = [&, iz = iz](int der) { return vals_point_.average(iz, der); };
                const auto jump_z = [&, iz = iz, nz = nz](int der) {
                    return vals_point_.jump(iz, der, nz);
                };

                const auto avg = eval_tensor_basis(val_x, val_y, avg_z);
                const auto jump = eval_tensor_basis(val_x, val_y, jump_z);

                return {avg, jump};
            }
        }
    };

    auto index_on_element(dof_index dof, element_index e) const noexcept -> dof_index {
        const auto [ex, ey, ez] = e;
        const auto [dx, dy, dz] = dof;

        const auto ix = space_x_.local_index(dx, ex);
        const auto iy = space_y_.local_index(dy, ey);
        const auto iz = space_z_.local_index(dz, ez);

        return {ix, iy, iz};
    }

    auto index_on_facet(dof_index dof, facet_index f) const noexcept -> dof_index {
        const auto [fx, fy, fz, dir] = f;
        const auto [dx, dy, dz] = dof;

        if (dir == orientation3::dir_x) {
            const auto ix = space_x_.facet_local_index(dx, fx);
            const auto iy = space_y_.local_index(dy, fy);
            const auto iz = space_z_.local_index(dz, fz);
            return {ix, iy, iz};
        } else if (dir == orientation3::dir_y) {
            const auto ix = space_x_.local_index(dx, fx);
            const auto iy = space_y_.facet_local_index(dy, fy);
            const auto iz = space_z_.local_index(dz, fz);
            return {ix, iy, iz};
        } else {
            assert(dir == orientation3::dir_z && "Invalid face orientation");
            const auto ix = space_x_.local_index(dx, fx);
            const auto iy = space_y_.local_index(dy, fy);
            const auto iz = space_z_.facet_local_index(dz, fz);
            return {ix, iy, iz};
        }
    }

    auto linearized(dof_index dof, std::array<int, 3> bounds) const noexcept -> simple_index {
        const auto [ix, iy, iz] = dof;
        const auto [nx, ny, nz] = bounds;
        const auto order = reverse_ordering<3>{{nx, ny, nz}};
        return order.linear_index(ix, iy, iz);
    }
};

class bspline_function3 {
private:
    const space3* space_;
    const double* coefficients_;
    mutable bspline::eval_ctx ctx_x_;
    mutable bspline::eval_ctx ctx_y_;
    mutable bspline::eval_ctx ctx_z_;
    mutable std::mutex ctx_lock_;

public:
    using point = space3::point;

    bspline_function3(const space3* space, const double* coefficients)
    : space_{space}
    , coefficients_{coefficients}
    , ctx_x_{space->space_x().degree()}
    , ctx_y_{space->space_y().degree()}
    , ctx_z_{space->space_z().degree()} { }

    auto operator()(point p) const noexcept -> double { return eval_(p); }

    auto operator()(point p, const regular_mesh3::face_data& face) const noexcept
        -> facet_value<double> {
        // TODO: implement this correctly
        if (face.type == facet_type::interior) {
            const auto [px, py, pz] = p;
            const auto [nx, ny, nz] = face.normal;
            const auto eps = 1e-16;
            const auto before = point{px - eps * nx, py - eps * ny, pz - eps * nz};
            const auto after = point{px + eps * nx, py + eps * ny, pz + eps * nz};

            const auto v0 = eval_(before);
            const auto v1 = eval_(after);
            return {0.5 * (v0 + v1), v0 - v1};
        } else {
            const auto v = eval_(p);
            return {v, v};
        }
    }

private:
    auto eval_(point p) const noexcept -> double {
        const auto [x, y, z] = p;

        auto coeffs = [this](int i, int j, int k) {
            const auto idx = space_->global_index({i, j, k});
            return coefficients_[idx];
        };

        const auto& bx = space_->space_x().basis();
        const auto& by = space_->space_y().basis();
        const auto& bz = space_->space_z().basis();

        std::scoped_lock guard{ctx_lock_};
        return bspline::eval(x, y, z, coeffs, bx, by, bz, ctx_x_, ctx_y_, ctx_z_);
    }
};

}  // namespace ads

template <typename Space, typename Quad, typename Out, typename Form>
auto assemble(const Space& space, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        const auto eval = space.dof_evaluator(e, points, 1);

        const auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        using dof_idx = typename Space::dof_index;
        using point_idx = typename decltype(points)::point_index;
        using value_type = decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>()));

        auto basis_vals = std::vector<value_type>(n);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs(e)) {
                const auto iloc = space.local_index(i, e);
                basis_vals[iloc] = eval(i, q);
            }
            for (int iloc = 0; iloc < n; ++iloc) {
                for (int jloc = 0; jloc < n; ++jloc) {
                    const auto& u = basis_vals[iloc];
                    const auto& v = basis_vals[jloc];
                    M(jloc, iloc) += form(u, v, x) * w;
                }
            }
        }

        for (auto i : space.dofs(e)) {
            const auto iloc = space.local_index(i, e);
            const auto I = space.global_index(i);
            for (auto j : space.dofs(e)) {
                const auto jloc = space.local_index(j, e);
                const auto J = space.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Trial, typename Test, typename Quad, typename Out, typename Form>
auto assemble(const Trial& trial, const Test& test, const Quad& quad, Out out, Form&& form)
    -> void {
    const auto& mesh = test.mesh();

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        const auto eval_trial = trial.dof_evaluator(e, points, 1);
        const auto eval_test = test.dof_evaluator(e, points, 1);

        const auto n_trial = trial.dof_count(e);
        const auto n_test = test.dof_count(e);
        auto M = ads::lin::tensor<double, 2>{{n_test, n_trial}};

        using point_idx = typename decltype(points)::point_index;
        using trial_dof_idx = typename Trial::dof_index;
        using test_dof_idx = typename Test::dof_index;
        using trial_value =
            decltype(eval_trial(std::declval<trial_dof_idx>(), std::declval<point_idx>()));
        using test_value =
            decltype(eval_test(std::declval<test_dof_idx>(), std::declval<point_idx>()));

        auto test_vals = std::vector<test_value>(n_test);
        auto trial_vals = std::vector<trial_value>(n_trial);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : trial.dofs(e)) {
                const auto iloc = trial.local_index(i, e);
                trial_vals[iloc] = eval_trial(i, q);
            }
            for (auto j : test.dofs(e)) {
                const auto jloc = test.local_index(j, e);
                test_vals[jloc] = eval_test(j, q);
            }
            for (int iloc = 0; iloc < n_trial; ++iloc) {
                for (int jloc = 0; jloc < n_test; ++jloc) {
                    const auto& u = trial_vals[iloc];
                    const auto& v = test_vals[jloc];
                    M(jloc, iloc) += form(u, v, x) * w;
                }
            }
        }

        for (auto i : trial.dofs(e)) {
            const auto iloc = trial.local_index(i, e);
            const auto I = trial.global_index(i);
            for (auto j : test.dofs(e)) {
                const auto jloc = test.local_index(j, e);
                const auto J = test.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Facets, typename Space, typename Quad, typename Out, typename Form>
auto assemble_facets(const Facets& facets, const Space& space, const Quad& quad, Out out,
                     Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto f : facets) {
        const auto facet = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval = space.dof_evaluator(f, points, 1);

        const auto n = space.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        using dof_idx = typename Space::dof_index;
        using point_idx = typename decltype(points)::point_index;
        using value_type =
            decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>(), facet.normal));

        auto basis_vals = std::vector<value_type>(n);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs_on_facet(f)) {
                const auto iloc = space.facet_local_index(i, f);
                basis_vals[iloc] = eval(i, q, facet.normal);
            }
            for (int iloc = 0; iloc < n; ++iloc) {
                for (int jloc = 0; jloc < n; ++jloc) {
                    const auto& u = basis_vals[iloc];
                    const auto& v = basis_vals[jloc];
                    M(jloc, iloc) += form(u, v, x, facet) * w;
                }
            }
        }

        for (auto i : space.dofs_on_facet(f)) {
            const auto iloc = space.facet_local_index(i, f);
            const auto I = space.global_index(i);
            for (auto j : space.dofs_on_facet(f)) {
                const auto jloc = space.facet_local_index(j, f);
                const auto J = space.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Facets, typename Trial, typename Test, typename Quad, typename Out,
          typename Form>
auto assemble_facets(const Facets& facets, const Trial& trial, const Test& test, const Quad& quad,
                     Out out, Form&& form) -> void {
    const auto& mesh = test.mesh();

    for (auto f : facets) {
        const auto facet = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval_trial = trial.dof_evaluator(f, points, 1);
        const auto eval_test = test.dof_evaluator(f, points, 1);

        const auto n_trial = trial.facet_dof_count(f);
        const auto n_test = test.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n_test, n_trial}};

        using point_idx = typename decltype(points)::point_index;
        using trial_dof_idx = typename Trial::dof_index;
        using test_dof_idx = typename Test::dof_index;
        using trial_value = decltype(eval_trial(std::declval<trial_dof_idx>(),
                                                std::declval<point_idx>(), facet.normal));
        using test_value = decltype(eval_test(std::declval<test_dof_idx>(),
                                              std::declval<point_idx>(), facet.normal));

        auto test_vals = std::vector<test_value>(n_test);
        auto trial_vals = std::vector<trial_value>(n_trial);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : trial.dofs_on_facet(f)) {
                const auto iloc = trial.facet_local_index(i, f);
                trial_vals[iloc] = eval_trial(i, q, facet.normal);
            }
            for (auto j : test.dofs_on_facet(f)) {
                const auto jloc = test.facet_local_index(j, f);
                test_vals[jloc] = eval_test(j, q, facet.normal);
            }
            for (int iloc = 0; iloc < n_trial; ++iloc) {
                for (int jloc = 0; jloc < n_test; ++jloc) {
                    const auto& u = trial_vals[iloc];
                    const auto& v = test_vals[jloc];
                    M(jloc, iloc) += form(u, v, x, facet) * w;
                }
            }
        }

        for (auto i : trial.dofs_on_facet(f)) {
            const auto iloc = trial.facet_local_index(i, f);
            const auto I = trial.global_index(i);
            for (auto j : test.dofs_on_facet(f)) {
                const auto jloc = test.facet_local_index(j, f);
                const auto J = test.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Space, typename Quad, typename Out, typename Form>
auto assemble_rhs(const Space& space, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        const auto eval = space.dof_evaluator(e, points, 1);

        const auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 1>{{n}};

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto j : space.dofs(e)) {
                const auto v = eval(j, q);
                const auto jloc = space.local_index(j, e);
                M(jloc) += form(v, x) * w;
            }
        }

        for (auto j : space.dofs(e)) {
            const auto jloc = space.local_index(j, e);
            const auto J = space.global_index(j);
            out(J, M(jloc));
        }
    }
}

template <typename Facets, typename Space, typename Quad, typename Out, typename Form>
auto assemble_rhs(const Facets& facets, const Space& space, const Quad& quad, Out out, Form&& form)
    -> void {
    const auto& mesh = space.mesh();

    for (auto f : facets) {
        const auto facet = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval = space.dof_evaluator(f, points, 1);

        const auto n = space.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 1>{{n}};

        for (auto q : points.indices()) {
            const auto [X, w] = points.data(q);
            for (auto j : space.dofs_on_facet(f)) {
                const auto v = avg(eval(j, q, facet.normal));
                const auto jloc = space.facet_local_index(j, f);
                M(jloc) += form(v, X, facet) * w;
            }
        }
        for (auto j : space.dofs_on_facet(f)) {
            const auto jloc = space.facet_local_index(j, f);
            const auto J = space.global_index(j);
            out(J, M(jloc));
        }
    }
}

template <typename Mesh, typename Quad, typename Function>
auto integrate(const Mesh& mesh, const Quad& quad, Function&& f) {
    using value_of_f = decltype(f(std::declval<typename Quad::point>()));
    auto a = value_of_f{};

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        for (auto q : points.indices()) {
            auto [x, w] = points.data(q);
            a += f(x) * w;
        }
    }
    return a;
}

template <typename Mesh, typename Quad, typename Norm, typename Function>
auto norm(const Mesh& mesh, const Quad& quad, Norm&& norm, Function&& f) -> double {
    const auto g = [&f, &norm](auto x) { return norm(f(x)); };
    const auto value = integrate(mesh, quad, g);
    return std::sqrt(value);
}

template <typename Mesh, typename Quad, typename Norm, typename Function, typename Exact>
auto error(const Mesh& mesh, const Quad& quad, Norm&& norm, Function&& f, Exact&& exact) -> double {
    const auto difference = [&f, &exact](auto x) { return f(x) - exact(x); };
    return ::norm(mesh, quad, std::forward<Norm>(norm), difference);
}

struct L2 {
    constexpr auto operator()(double v) const noexcept -> double { return v * v; }

    constexpr auto operator()(const ads::value_type& v) const noexcept -> double {
        const auto& self = *this;
        return self(v.val);
    }
};

namespace detail {

template <typename F>
struct validity_checker {
    template <typename... Ts>
    constexpr auto operator()(Ts&&...) const {
        return std::is_invocable<F, Ts...>{};
    }
};

template <typename F>
constexpr auto is_valid(F) {
    return validity_checker<F>{};
}

template <typename T>
inline constexpr bool has_val = detail::is_valid([](auto&& x) -> decltype(x.val) {})(T{});

}  // namespace detail

// derivative of order X
template <int K>
struct derivatives_of_order { };

// edge value
template <typename Val>
struct facet_discontinuity { };

struct element_integral { };
struct boundary_integral { };
struct facet_integral { };
struct interior_facet_integral { };

struct H10 {
    using required_data = derivatives_of_order<1>;
    using computation_method = element_integral;

    auto operator()(const ads::value_type& v) const noexcept -> double {
        using ads::dot, ads::grad;
        return dot(grad(v), grad(v));
    }
};

struct H1 {
    using required_data = derivatives_of_order<1>;
    using computation_method = element_integral;

    auto operator()(const ads::value_type& v) const noexcept -> double {
        using ads::dot, ads::grad;
        return v.val * v.val + dot(grad(v), grad(v));
    }
};

struct L2_2 {
    using required_data = derivatives_of_order<0>;
    using computation_method = element_integral;

    auto operator()(double v) const noexcept -> double { return v * v; }
};

template <typename Function, typename Arg, typename InputTag>
auto eval(Function&& f, const Arg& x, InputTag) {
    return f(x);
}

template <                                                             //
    typename Function,                                                 //
    typename Arg,                                                      //
    typename = std::enable_if_t<                                       //
        !std::is_convertible_v<                                        //
            std::invoke_result_t<Function, Arg>,                       //
            double                                                     //
            > && detail::has_val<std::invoke_result_t<Function, Arg>>  //
        >                                                              //
    >
auto eval(Function&& f, const Arg& x, derivatives_of_order<0>) -> double {
    const auto v = f(x);
    return v.val;
}

template <typename Mesh, typename Quad, typename Norm, typename Function>
auto alt_norm(const Mesh& mesh, const Quad& quad, Norm&& norm, Function&& f) -> double {
    using norm_input = typename Norm::required_data;

    const auto g = [&](auto x) {
        auto value = eval(std::forward<Function>(f), x, norm_input{});
        static_assert(std::is_invocable_v<Norm, decltype(value)>,
                      "Invalid value type for the specified norm");
        return norm(value);
    };
    const auto value = integrate(mesh, quad, g);
    return std::sqrt(value);
}

template <typename Facets, typename Mesh, typename Quad, typename Function>
auto integrate_facets(const Facets& facets, const Mesh& mesh, const Quad& quad, Function&& fun) {
    using value_of_f = decltype(fun(std::declval<typename Quad::point>(),
                                    std::declval<typename Mesh::facet_data>()));
    auto a = value_of_f{};

    for (auto f : facets) {
        const auto facet = mesh.facet(f);
        const auto points = quad.coordinates(f);
        for (auto q : points.indices()) {
            auto [x, w] = points.data(q);
            a += fun(x, facet) * w;
        }
    }
    return a;
}

template <typename... Args>
auto sum_norms(Args... args) noexcept -> double {
    const auto sum = ((args * args) + ...);
    return std::sqrt(sum);
}

template <typename... Funs>
auto save_to_file(const std::string& path, Funs&&... funs) -> void {
    constexpr auto res = 200;
    auto out = fmt::output_file(path);

    for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            out.print("{} {}", x, y);
            (out.print(" {:.7}", funs({x, y})), ...);
            out.print("\n");
        }
    }
}

#endif  // ADS_EXPERIMENTAL_ALL_HPP
