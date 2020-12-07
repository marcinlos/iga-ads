#include <vector>
#include <tuple>
#include <optional>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <type_traits>
#include <chrono>
#include <mutex>
#include <fmt/chrono.h>
#include <fmt/os.h>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/counting_range.hpp>

#include "ads/util.hpp"
#include "ads/util/iter/product.hpp"
#include "ads/bspline/bspline.hpp"
#include "ads/bspline/eval.hpp"
#include "ads/quad/gauss.hpp"
#include "ads/lin/tensor.hpp"
#include "ads/util/function_value.hpp"
#include "mumps.hpp"
#include "ads/executor/galois.hpp"
#include "ads/executor/sequential.hpp"


namespace ads {

    template <typename T>
    auto as_signed(T a) {
        return std::make_signed_t<T>(a);
    }

    template <typename T>
    auto as_unsigned(T a) {
        return std::make_unsigned_t<T>(a);
    }


    using partition = std::vector<double>;

    using simple_index          = int;
    using simple_index_iterator = boost::counting_iterator<simple_index>;
    using simple_index_range    = boost::iterator_range<simple_index_iterator>;

    auto range(int start, int past_end) noexcept -> simple_index_range {
        return boost::counting_range(start, past_end);
    }

    auto empty_range() noexcept -> simple_index_range {
        return range(0, 0);
    }

    struct index_types {
        using index          = std::tuple<int, int>;
        using index_iterator = util::iter_product2<simple_index_iterator, index>;
        using index_range    = boost::iterator_range<index_iterator>;
    };


    struct interval {
        double left;
        double right;

        constexpr interval(double left, double right) noexcept
        : left{left}, right{right}
        { }
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

    auto operator <<(std::ostream& os, interval s) -> std::ostream& {
        return os << "[" << s.left << ", " << s.right << "]";
    }

    enum class facet_type {
        interior = true,
        boundary = false
    };

    class interval_mesh {
    private:
        partition points_;

    public:
        using point            = double;
        using element_index    = simple_index;
        using element_iterator = simple_index_iterator;
        using element_range    = simple_index_range;
        using facet_index      = simple_index;
        using facet_range      = simple_index_range;

        explicit interval_mesh(partition points) noexcept
        : points_{std::move(points)}
        { }

        auto elements() const noexcept -> element_range {
            return range(0, element_count());
        }

        auto element_count() const noexcept -> int {
            return points_.size() - 1;
        }

        auto subinterval(element_index e) const noexcept -> interval {
            return ads::subinterval(points_, e);
        }

        struct point_data {
            point        position;
            facet_type   type;
            double       normal;
        };

        auto facets() const noexcept -> facet_range {
            const auto facet_count = points_.size();
            return range(0, facet_count);
        }

        auto boundary_facets() const noexcept -> std::array<facet_index, 2> {
            const auto last = points_.size() - 1;
            return {0, last};
        }

        auto interior_facets() const noexcept -> facet_range {
            const auto facet_count = points_.size();
            return range(1, facet_count - 1);
        }

        auto facet(facet_index i) const noexcept -> point_data {
            assert(i < as_signed(points_.size()) && "Point index out of range");
            // all points are positive except for the first one
            const auto normal = i > 0 ? 1.0 : -1.0;
            const auto type   = (i > 0 && i < as_signed(points_.size()) - 1) ? facet_type::interior
                                                                             : facet_type::boundary;
            return {points_[i], type, normal};
        }

    };

    enum class orientation {
        horizontal,
        vertical
    };

    class regular_mesh {
    private:
        interval_mesh mesh_x_;
        interval_mesh mesh_y_;

    public:
        using point            = std::tuple<double, double>;
        using element_index    = index_types::index;
        using element_iterator = index_types::index_iterator;
        using element_range    = index_types::index_range;

        struct edge_index {
            simple_index ix;
            simple_index iy;
            orientation  dir;
        };

        using facet_index      = edge_index;
        // using facet_iterator   = void;
        // using facet_range      = void;

        regular_mesh(partition xs, partition ys) noexcept
        : mesh_x_{std::move(xs)}
        , mesh_y_{std::move(ys)}
        { }

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
            interval    span;
            double      position;
            orientation direction;
            facet_type  type;
            point       normal;
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
                const auto span          = mesh_x_.subinterval(ix);
                const auto normal        = point{0, ny};
                return {span, y, dir, type, normal};
            } else {
                assert(dir == orientation::vertical && "Invalid edge orientation");
                const auto [x, type, nx] = mesh_x_.facet(ix);
                const auto span          = mesh_y_.subinterval(iy);
                const auto normal        = point{nx, 0};
                return {span, x, dir, type, normal};
            }
        }

    };

    auto operator <<(std::ostream& os, const regular_mesh::edge_index& edge) -> std::ostream& {
        const auto [ix, iy, dir] = edge;
        const char* sign = dir == orientation::horizontal ? "-" : "|";
        return os << "(" << ix << ", " << iy << ")[" << sign << "]";
    }


    using local_dof  = simple_index;
    using global_dof = simple_index;


    auto spans_for_elements(const bspline::basis& b) -> std::vector<int> {
        auto spans = std::vector<int>{};
        spans.reserve(b.elements());

        for (int i = 0; i + 1 < as_signed(b.knot_size()); ++ i) {
            if (b.knot[i] != b.knot[i + 1]) {
                spans.push_back(i);
            }
        }
        assert(as_signed(spans.size()) == b.elements());
        return spans;
    }


    class bspline_space {
    private:
        bspline::basis   basis_;
        std::vector<int> first_dofs_;
        std::vector<int> spans_;

    public:
        using point         = double;
        using element_index = simple_index;
        using facet_index   = simple_index;
        using dof_index     = simple_index;
        using dof_iterator  = simple_index_iterator;
        using dof_range     = simple_index_range;

        explicit bspline_space(bspline::basis basis)
        : basis_{std::move(basis)}
        , first_dofs_{first_nonzero_dofs(basis_)}
        , spans_{spans_for_elements(basis_)}
        { }

        auto basis() const noexcept -> const bspline::basis& {
            return basis_;
        }

        auto degree() const noexcept -> int {
            return basis_.degree;
        }

        auto dofs_per_element() const noexcept -> int {
            return degree() + 1;
        }

        auto dof_count() const noexcept -> int {
            return basis_.dofs();
        }

        auto dof_count(element_index) const noexcept -> int {
            return dofs_per_element();
        }

        auto facet_dof_count(facet_index f) const noexcept -> int {
            auto dofs = dofs_on_facet(f);
            using std::begin;
            using std::end;
            return std::distance(begin(dofs), end(dofs));
        }

        auto dofs() const noexcept -> dof_range {
            return range(0, dof_count());
        }

        auto dofs(element_index e) const noexcept -> dof_range {
            const auto first = first_dofs_[e];
            return range(first, first + dofs_per_element());
        }

        auto local_index(dof_index dof, element_index e) const noexcept -> local_dof {
            const auto first = first_dofs_[e];
            return dof - first;
        }

        auto first_dof(element_index e) const noexcept -> global_dof {
            return first_dofs_[e];
        }

        auto last_dof(element_index e) const noexcept -> global_dof {
            return first_dof(e) + dofs_per_element() - 1;
        }

        auto dofs_on_facet(facet_index f) const noexcept -> dof_range {
            const auto last_element  = basis_.elements() - 1;
            const auto elem_left     = std::max(f - 1, 0);
            const auto elem_right    = std::min(f, last_element);
            const auto first         = first_dofs_[elem_left];
            const auto one_past_last = first_dofs_[elem_right] + dofs_per_element();
            return range(first, one_past_last);
        }

        auto facet_local_index(dof_index dof, facet_index f) const noexcept -> local_dof {
            const auto elem_left = std::max(f - 1, 0);
            const auto first     = first_dofs_[elem_left];
            return dof - first;
        }

        auto span(element_index e) const noexcept -> int {
            return spans_[e];
        }
    };


    class bspline_basis_values {
    private:
        std::vector<double>   buffer_;
        std::vector<double**> point_data_; // indexed by point
        std::vector<double*>  dof_data_;   // indexed by point and derivative

    public:

        bspline_basis_values(int points, int dofs, int ders)
        : buffer_(points * dofs * (ders + 1))
        , point_data_(points)
        , dof_data_(points * (ders + 1)) {

            for (int i = 0; i < points * (ders + 1); ++ i) {
                dof_data_[i] = &buffer_[i * dofs];
            }
            for (int i = 0; i < points; ++ i) {
                point_data_[i] = &dof_data_[i * (ders + 1)];
            }
        }

        auto operator ()(int point, local_dof i, int der) const -> double {
            return point_data_[point][der][i];
        }

        auto point_buffer(int point) noexcept -> double** {
            return point_data_[point];
        }
    };


    auto evaluate_basis(const std::vector<double>& points, const bspline_space& space, int ders) -> bspline_basis_values {
        const auto point_count = static_cast<int>(points.size());
        const auto dof_count   = space.dofs_per_element();

        auto values  = bspline_basis_values{point_count, dof_count, ders};
        auto context = bspline::eval_ctx{space.degree()};

        for (int q = 0; q < as_signed(points.size()); ++ q) {
            const auto x      = points[q];
            const auto buffer = values.point_buffer(q);
            const auto span   = find_span(x, space.basis());
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
                                       std::optional<bspline_basis_values> right,
                                       local_dof left_last, local_dof right_first) noexcept
        : left_{std::move(left)}
        , right_{std::move(right)}
        , left_last_{left_last}
        , right_first_{right_first} {
            assert((left_ || right_) && "Neither left nor right adjacent element data specified");
        }

        auto operator ()(local_dof i, int der) const noexcept -> double {
            if (left_ && i <= left_last_) {
                return left(i, der);
            } else { // right_ has value
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
            const auto left_val  = left(i, der);
            const auto right_val = right(i, der);
            return normal * (left_val - right_val);
        }

        auto average(local_dof i, int der) const noexcept -> double {
            const auto left_val  = left(i, der);
            const auto right_val = right(i, der);
            const auto sum       = left_val + right_val;
            if (left_ && right_) {
                return sum / 2;
            } else {
                // one of these is zero
                return sum;
            }
        }
    };


    auto evaluate_basis_at_point(double x, const bspline_space& space, int ders, int span) -> bspline_basis_values {
        const auto dof_count = space.dofs_per_element();

        auto values  = bspline_basis_values{1, dof_count, ders};
        auto context = bspline::eval_ctx{space.degree()};

        const auto buffer = values.point_buffer(0);

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
        const auto  x     = basis.points[f];

        const auto maybe_elem_left  = element_left(f, space.basis());
        const auto maybe_elem_right = element_right(f, space.basis());

        if (maybe_elem_left && maybe_elem_right) {
            const auto elem_left       = maybe_elem_left.value();
            const auto elem_right      = maybe_elem_right.value();
            const auto span_left       = space.span(elem_left);
            const auto span_right      = space.span(elem_right);
            const auto left_last       = space.last_dof(elem_left);
            const auto right_first     = space.first_dof(elem_right);
            const auto left_last_loc   = space.facet_local_index(left_last, f);
            const auto right_first_loc = space.facet_local_index(right_first, f);

            auto vals_left  = evaluate_basis_at_point(x, space, ders, span_left);
            auto vals_right = evaluate_basis_at_point(x, space, ders, span_right);

            return {std::move(vals_left), std::move(vals_right), left_last_loc, right_first_loc};

        } else if (maybe_elem_left) {
            const auto elem_left     = maybe_elem_left.value();
            const auto span_left     = space.span(elem_left);
            const auto left_last     = space.last_dof(elem_left);
            const auto left_last_loc = space.facet_local_index(left_last, f);

            auto vals_left = evaluate_basis_at_point(x, space, ders, span_left);

            return {std::move(vals_left), {}, left_last_loc, {}};

        } else { // maybe_elem_right
            assert(maybe_elem_right && "No elements adjacent to specified face");
            const auto elem_right      = maybe_elem_right.value();
            const auto span_right      = space.span(elem_right);
            const auto right_first     = space.first_dof(elem_right);
            const auto right_first_loc = space.facet_local_index(right_first, f);

            auto vals_right = evaluate_basis_at_point(x, space, ders, span_right);

            return {{}, std::move(vals_right), {}, right_first_loc};
        }
    }


    class interval_quadrature_points {
    private:
        std::vector<double> points_;
        const double*       weights_;
        double              scale_;

    public:
        using point          = double;
        using point_index    = simple_index;
        using point_iterator = simple_index_iterator;
        using point_range    = simple_index_range;

        interval_quadrature_points(std::vector<double> points, const double* weights, double scale)
        : points_{std::move(points)}
        , weights_{weights}
        , scale_{scale}
        { }

        auto points() const noexcept -> const std::vector<double>& {
            return points_;
        }

        auto indices() const noexcept -> point_range {
            return range(0, points_.size());
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
            point  x;
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
        using point          = std::tuple<double, double>;
        using point_index    = index_types::index;
        using point_iterator = index_types::index_iterator;
        using point_range    = index_types::index_range;

        tensor_quadrature_points(interval_quadrature_points ptx,
                                 interval_quadrature_points pty) noexcept
        : ptx_{std::move(ptx)}
        , pty_{std::move(pty)}
        { }

        auto xs() const noexcept -> const std::vector<double>& {
            return ptx_.points();
        }

        auto ys() const noexcept -> const std::vector<double>& {
            return pty_.points();
        }

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
            point  x;
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
        double                     position_;
        orientation                direction_;

    public:
        using point          = std::tuple<double, double>;
        using point_index    = interval_quadrature_points::point_index;
        using point_iterator = interval_quadrature_points::point_iterator;
        using point_range    = interval_quadrature_points::point_range;

        edge_quadrature_points(interval_quadrature_points points, double position,
                               orientation direction) noexcept
        : points_{std::move(points)}
        , position_{position}
        , direction_{direction}
        { }

        auto points() const noexcept -> const std::vector<double>& {
            return points_.points();
        }

        auto position() const noexcept -> double {
            return position_;
        }

        auto indices() const noexcept -> point_range {
            return points_.indices();
        }

        auto coords(point_index q) const noexcept -> point {
            const auto s = points_.coords(q);
            if (direction_ == orientation::horizontal) {
                return {s, position_};
            } else {
                return {position_, s};
            }
        }

        auto weight(point_index q) const noexcept -> double {
            return points_.weight(q);
        }

        struct point_data {
            point  x;
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
        using point          = regular_mesh::point;
        using element_index  = regular_mesh::element_index;
        using facet_index    = regular_mesh::facet_index;
        using point_set      = tensor_quadrature_points;
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
            const auto size    = length(target);
            const auto scale   = size / 2;        // Gauss quadrature is defined for [-1, 1]
            const auto weights = quad::gauss::Ws[point_count_];

            return {transform_points(target), weights, scale};
        }

        auto transform_points(interval target) const -> std::vector<double> {
            auto points = std::vector<double>(point_count_);

            for (int i = 0; i < point_count_; ++ i) {
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
        const auto  Bx = vals_x(0);
        const auto dBx = vals_x(1);
        const auto  By = vals_y(0);
        const auto dBy = vals_y(1);

        const auto v  =  Bx *  By;
        const auto dx = dBx *  By;
        const auto dy =  Bx * dBy;

        return {v, dx, dy};
    }

    class space {
    private:
        regular_mesh* mesh_;
        bspline_space space_x_;
        bspline_space space_y_;
        global_dof    dof_offset_;

        class evaluator;
        class edge_evaluator;

    public:
        using point         = regular_mesh::point;
        using element_index = regular_mesh::element_index;
        using facet_index   = regular_mesh::facet_index;
        using dof_index     = index_types::index;
        using dof_iterator  = index_types::index_iterator;
        using dof_range     = index_types::index_range;

        space(regular_mesh* mesh, bspline::basis bx, bspline::basis by, global_dof dof_offset = 0) noexcept
        : mesh_{mesh}
        , space_x_{std::move(bx)}
        , space_y_{std::move(by)}
        , dof_offset_{dof_offset}
        { }

        auto mesh() const noexcept -> const regular_mesh& {
            return *mesh_;
        }

        auto space_x() const noexcept -> const bspline_space& {
            return space_x_;
        }

        auto space_y() const noexcept -> const bspline_space& {
            return space_y_;
        }

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


        auto dof_evaluator(element_index e, const tensor_quadrature_points& points, int ders) const -> evaluator {
            auto data_x = evaluate_basis(points.xs(), space_x_, ders);
            auto data_y = evaluate_basis(points.ys(), space_y_, ders);

            return evaluator{this, e, ders, std::move(data_x), std::move(data_y)};
        }

        auto dof_evaluator(facet_index f, const edge_quadrature_points& points, int ders) const -> edge_evaluator {
            const auto [fx, fy, dir] = f;
            if (dir == orientation::horizontal) {
                auto data_x = evaluate_basis(points.points(), space_x_, ders);
                auto data_y = evaluate_basis(fy, space_y_, ders);
                return edge_evaluator{this, f, ders, std::move(data_x), std::move(data_y)};
            } else {
                assert(dir == orientation::vertical && "Invalid edge orientation");
                auto data_x = evaluate_basis(fx, space_x_, ders);
                auto data_y = evaluate_basis(points.points(), space_y_, ders);
                return edge_evaluator{this, f, ders, std::move(data_y), std::move(data_x)};
            }
        }

    private:

        class evaluator {
        private:
            const space*         space_;
            element_index        element_;
            int                  derivatives_;
            bspline_basis_values vals_x_;
            bspline_basis_values vals_y_;

        public:
            using point_index = tensor_quadrature_points::point_index;

            evaluator(const space* space, element_index element, int derivatives,
                      bspline_basis_values vals_x, bspline_basis_values vals_y) noexcept
            : space_{space}
            , element_{element}
            , derivatives_{derivatives}
            , vals_x_{std::move(vals_x)}
            , vals_y_{std::move(vals_y)}
            { }

            auto operator ()(dof_index dof, point_index q) const noexcept -> value_type {
                const auto [qx, qy] = q;
                const auto [ix, iy] = space_->index_on_element(dof, element_);

                return eval_tensor_basis([&,qx=qx,ix=ix](int der) { return vals_x_(qx, ix, der); },
                                         [&,qy=qy,iy=iy](int der) { return vals_y_(qy, iy, der); });
            }
        };

        class edge_evaluator {
        private:
            const space*                   space_;
            facet_index                    facet_;
            int                            derivatives_;
            bspline_basis_values           vals_interval_;
            bspline_basis_values_on_vertex vals_point_;

        public:
            using point_index = edge_quadrature_points::point_index;

            edge_evaluator(const space* space, facet_index facet, int derivatives,
                           bspline_basis_values vals_interval,
                           bspline_basis_values_on_vertex vals_point) noexcept
            : space_{space}
            , facet_{facet}
            , derivatives_{derivatives}
            , vals_interval_{std::move(vals_interval)}
            , vals_point_{std::move(vals_point)} { }

            auto operator ()(dof_index dof, point_index q, point normal) const noexcept -> facet_value<value_type> {
                const auto [ix, iy] = space_->index_on_facet(dof, facet_);
                const auto [nx, ny] = normal;

                if (facet_.dir == orientation::horizontal) {
                    const auto avg = eval_tensor_basis(
                            [&,ix=ix](int der) { return vals_interval_(q, ix, der); },
                            [&,iy=iy](int der) { return vals_point_.average(iy, der); });

                    const auto jump = eval_tensor_basis(
                            [&,ix=ix]      (int der) { return vals_interval_(q, ix, der); },
                            [&,iy=iy,ny=ny](int der) { return vals_point_.jump(iy, der, ny); });

                    return {avg, jump};
                } else {
                    const auto avg = eval_tensor_basis(
                            [&,ix=ix](int der) { return vals_point_.average(ix, der); },
                            [&,iy=iy](int der) { return vals_interval_(q, iy, der); });

                    const auto jump = eval_tensor_basis(
                            [&,ix=ix,nx=nx](int der) { return vals_point_.jump(ix, der, nx); },
                            [&,iy=iy]      (int der) { return vals_interval_(q, iy, der); });

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
            const auto [dx, dy]      = dof;

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
            const auto order = standard_ordering<2>{{nx, ny}};
            return order.linear_index(ix, iy);
        }
    };


    class bspline_function {
    private:
        const space*              space_;
        const double*             coefficients_;
        mutable bspline::eval_ctx ctx_x_;
        mutable bspline::eval_ctx ctx_y_;
        mutable std::mutex        ctx_lock_;

    public:
        using point  = space::point;

        bspline_function(const space* space, const double* coefficients)
        : space_{space}
        , coefficients_{coefficients}
        , ctx_x_{space->space_x().degree()}
        , ctx_y_{space->space_y().degree()}
        { }

        auto operator ()(point p) const noexcept -> double {
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

        auto operator ()(point p, const regular_mesh::edge_data& edge) const noexcept -> facet_value<double> {
            // TODO: implement this correctly
            if (edge.type == facet_type::interior) {
                const auto [px, py] = p;
                const auto [nx, ny] = edge.normal;
                const auto eps = 1e-16;
                const auto before = point{px - eps * nx, py - eps * ny};
                const auto after  = point{px + eps * nx, py + eps * ny};

                const auto v0 = (*this)(before);
                const auto v1 = (*this)(after);
                return {0.5 * (v0 + v1), v0 - v1};
            } else {
                const auto v = (*this)(p);
                return {v, v};
            }
        }
    };


    auto evenly_spaced(double a, double b, int elems) -> partition {
        assert(elems > 0 && "Invalid number of partition elements");
        auto xs = partition(elems + 1);

        for (int i = 0; i <= elems; ++ i) {
            xs[i] = lerp(i, elems, a, b);
        }
        return xs;
    }

    auto make_bspline_basis(const partition& points, int p, int c) -> bspline::basis {
        const auto n = as_signed(points.size());
        const auto r = p - c;
        auto knot = bspline::knot_vector{};
        knot.reserve((p + 1) + (n - 1) * r + (p + 1));

        auto append = [&knot](int k, double x) {
            std::fill_n(back_inserter(knot), k, x);
        };

        append(p + 1, points[0]);
        for (int i = 1; i < n - 1; ++ i) {
            append(r, points[i]);
        }
        append(p + 1, points[n - 1]);

        return bspline::basis{std::move(knot), p};
    }
}


template <typename Space, typename Quad, typename Out, typename Form>
auto assemble(const Space& space, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto e : mesh.elements()) {
        const auto points = quad.coordinates(e);
        const auto eval   = space.dof_evaluator(e, points, 1);

        const auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        using dof_idx    = typename Space::dof_index;
        using point_idx  = typename decltype(points)::point_index;
        using value_type = decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>()));

        auto basis_vals = std::vector<value_type>(n);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs(e)) {
                const auto iloc = space.local_index(i, e);
                basis_vals[iloc] = eval(i, q);
            }
            for (auto i : space.dofs(e)) {
                for (auto j : space.dofs(e)) {
                    const auto iloc = space.local_index(i, e);
                    const auto jloc = space.local_index(j, e);
                    const auto& u = basis_vals[iloc];
                    const auto& v = basis_vals[jloc];
                    M(jloc, iloc) += form(u, v, x) * w;
                }
            }
        }

        for (auto i : space.dofs(e)) {
            const auto iloc = space.local_index(i, e);
            const auto I    = space.global_index(i);
            for (auto j : space.dofs(e)) {
                const auto jloc = space.local_index(j, e);
                const auto J    = space.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Trial, typename Test, typename Quad, typename Out, typename Form>
auto assemble(const Trial& trial, const Test& test, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = test.mesh();

    for (auto e : mesh.elements()) {
        const auto points     = quad.coordinates(e);
        const auto eval_trial = trial.dof_evaluator(e, points, 1);
        const auto eval_test  = test.dof_evaluator(e, points, 1);

        const auto n_trial = trial.dof_count(e);
        const auto n_test  = test.dof_count(e);
        auto M = ads::lin::tensor<double, 2>{{n_test, n_trial}};

        using point_idx     = typename decltype(points)::point_index;
        using trial_dof_idx = typename Trial::dof_index;
        using test_dof_idx  = typename Test::dof_index;
        using trial_value   = decltype(eval_trial(std::declval<trial_dof_idx>(), std::declval<point_idx>()));
        using test_value    = decltype(eval_test(std::declval<test_dof_idx>(), std::declval<point_idx>()));

        auto test_vals  = std::vector<test_value>(n_test);
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
            for (auto i : trial.dofs(e)) {
                for (auto j : test.dofs(e)) {
                    const auto iloc = trial.local_index(i, e);
                    const auto jloc = test.local_index(j, e);
                    const auto& u = trial_vals[iloc];
                    const auto& v = test_vals[jloc];
                    M(jloc, iloc) += form(u, v, x) * w;
                }
            }
        }

        for (auto i : trial.dofs(e)) {
            const auto iloc = trial.local_index(i, e);
            const auto I    = trial.global_index(i);
            for (auto j : test.dofs(e)) {
                const auto jloc = test.local_index(j, e);
                const auto J    = test.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Facets, typename Space, typename Quad, typename Out, typename Form>
auto assemble_facets(const Facets& facets, const Space& space, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto f : facets) {
        const auto facet  = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval   = space.dof_evaluator(f, points, 1);

        const auto n = space.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        using dof_idx    = typename Space::dof_index;
        using point_idx  = typename decltype(points)::point_index;
        using value_type = decltype(eval(std::declval<dof_idx>(), std::declval<point_idx>(), facet.normal));

        auto basis_vals = std::vector<value_type>(n);

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs_on_facet(f)) {
                const auto iloc = space.facet_local_index(i, f);
                basis_vals[iloc] = eval(i, q, facet.normal);
            }
            for (auto i : space.dofs_on_facet(f)) {
                for (auto j : space.dofs_on_facet(f)) {
                    const auto iloc = space.facet_local_index(i, f);
                    const auto jloc = space.facet_local_index(j, f);
                    const auto& u = basis_vals[iloc];
                    const auto& v = basis_vals[jloc];
                    M(jloc, iloc) += form(u, v, x, facet) * w;
                }
            }
        }

        for (auto i : space.dofs_on_facet(f)) {
            const auto iloc = space.facet_local_index(i, f);
            const auto I    = space.global_index(i);
            for (auto j : space.dofs_on_facet(f)) {
                const auto jloc = space.facet_local_index(j, f);
                const auto J    = space.global_index(j);
                out(J, I, M(jloc, iloc));
            }
        }
    }
}

template <typename Facets, typename Trial, typename Test, typename Quad, typename Out, typename Form>
auto assemble_facets(const Facets& facets, const Trial& trial, const Test& test, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = test.mesh();

    for (auto f : facets) {
        const auto facet      = mesh.facet(f);
        const auto points     = quad.coordinates(f);
        const auto eval_trial = trial.dof_evaluator(f, points, 1);
        const auto eval_test  = test.dof_evaluator(f, points, 1);

        const auto n_trial = trial.facet_dof_count(f);
        const auto n_test  = test.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 2>{{n_test, n_trial}};

        using point_idx     = typename decltype(points)::point_index;
        using trial_dof_idx = typename Trial::dof_index;
        using test_dof_idx  = typename Test::dof_index;
        using trial_value   = decltype(eval_trial(std::declval<trial_dof_idx>(),
                                                  std::declval<point_idx>(), facet.normal));
        using test_value    = decltype(eval_test(std::declval<test_dof_idx>(),
                                                 std::declval<point_idx>(), facet.normal));

        auto test_vals  = std::vector<test_value>(n_test);
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
            for (auto i : trial.dofs_on_facet(f)) {
                for (auto j : test.dofs_on_facet(f)) {
                    const auto iloc = trial.facet_local_index(i, f);
                    const auto jloc = test.facet_local_index(j, f);
                    const auto& u = trial_vals[iloc];
                    const auto& v = test_vals[jloc];
                    M(jloc, iloc) += form(u, v, x, facet) * w;
                }
            }
        }

        for (auto i : trial.dofs_on_facet(f)) {
            const auto iloc = trial.facet_local_index(i, f);
            const auto I    = trial.global_index(i);
            for (auto j : test.dofs_on_facet(f)) {
                const auto jloc = test.facet_local_index(j, f);
                const auto J    = test.global_index(j);
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
        const auto eval   = space.dof_evaluator(e, points, 1);

        const auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 1>{{n}};

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto j : space.dofs(e)) {
                const auto v    = eval(j, q);
                const auto jloc = space.local_index(j, e);
                M(jloc) += form(v, x) * w;
            }
        }

        for (auto j : space.dofs(e)) {
            const auto jloc = space.local_index(j, e);
            const auto J    = space.global_index(j);
            out(J, M(jloc));
        }
    }
}

template <typename Facets, typename Space, typename Quad, typename Out, typename Form>
auto assemble_rhs(const Facets& facets, const Space& space, const Quad& quad, Out out, Form&& form) -> void {
    const auto& mesh = space.mesh();

    for (auto f : facets) {
        const auto facet  = mesh.facet(f);
        const auto points = quad.coordinates(f);
        const auto eval   = space.dof_evaluator(f, points, 1);

        const auto n = space.facet_dof_count(f);
        auto M = ads::lin::tensor<double, 1>{{n}};

        for (auto q : points.indices()) {
            const auto [X, w] = points.data(q);
            for (auto j : space.dofs_on_facet(f)) {
                const auto v    = avg(eval(j, q, facet.normal));
                const auto jloc = space.facet_local_index(j, f);
                M(jloc) += form(v, X, facet) * w;
            }
        }
        for (auto j : space.dofs_on_facet(f)) {
            const auto jloc = space.facet_local_index(j, f);
            const auto J    = space.global_index(j);
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
    const auto g     = [&f,&norm](auto x) { return norm(f(x)); };
    const auto value = integrate(mesh, quad, g);
    return std::sqrt(value);
}

template <typename Mesh, typename Quad, typename Norm, typename Function, typename Exact>
auto error(const Mesh& mesh, const Quad& quad, Norm&& norm, Function&& f, Exact&& exact) -> double {
    const auto difference = [&f,&exact](auto x) { return f(x) - exact(x); };
    return ::norm(mesh, quad, std::forward<Norm>(norm), difference);
}

struct L2 {

    constexpr auto operator ()(double v) const noexcept -> double {
        return v * v;
    }

    constexpr auto operator ()(const ads::value_type& v) const noexcept -> double {
        return (*this)(v.val);
    }
};

namespace detail {

    template <typename F>
    struct validity_checker {
        template <typename... Ts>
        constexpr auto operator ()(Ts&&...) const {
            return std::is_invocable<F, Ts...>{};
        }
    };

    template <typename F>
    constexpr auto is_valid(F) {
        return validity_checker<F>{};
    }

    template <typename T>
    inline constexpr bool has_val = detail::is_valid([](auto&& x) -> decltype(x.val){})(T{});
}

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
    using required_data      = derivatives_of_order<1>;
    using computation_method = element_integral;

    auto operator ()(const ads::value_type& v) const noexcept -> double {
        using namespace ads;
        return dot(grad(v), grad(v));
    }
};

struct H1 {
    using required_data      = derivatives_of_order<1>;
    using computation_method = element_integral;

    auto operator ()(const ads::value_type& v) const noexcept -> double {
        using namespace ads;
        return v.val * v.val + dot(grad(v), grad(v));
    }
};

struct L2_2 {
    using required_data      = derivatives_of_order<0>;
    using computation_method = element_integral;

    auto operator ()(double v) const noexcept -> double {
        return v * v;
    }
};


template <typename Function, typename Arg, typename InputTag>
auto eval(Function&& f, const Arg& x, InputTag) {
    return f(x);
}

template <
    typename Function,
    typename Arg,
    typename = std::enable_if_t<
        !std::is_convertible_v<
            std::invoke_result_t<Function, Arg>,
            double
        > && detail::has_val<std::invoke_result_t<Function, Arg>>
    >
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
        static_assert(std::is_invocable_v<Norm, decltype(value)>, "Invalid value type for the specified norm");
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
        const auto facet  = mesh.facet(f);
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
    auto buff = fmt::memory_buffer{};

    for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            fmt::format_to(buff, "{} {}", x, y);
            (fmt::format_to(buff, " {:.7}", funs({x, y})), ...);
            fmt::format_to(buff, "\n");
        }
    }

    auto out = std::ofstream{path};
    out << to_string(buff);
}


constexpr double pi = M_PI;
using std::sin;
using std::cos;

/////////////

template <typename Concrete>
class poisson_base {
private:
    auto self() const -> const Concrete* {
        return static_cast<const Concrete*>(this);
    }

public:
    auto u() const noexcept {
        return [this](auto x) { return self()->u(x); };
    }

    auto f() const noexcept {
        return [this](auto x) { return self()->f(x); };
    }

    auto g() const noexcept {
        return [this](auto x) { return self()->g(x); };
    }
};


class poisson_type1 : private poisson_base<poisson_type1> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        const auto v = sin(pi * x) * sin(pi * y);
        return v * v;
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        const auto sx = sin(pi * x);
        const auto c2x = cos(2 * pi * x);
        const auto sy = sin(pi * y);
        const auto c2y = cos(2 * pi * y);
        return -2 * pi * pi * (c2x * sy * sy + c2y * sx * sx);
    }

    auto g(point /*X*/) const noexcept -> double {
        return 0.0;
    }
};


class poisson_type2 : private poisson_base<poisson_type2> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        return 1 + sin(pi * x) * sin(pi * y);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        return 2 * pi * pi * sin(pi * x) * sin(pi * y);
    }

    auto g(point /*X*/) const noexcept -> double {
        return 1.0;
    }
};


class poisson_type3 : private poisson_base<poisson_type3> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::point_t;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y] = X;
        return x*x + 0.5 * y*y + sin(pi * x) * sin(pi * y);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y] = X;
        return -3 + 2 * pi * pi * sin(pi * x) * sin(pi * y);
    }

    auto g(point X) const noexcept -> double {
        const auto [x, y] = X;
        return x*x + 0.5 * y*y;
    }
};

/////////////

void DG_poisson();
void DG_stokes();
void DGiGRM_stokes();

void DG_poisson_3D();
void DG_stokes_3D();
void DGiGRM_stokes_3D();


int main() {
    // DG_poisson();
    // DG_stokes();
    // DGiGRM_stokes();

    // DG_poisson_3D();
    // DG_stokes_3D();
    DGiGRM_stokes_3D();
}

void DG_poisson() {
    auto elems = 128;
    auto p     = 3;
    auto c     = -1;
    // auto eta   = 1000000.0;  // good for Nitsche BC
    auto eta   = 10.0;

    auto poisson = poisson_type1{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);

    auto mesh  = ads::regular_mesh{xs, ys};
    auto space = ads::space{&mesh, bx, by};
    auto quad  = ads::quadrature{&mesh, p + 1};

    auto n = space.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F       = std::vector<double>(n);
    auto problem = mumps::problem{F.data(), n};
    auto solver  = mumps::solver{};

    auto out = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int J, double val) { F.data()[J] += val; };
    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(space, quad, out, [](auto u, auto v, auto /*x*/) {
        return dot(grad(u), grad(v));
    });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), space, quad, out, [eta](auto u, auto v, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(v)), n) * jump(u).val
               - dot(grad(avg(u)), n) * jump(v).val
               + eta/h * jump(u).val * jump(v).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(space, quad, rhs, [&poisson](auto v, auto x) {
        return v.val * poisson.f(x);
    });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), space, quad, rhs, [eta,&poisson](auto v, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = poisson.g(x);
        return - dot(grad(v), n) * g
               + eta/h * g * v.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    fmt::print("Computing error\n");

    auto u = ads::bspline_function(&space, F.data());

    auto t_before_err = std::chrono::steady_clock::now();
    auto err = error(mesh, quad, L2{}, u, poisson.u());
    auto t_after_err = std::chrono::steady_clock::now();

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", u);
    auto t_after_output = std::chrono::steady_clock::now();

    fmt::print("error = {:.6}\n", err);

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}


template <typename Concrete>
class stokes_base {
private:
    auto self() const -> const Concrete* {
        return static_cast<const Concrete*>(this);
    }

public:
    auto vx() const noexcept {
        return [this](auto x) { return self()->vx(x); };
    }

    auto vy() const noexcept {
        return [this](auto x) { return self()->vy(x); };
    }

    auto vz() const noexcept {
        return [this](auto x) { return self()->vz(x); };
    }

    auto p(double avg = 0) const noexcept {
        return [this,avg](auto x) { return avg + self()->p(x); };
    }

    auto fx() const noexcept {
        return [this](auto x) { return self()->fx(x); };
    }

    auto fy() const noexcept {
        return [this](auto x) { return self()->fy(x); };
    }

    auto fz() const noexcept {
        return [this](auto x) { return self()->fz(x); };
    }
};


class stokes_polynomial : private stokes_base<stokes_polynomial> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p(point p) const noexcept -> double {
        const auto [x, y] = p;
        return x * (1 - x) - 1./6;
    }

    auto vx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return x*x * (1 - x) * (1 - x) * (2 * y - 6 * y*y + 4 * y*y*y);
    }

    auto vy(point p) const noexcept -> double {
        const auto [x, y] = p;
        return -y*y * (1 - y) * (1 - y) * (2 * x - 6 * x*x + 4 * x*x*x);
    }

    auto fx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return (12 - 24 * y) * x*x*x*x +
               (-24 + 48 * y) * x*x*x +
               (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
               (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
               1 - 4 * y + 12 * y*y - 8 * y*y*y;
    }

    auto fy(point p) const noexcept -> double {
        const auto [x, y] = p;
        return (8 - 48 * y + 48 * y*y) * x*x*x
             + (-12 + 72 * y - 72 * y*y) * x*x
             + (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x
             - 12 * y*y + 24 * y*y*y - 12 * y*y*y*y;
    }
};


class stokes_nonpoly : private stokes_base<stokes_nonpoly> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto xx = x * x;
        const auto yy = y * y;
        const auto ex = exp(x);
        const auto e = exp(1);

        return -424 + 156*e + (yy - y) * (-456 + ex * (456 + xx * (228 - 5 * (yy - y)) + 2 * x * (-228 + (yy - y)) + 2 * pow(x,3) * (-36 + (yy - y)) + pow(x,4) * (12 + (yy - y))));
    };

    auto vx(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto ex = exp(x);

        return 2*ex * pow(-1 + x, 2) * x * x * (y * y - y) * (-1 + 2*y);
    }

    auto vy(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto ex = exp(x);

        return -ex * (-1 + x) * x* (-2 + x * (3 + x)) * pow(-1 + y, 2) * y * y;
    }

    auto fx(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto xx = x * x;
        const auto yy = y * y;
        const auto ex = exp(x);

        const auto px = ex * (y - 1) * y * (pow(x,4) * (yy - y + 12) + 6 * pow(x,3) * (yy - y - 4) + xx * (yy - y + 12) - 8 * x * (y - 1) * y + 2 * (y - 1) * y);
        const auto Lux = 2 * ex * (pow(x,4) * (2 * pow(y,3) - 3 * yy + 13 * y - 6) + 6 * pow(x,3) * (2 * pow(y,3) - 3 * yy - 3 * y + 2) +
                                xx * (2 * pow(y,3) - 3 * yy + 13 * y - 6) - 8 * x * y * (2 * yy - 3 * y + 1) + 2 * y * (2 * yy - 3 * y + 1));

        return -Lux + px;
    }

    auto fy(point p) const noexcept -> double {
        using std::exp;
        using std::pow;
        const auto [x, y] = p;
        const auto xx = x * x, yy = y * y;
        const auto ex = exp(x);

        const auto py = 2 * (2 * y - 1) * (ex * (pow(x,4) * (yy - y + 6) + 2 * pow(x,3) * (yy - y - 18) + xx * (-5 * yy + 5 * y + 114) + 2 * x * (yy - y - 114) + 228) - 228);
        const auto Luy = -ex * (pow(x,4) * (pow(y,4) - 2 * pow(y,3) + 13 * yy - 12 * y + 2) + 2 * pow(x,3) * (5 * pow(y,4) - 10 * pow(y,3) + 17 * yy - 12 * y + 2) +
                            xx * (19 * pow(y,4) - 38 * pow(y,3) - 41 * yy + 60 * y - 10) + x * (-6 * pow(y,4) + 12 * pow(y,3) + 18 * yy - 24 * y + 4) - 6 * pow(y - 1, 2) * yy);

        return -Luy + py;
    }
};

class stokes_cavity : private stokes_base<stokes_cavity> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::point_t;
    using base::p, base::vx, base::vy, base::fx, base::fy;

    auto p([[maybe_unused]] point p) const noexcept -> double {
        return 0.0;
    }

    auto vx(point p) const noexcept -> double {
        const auto [x, y] = p;
        return std::abs(y - 1) < 1e-5;
    }

    auto vy([[maybe_unused]] point p) const noexcept -> double {
        return 0.0;
    }

    auto fx([[maybe_unused]] point p) const noexcept -> double {
        return 0.0;
    }

    auto fy([[maybe_unused]] point p) const noexcept -> double {
        return 0.0;
    }
};

class space_factory {
private:
    ads::global_dof offset_ = 0;

public:
    template <typename Space, typename... Args>
    auto next(Args&&... args) -> Space {
        auto space = Space{std::forward<Args>(args)..., offset_};
        offset_ += space.dof_count();
        return space;
    }
};


void DG_stokes() {
    auto elems = 32;
    auto p     = 3;
    auto c     = -1;
    auto eta   = 10.0 * (p + 1) * (p + 2);

    // auto stokes = stokes_cavity{};
    // auto stokes = stokes_polynomial{};
    auto stokes = stokes_nonpoly{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);

    auto mesh  = ads::regular_mesh{xs, ys};
    auto quad  = ads::quadrature{&mesh, std::max(p + 1, 2)};

    auto spaces = space_factory{};

    auto Vx = spaces.next<ads::space>(&mesh, bx, by);
    auto Vy = spaces.next<ads::space>(&mesh, bx, by);
    auto P  = spaces.next<ads::space>(&mesh, bx, by);

    auto n = Vx.dof_count() + Vy.dof_count() + P.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F       = std::vector<double>(n);
    auto problem = mumps::problem{F.data(), n};
    auto solver  = mumps::solver{};

    auto M = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(Vx,    quad, M, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy,    quad, M, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(P, Vx, quad, M, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P, Vy, quad, M, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(Vx, P, quad, M, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy, P, quad, M, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), Vx, quad, M, [eta](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, quad, M, [eta](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), P, Vx, quad, M, [](auto p, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{jump(vx).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vy, quad, M, [](auto p, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{0, jump(vy).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, P, quad, M, [](auto ux, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{jump(ux).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, P, quad, M, [](auto uy, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{0, jump(uy).val};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.interior_facets(), P, quad, M, [](auto p, auto q, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return h * jump(p).val * jump(q).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Vx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Vy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), Vx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function(&Vx, F.data());
    auto sol_vy = ads::bspline_function(&Vy, F.data());
    auto sol_p  = ads::bspline_function(&P,  F.data());

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean   = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_p  = error(mesh, quad, L2{},  sol_p,  stokes.p(mean));
    auto err    = sum_norms(err_vx, err_vy, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    const auto vx_J = integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        const auto d = jump(sol_vx(x, edge));
        return 1/h * d * d;
    });
    const auto vy_J = integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        const auto d = jump(sol_vy(x, edge));
        return 1/h * d * d;
    });
    const auto r_P = integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        auto d = jump(sol_p(x, edge));
        return h * d * d;
    });
    fmt::print("vx seminorm = {:.6}\n", std::sqrt(vx_J));
    fmt::print("vy seminorm = {:.6}\n", std::sqrt(vy_J));
    fmt::print("p  seminorm = {:.6}\n", std::sqrt(r_P));

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", sol_vx, sol_vy, sol_p);
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}

void DGiGRM_stokes() {
    auto elems = 64;
    auto p     = 4;
    auto eta   = 10.0 * (p + 1) * (p + 2);

    // auto stokes = stokes_cavity{};
    // auto stokes = stokes_polynomial{};
    auto stokes = stokes_nonpoly{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto mesh  = ads::regular_mesh{xs, ys};
    auto quad  = ads::quadrature{&mesh, 5};

    // Test
    auto p_test =  4;
    auto c_test = -1;

    auto Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto By = ads::make_bspline_basis(xs, p_test, c_test);

    auto tests = space_factory{};

    auto Wx = tests.next<ads::space>(&mesh, Bx, By);
    auto Wy = tests.next<ads::space>(&mesh, Bx, By);
    auto Q  = tests.next<ads::space>(&mesh, Bx, By);

    auto N = Wx.dof_count() + Wy.dof_count() + Q.dof_count();
    fmt::print("Test  DoFs: {:10L}\n", N);

    // Trial
    auto p_trial = 4;
    auto c_trial = 0;   // >= 0

    auto bx = ads::make_bspline_basis(xs, p_trial, c_trial);
    auto by = ads::make_bspline_basis(ys, p_trial, c_trial);

    auto trials = space_factory{};

    auto Vx = trials.next<ads::space>(&mesh, bx, by);
    auto Vy = trials.next<ads::space>(&mesh, bx, by);
    auto P  = trials.next<ads::space>(&mesh, bx, by);

    auto n = Vx.dof_count() + Vy.dof_count() + P.dof_count();
    fmt::print("Trial DoFs: {:10L}\n", n);
    fmt::print("Total:      {:10L}\n", N + n);

    auto F       = std::vector<double>(N + n);
    auto problem = mumps::problem{F.data(), F.size()};
    auto solver  = mumps::solver{};

    auto G = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto B = [&problem,N](int row, int col, double val) {
        if (val != 0) {
            problem.add(    row + 1, N + col + 1, val);
            problem.add(N + col + 1,     row + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(Wx, quad, G, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Wy, quad, G, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Q,  quad, G, [](auto p,  auto q,  auto /*x*/) { return p.val * q.val;           });

    assemble(Vx, Wx, quad, B, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy, Wy, quad, B, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(P,  Wx, quad, B, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P,  Wy, quad, B, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(Vx,  Q, quad, B, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy,  Q, quad, B, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), Wx, quad, G, [](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return 1/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Wy, quad, G, [](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return 1/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.interior_facets(), Q, quad, G, [](auto p, auto q, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return h * jump(p).val * jump(q).val;
    });

    assemble_facets(mesh.facets(), Vx, Wx, quad, B, [eta](auto ux, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, Wy, quad, B, [eta](auto uy, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), P, Wx, quad, B, [](auto p, auto vx, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{jump(vx).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wy, quad, B, [](auto p, auto vy, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  v = ads::point_t{0, jump(vy).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.boundary_facets(), Vx, Q, quad, B, [](auto ux, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{jump(ux).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.boundary_facets(), Vy, Q, quad, B, [](auto uy, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{0, jump(uy).val};
        return - dot(u, n) * avg(q).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Wx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Wy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), Wx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& edge) {
        const auto& n = edge.normal;
        const auto  h = length(edge.span);
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function(&Vx, F.data() + N);
    auto sol_vy = ads::bspline_function(&Vy, F.data() + N);
    auto sol_p  = ads::bspline_function(&P,  F.data() + N);

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean   = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_p  = error(mesh, quad, L2{}, sol_p,  stokes.p(mean));
    auto err    = sum_norms(err_vx, err_vy, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto r_vx = ads::bspline_function(&Wx, F.data());
    auto r_vy = ads::bspline_function(&Wy, F.data());
    auto r_p  = ads::bspline_function(&Q,  F.data());

    auto norm_r_vx = norm(mesh, quad, L2{}, r_vx);
    auto J_r_vx = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        auto d = jump(r_vx(x, edge));
        return 1/h * d * d;
    }));
    auto norm_r_vy = norm(mesh, quad, L2{}, r_vy);
    auto J_r_vy = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        auto d = jump(r_vy(x, edge));
        return 1/h * d * d;
    }));
    auto norm_r_p = norm(mesh, quad, L2{}, r_p);
    auto q_r_p = std::sqrt(integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& edge) {
        const auto  h = length(edge.span);
        auto d = jump(r_p(x, edge));
        return h * d * d;
    }));
    auto res_norm = sum_norms(norm_r_vx, J_r_vx, norm_r_vy, J_r_vy, norm_r_p, q_r_p);
    fmt::print("||r_vx|| = {:.6}\n", norm_r_vx);
    fmt::print("||r_vy|| = {:.6}\n", norm_r_vy);
    fmt::print("||r_p||  = {:.6}\n", norm_r_p);
    fmt::print("|r_vx|   = {:.6}\n", J_r_vx);
    fmt::print("|r_vy|   = {:.6}\n", J_r_vy);
    fmt::print("|r_p|    = {:.6}\n", q_r_p);
    fmt::print("res norm = {:.6}\n", res_norm);

    auto val = alt_norm(mesh, quad, L2_2{}, [](auto X) {
        auto [x, y] = X;
        return ads::value_type{x*x - y*y, 2*x, -2*y};
        // return 1;
        // return ads::facet_value<double>{0, 1};
    });
    fmt::print("H1 norm : {}\n", val);

    auto t_before_output = std::chrono::steady_clock::now();
    save_to_file("result.data", sol_vx, sol_vy, sol_p, r_vx, r_vy, r_p);
    auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}

namespace ads {


    struct index_types3 {
        using index          = std::tuple<int, int, int>;
        using index_iterator = util::iter_product3<simple_index_iterator, index>;
        using index_range    = boost::iterator_range<index_iterator>;
    };

    enum class orientation3 {
        dir_x,
        dir_y,
        dir_z
    };

    class regular_mesh3 {
    private:
        interval_mesh mesh_x_;
        interval_mesh mesh_y_;
        interval_mesh mesh_z_;

    public:
        using point            = std::tuple<double, double, double>;
        using element_index    = index_types3::index;
        using element_iterator = index_types3::index_iterator;
        using element_range    = index_types3::index_range;

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
        , mesh_z_{std::move(zs)}
        { }

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
            interval     span1;
            interval     span2;
            double       position;
            double       diameter;
            orientation3 direction;
            facet_type   type;
            point        normal;
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
                const auto span_y        = mesh_y_.subinterval(iy);
                const auto span_z        = mesh_z_.subinterval(iz);
                const auto normal        = point{nx, 0, 0};
                const auto diam          = diameter(span_y, span_z);
                return {span_y, span_z, x, diam, dir, type, normal};
            } else if (dir == orientation3::dir_y) {
                // 1 - x, 2 - z
                const auto span_x        = mesh_x_.subinterval(ix);
                const auto [y, type, ny] = mesh_y_.facet(iy);
                const auto span_z        = mesh_z_.subinterval(iz);
                const auto normal        = point{0, ny, 0};
                const auto diam          = diameter(span_x, span_z);
                return {span_x, span_z, y, diam, dir, type, normal};
            } else {
                assert(dir == orientation3::dir_z && "Invalid face orientation");
                // 1 - x, 2 - y
                const auto span_x        = mesh_x_.subinterval(ix);
                const auto span_y        = mesh_y_.subinterval(iy);
                const auto [z, type, nz] = mesh_z_.facet(iz);
                const auto normal        = point{0, 0, nz};
                const auto diam          = diameter(span_x, span_y);
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
        using point          = std::tuple<double, double, double>;
        using point_index    = index_types3::index;
        using point_iterator = index_types3::index_iterator;
        using point_range    = index_types3::index_range;

        tensor_quadrature_points3(interval_quadrature_points ptx,
                                  interval_quadrature_points pty,
                                  interval_quadrature_points ptz) noexcept
        : ptx_{std::move(ptx)}
        , pty_{std::move(pty)}
        , ptz_{std::move(ptz)}
        { }

        auto xs() const noexcept -> const std::vector<double>& {
            return ptx_.points();
        }

        auto ys() const noexcept -> const std::vector<double>& {
            return pty_.points();
        }

        auto zs() const noexcept -> const std::vector<double>& {
            return ptz_.points();
        }

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
            point  x;
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
        double                     position_;
        orientation3               direction_;

    public:
        using point          = std::tuple<double, double, double>;
        using point_index    = index_types::index;
        using point_iterator = index_types::index_iterator;
        using point_range    = index_types::index_range;

        face_quadrature_points(interval_quadrature_points points1, interval_quadrature_points points2,
                               double position, orientation3 direction) noexcept
        : points1_{std::move(points1)}
        , points2_{std::move(points2)}
        , position_{position}
        , direction_{direction}
        { }

        auto points1() const noexcept -> const std::vector<double>& {
            return points1_.points();
        }

        auto points2() const noexcept -> const std::vector<double>& {
            return points2_.points();
        }

        auto position() const noexcept -> double {
            return position_;
        }

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
            point  x;
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
        using point          = regular_mesh3::point;
        using element_index  = regular_mesh3::element_index;
        using facet_index    = regular_mesh3::facet_index;
        using point_set      = tensor_quadrature_points3;
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

            return{std::move(pts1), std::move(pts2), face.position, face.direction};
        }

    private:
        // TODO: remove duplication
        auto data_for_interval(interval target) const -> interval_quadrature_points {
            const auto size    = length(target);
            const auto scale   = size / 2;        // Gauss quadrature is defined for [-1, 1]
            const auto weights = quad::gauss::Ws[point_count_];

            return {transform_points(target), weights, scale};
        }

        auto transform_points(interval target) const -> std::vector<double> {
            auto points = std::vector<double>(point_count_);

            for (int i = 0; i < point_count_; ++ i) {
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
        return std::get<0>(a) * std::get<0>(b)
             + std::get<1>(a) * std::get<1>(b)
             + std::get<2>(a) * std::get<2>(b);
    }

    auto grad(const value_type3& v) noexcept -> point3_t {
        return {v.dx, v.dy, v.dz};
    }

    template <typename ValsX, typename ValsY, typename ValsZ>
    inline auto eval_tensor_basis(const ValsX& vals_x, const ValsY& vals_y, const ValsZ& vals_z) noexcept
    -> value_type3 {
        const auto  Bx = vals_x(0);
        const auto dBx = vals_x(1);
        const auto  By = vals_y(0);
        const auto dBy = vals_y(1);
        const auto  Bz = vals_z(0);
        const auto dBz = vals_z(1);

        const auto v  =  Bx *  By *  Bz;
        const auto dx = dBx *  By *  Bz;
        const auto dy =  Bx * dBy *  Bz;
        const auto dz =  Bx *  By * dBz;

        return {v, dx, dy, dz};
    }

    class space3 {
    private:
        regular_mesh3* mesh_;
        bspline_space  space_x_;
        bspline_space  space_y_;
        bspline_space  space_z_;
        global_dof     dof_offset_;

        class evaluator;
        class face_evaluator;

    public:
        using point         = regular_mesh3::point;
        using element_index = regular_mesh3::element_index;
        using facet_index   = regular_mesh3::facet_index;
        using dof_index     = index_types3::index;
        using dof_iterator  = index_types3::index_iterator;
        using dof_range     = index_types3::index_range;

        space3(regular_mesh3* mesh, bspline::basis bx, bspline::basis by, bspline::basis bz,
               global_dof dof_offset = 0) noexcept
        : mesh_{mesh}
        , space_x_{std::move(bx)}
        , space_y_{std::move(by)}
        , space_z_{std::move(bz)}
        , dof_offset_{dof_offset}
        { }

        auto mesh() const noexcept -> const regular_mesh3& {
            return *mesh_;
        }

        auto space_x() const noexcept -> const bspline_space& {
            return space_x_;
        }

        auto space_y() const noexcept -> const bspline_space& {
            return space_y_;
        }

        auto space_z() const noexcept -> const bspline_space& {
            return space_z_;
        }

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

        auto dof_evaluator(element_index e, const tensor_quadrature_points3& points, int ders) const -> evaluator {
            auto data_x = evaluate_basis(points.xs(), space_x_, ders);
            auto data_y = evaluate_basis(points.ys(), space_y_, ders);
            auto data_z = evaluate_basis(points.zs(), space_z_, ders);

            return evaluator{this, e, ders, std::move(data_x), std::move(data_y), std::move(data_z)};
        }

        auto dof_evaluator(facet_index f, const face_quadrature_points& points, int ders) const -> face_evaluator {
            const auto [fx, fy, fz, dir] = f;

            if (dir == orientation3::dir_x) {
                // 1 - y, 2 - z
                auto data_x = evaluate_basis(fx, space_x_, ders);
                auto data_y = evaluate_basis(points.points1(), space_y_, ders);
                auto data_z = evaluate_basis(points.points2(), space_z_, ders);
                return face_evaluator{this, f, ders, std::move(data_y), std::move(data_z), std::move(data_x)};
            } else if (dir == orientation3::dir_y) {
                // 1 - x, 2 - z
                auto data_x = evaluate_basis(points.points1(), space_x_, ders);
                auto data_y = evaluate_basis(fy, space_y_, ders);
                auto data_z = evaluate_basis(points.points2(), space_z_, ders);
                return face_evaluator{this, f, ders, std::move(data_x), std::move(data_z), std::move(data_y)};
            } else {
                assert(dir == orientation3::dir_z && "Invalid face orientation");
                // 1 - x, 2 - y
                auto data_x = evaluate_basis(points.points1(), space_x_, ders);
                auto data_y = evaluate_basis(points.points2(), space_y_, ders);
                auto data_z = evaluate_basis(fz, space_z_, ders);
                return face_evaluator{this, f, ders, std::move(data_x), std::move(data_y), std::move(data_z)};
            }
        }

    private:

        class evaluator {
        private:
            const space3*        space_;
            element_index        element_;
            int                  derivatives_;
            bspline_basis_values vals_x_;
            bspline_basis_values vals_y_;
            bspline_basis_values vals_z_;

        public:
            using point_index = tensor_quadrature_points3::point_index;

            evaluator(const space3* space, element_index element, int derivatives,
                      bspline_basis_values vals_x, bspline_basis_values vals_y, bspline_basis_values vals_z) noexcept
            : space_{space}
            , element_{element}
            , derivatives_{derivatives}
            , vals_x_{std::move(vals_x)}
            , vals_y_{std::move(vals_y)}
            , vals_z_{std::move(vals_z)}
            { }

            auto operator ()(dof_index dof, point_index q) const noexcept -> value_type3 {
                const auto [qx, qy, qz] = q;
                const auto [ix, iy, iz] = space_->index_on_element(dof, element_);

                return eval_tensor_basis([&,qx=qx,ix=ix](int der) { return vals_x_(qx, ix, der); },
                                         [&,qy=qy,iy=iy](int der) { return vals_y_(qy, iy, der); },
                                         [&,qz=qz,iz=iz](int der) { return vals_z_(qz, iz, der); });
            }
        };

        class face_evaluator {
        private:
            const space3* space_;
            facet_index   facet_;
            int           derivatives_;

            // dir x: 1 - y, 2 - z
            // dir y: 1 - x, 2 - z
            // dir z: 1 - x, 2 - y
            bspline_basis_values vals_interval1_;
            bspline_basis_values vals_interval2_;

            bspline_basis_values_on_vertex vals_point_;

        public:
            // using point_index = face_quadrature_points::point_index;
            using point_index = std::tuple<int, int>;

            face_evaluator(const space3* space, facet_index facet, int derivatives,
                           bspline_basis_values vals_interval1, bspline_basis_values vals_interval2,
                           bspline_basis_values_on_vertex vals_point) noexcept
            : space_{space}
            , facet_{facet}
            , derivatives_{derivatives}
            , vals_interval1_{std::move(vals_interval1)}
            , vals_interval2_{std::move(vals_interval2)}
            , vals_point_{std::move(vals_point)}
            { }

            auto operator ()(dof_index dof, point_index q, point normal) const noexcept -> facet_value<value_type3> {
                const auto [ix, iy, iz] = space_->index_on_facet(dof, facet_);
                const auto [nx, ny, nz] = normal;
                const auto [q1, q2] = q;

                if (facet_.dir == orientation3::dir_x) {
                    // 1 - y, 2 - z
                    const auto val_y  = [&,iy=iy,q1=q1](int der) { return vals_interval1_(q1, iy, der); };
                    const auto val_z  = [&,iz=iz,q2=q2](int der) { return vals_interval2_(q2, iz, der); };
                    const auto avg_x  = [&,ix=ix]      (int der) { return vals_point_.average(ix, der); };
                    const auto jump_x = [&,ix=ix,nx=nx](int der) { return vals_point_.jump(ix, der, nx); };

                    const auto avg  = eval_tensor_basis(avg_x,  val_y, val_z);
                    const auto jump = eval_tensor_basis(jump_x, val_y, val_z);

                    return {avg, jump};
                } else if (facet_.dir == orientation3::dir_y) {
                    // dir y: 1 - x, 2 - z
                    const auto val_x  = [&,ix=ix,q1=q1](int der) { return vals_interval1_(q1, ix, der); };
                    const auto val_z  = [&,iz=iz,q2=q2](int der) { return vals_interval2_(q2, iz, der); };
                    const auto avg_y  = [&,iy=iy]      (int der) { return vals_point_.average(iy, der); };
                    const auto jump_y = [&,iy=iy,ny=ny](int der) { return vals_point_.jump(iy, der, ny); };

                    const auto avg  = eval_tensor_basis(val_x, avg_y,  val_z);
                    const auto jump = eval_tensor_basis(val_x, jump_y, val_z);

                    return {avg, jump};
                } else {
                    // dir z: 1 - x, 2 - y
                    const auto val_x  = [&,ix=ix,q1=q1](int der) { return vals_interval1_(q1, ix, der); };
                    const auto val_y  = [&,iy=iy,q2=q2](int der) { return vals_interval2_(q2, iy, der); };
                    const auto avg_z  = [&,iz=iz]      (int der) { return vals_point_.average(iz, der); };
                    const auto jump_z = [&,iz=iz,nz=nz](int der) { return vals_point_.jump(iz, der, nz); };

                    const auto avg  = eval_tensor_basis(val_x, val_y, avg_z);
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
            const auto [dx, dy, dz]      = dof;

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
            const auto order = standard_ordering<3>{{nx, ny, nz}};
            return order.linear_index(ix, iy, iz);
        }
    };


    class bspline_function3 {
    private:
        const space3*             space_;
        const double*             coefficients_;
        mutable bspline::eval_ctx ctx_x_;
        mutable bspline::eval_ctx ctx_y_;
        mutable bspline::eval_ctx ctx_z_;
        mutable std::mutex        ctx_lock_;

    public:
        using point = space3::point;

        bspline_function3(const space3* space, const double* coefficients)
        : space_{space}
        , coefficients_{coefficients}
        , ctx_x_{space->space_x().degree()}
        , ctx_y_{space->space_y().degree()}
        , ctx_z_{space->space_z().degree()}
        { }

        auto operator ()(point p) const noexcept -> double {
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

        auto operator ()(point p, const regular_mesh3::face_data& face) const noexcept -> facet_value<double> {
            // TODO: implement this correctly
            if (face.type == facet_type::interior) {
                const auto [px, py, pz] = p;
                const auto [nx, ny, nz] = face.normal;
                const auto eps = 1e-16;
                const auto before = point{px - eps * nx, py - eps * ny, pz - eps * nz};
                const auto after  = point{px + eps * nx, py + eps * ny, pz + eps * nz};

                const auto v0 = (*this)(before);
                const auto v1 = (*this)(after);
                return {0.5 * (v0 + v1), v0 - v1};
            } else {
                const auto v = (*this)(p);
                return {v, v};
            }
        }
    };

}


class poisson3_type1 : private poisson_base<poisson3_type1> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto g([[maybe_unused]] point X) const noexcept -> double {
        return 0.0;
    }
};

class poisson3_type2 : private poisson_base<poisson3_type2> {
private:
    using base = poisson_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::u, base::f, base::g;

    auto u(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x*x + 0.5 * y*y + 0.3 * z*z + sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto f(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return -3.6 + 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z);
    }

    auto g(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x*x + 0.5 * y*y + 0.3 * z*z;
    }
};

void DG_poisson_3D() {
    auto elems = 8;
    auto p     = 4;
    auto c     = -1;
    auto eta   = 10.0 * (p + 1) * (p + 2);

    auto poisson = poisson3_type2{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);
    auto bz = ads::make_bspline_basis(zs, p, c);

    auto mesh  = ads::regular_mesh3{xs, ys, zs};
    auto quad  = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto space = ads::space3{&mesh, bx, by, bz};

    auto n = space.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F       = std::vector<double>(n);
    auto problem = mumps::problem{F.data(), n};
    auto solver  = mumps::solver{};

    auto out = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int J, double val) { F.data()[J] += val; };
    using ads::dot;
    using ads::grad;

    // auto executor = ads::galois_executor{12};
    // auto executor = ads::sequential_executor{};

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(space, quad, out, [](auto u, auto v, auto /*x*/) {
        return dot(grad(u), grad(v));
    });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), space, quad, out, [eta](auto u, auto v, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(v)), n) * jump(u).val
               - dot(grad(avg(u)), n) * jump(v).val
               + eta/h * jump(u).val * jump(v).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(space, quad, rhs, [&poisson](auto v, auto x) {
        return v.val * poisson.f(x);
    });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), space, quad, rhs, [eta,&poisson](auto v, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = poisson.g(x);
        return - dot(grad(v), n) * g
               + eta/h * g * v.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    fmt::print("Computing error\n");

    auto u = ads::bspline_function3(&space, F.data());

    auto t_before_err = std::chrono::steady_clock::now();
    auto err = error(mesh, quad, L2{}, u, poisson.u());
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("error = {:.6}\n", err);

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    // fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}

class stokes3_type1 : private stokes_base<stokes3_type1> {
private:
    using base = stokes_base;
    friend base;

public:
    using point = ads::regular_mesh3::point;
    using base::p, base::vx, base::vy, base::vz, base::fx, base::fy, base::fz;

    auto p(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x * (1 - x) - 1./6;
    }

    auto vx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return x*x * (1 - x) * (1 - x) * (2 * y - 6 * y*y + 4 * y*y*y);
    }

    auto vy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return -y*y * (1 - y) * (1 - y) * (2 * x - 6 * x*x + 4 * x*x*x);
    }

    auto vz(point /*X*/) const noexcept -> double {
        return 0.0;
    }

    auto fx(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return (12 - 24 * y) * x*x*x*x +
               (-24 + 48 * y) * x*x*x +
               (-48 * y + 72 * y*y - 48 * y*y*y + 12) * x*x +
               (-2 + 24*y - 72 * y*y + 48 * y*y*y) * x +
               1 - 4 * y + 12 * y*y - 8 * y*y*y;
    }

    auto fy(point X) const noexcept -> double {
        const auto [x, y, z] = X;
        return (8 - 48 * y + 48 * y*y) * x*x*x
             + (-12 + 72 * y - 72 * y*y) * x*x
             + (4 - 24 * y + 48 * y*y - 48 * y*y*y + 24 * y*y*y*y) * x
             - 12 * y*y + 24 * y*y*y - 12 * y*y*y*y;
    }

    auto fz(point /*X*/) const noexcept -> double {
        return 0.0;
    }
};

template <typename Vx, typename Vy, typename Vz, typename P>
auto save_to_file3(const std::string& path, Vx&& vx, Vy&& vy, Vz&& vz, P&& pressure) -> void {
    constexpr auto res = 50;
    auto extent = fmt::format("0 {0} 0 {0} 0 {0}", res);

    auto out = fmt::output_file(path);
    out.print("<?xml version=\"1.0\"?>\n");
    out.print("<VTKFile type=\"ImageData\" version=\"0.1\">\n");
    out.print("  <ImageData WholeExtent=\"{}\" origin=\"0 0 0\" spacing=\"1 1 1\">\n", extent);
    out.print("    <Piece Extent=\"{}\">\n", extent);
    out.print("      <PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n", extent);

    out.print("        <DataArray Name=\"Velocity\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    for (auto z : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
                const auto X = ads::point3_t{x, y, z};
                out.print("{:.7} {:.7} {:.7}\n", vx(X), vy(X), vz(X));
            }
        }
    }
    out.print("        </DataArray>\n");

    out.print("        <DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\">\n");
    for (auto z : ads::evenly_spaced(0.0, 1.0, res)) {
        for (auto y : ads::evenly_spaced(0.0, 1.0, res)) {
            for (auto x : ads::evenly_spaced(0.0, 1.0, res)) {
                const auto X = ads::point3_t{x, y, z};
                out.print("{:.7}\n", pressure(X));
            }
        }
    }
    out.print("        </DataArray>\n");

    out.print("      </PointData>\n");
    out.print("    </Piece>\n");
    out.print("  </ImageData>\n");
    out.print("</VTKFile>\n");
}


void DG_stokes_3D() {
    auto elems = 4;
    auto p     = 4;
    auto c     = -1;
    auto eta   = 10.0 * (p + 1) * (p + 2);

    auto stokes = stokes3_type1{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);
    auto bz = ads::make_bspline_basis(zs, p, c);

    auto mesh  = ads::regular_mesh3{xs, ys, zs};
    auto quad  = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    auto spaces = space_factory{};

    auto Vx = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto Vy = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto Vz = spaces.next<ads::space3>(&mesh, bx, by, bz);
    auto P  = spaces.next<ads::space3>(&mesh, bx, by, bz);

    auto n = Vx.dof_count() + Vy.dof_count() + Vz.dof_count() + P.dof_count();
    fmt::print("DoFs: {}\n", n);

    auto F       = std::vector<double>(n);
    auto problem = mumps::problem{F.data(), n};
    auto solver  = mumps::solver{};

    auto M = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(Vx,    quad, M, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy,    quad, M, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Vz,    quad, M, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(P, Vx, quad, M, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P, Vy, quad, M, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(P, Vz, quad, M, [](auto p,  auto vz, auto /*x*/) { return - p.val * vz.dz;         });
    assemble(Vx, P, quad, M, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy, P, quad, M, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    assemble(Vz, P, quad, M, [](auto uz, auto  q, auto /*x*/) { return   uz.dz * q.val;         });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), Vx, quad, M, [eta](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, quad, M, [eta](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Vz, quad, M, [eta](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vz)), n) * jump(uz).val
               - dot(grad(avg(uz)), n) * jump(vz).val
               + eta/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.facets(), P, Vx, quad, M, [](auto p, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{jump(vx).val, 0, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vy, quad, M, [](auto p, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, jump(vy).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Vz, quad, M, [](auto p, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, 0, jump(vz).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, P, quad, M, [](auto ux, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{jump(ux).val, 0, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, P, quad, M, [](auto uy, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, jump(uy).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vz, P, quad, M, [](auto uz, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, 0, jump(uz).val};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.interior_facets(), P, quad, M, [](auto p, auto q, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return h * jump(p).val * jump(q).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Vx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Vy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    assemble_rhs(Vz, quad, rhs, [&stokes](auto vz, auto x) { return vz.val * stokes.fz(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), Vx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    assemble_rhs(mesh.boundary_facets(), Vz, quad, rhs, [eta,&stokes](auto vz, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vz(x);
        return - dot(grad(vz), n) * g
               + eta/h * g * vz.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function3(&Vx, F.data());
    auto sol_vy = ads::bspline_function3(&Vy, F.data());
    auto sol_vz = ads::bspline_function3(&Vz, F.data());
    auto sol_p  = ads::bspline_function3(&P,  F.data());

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean   = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_vz = error(mesh, quad, L2{}, sol_vz, stokes.vz());
    auto err_p  = error(mesh, quad, L2{},  sol_p,  stokes.p(mean));
    auto err    = sum_norms(err_vx, err_vy, err_vz, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("vz error = {:.6}\n", err_vz);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto vf = integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto  h = face.diameter;
        auto d = jump(sol_p(x, face));
        return h * d * d;
    });
    fmt::print("Pressure jump seminorm = {:.6}\n", std::sqrt(vf));

    // auto t_before_output = std::chrono::steady_clock::now();
    // save_to_file("result.data", sol_vx, sol_vy, sol_p);
    // auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    // fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}


void DGiGRM_stokes_3D() {
    auto elems = 4;
    auto p     = 4;
    auto eta   = 10.0 * (p + 1) * (p + 2);

    auto stokes = stokes3_type1{};

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);
    auto zs = ads::evenly_spaced(0.0, 1.0, elems);

    auto mesh  = ads::regular_mesh3{xs, ys, zs};
    auto quad  = ads::quadrature3{&mesh, std::max(p + 1, 2)};

    // Test
    auto p_test =  4;
    auto c_test = -1;

    auto Bx = ads::make_bspline_basis(xs, p_test, c_test);
    auto By = ads::make_bspline_basis(ys, p_test, c_test);
    auto Bz = ads::make_bspline_basis(zs, p_test, c_test);

    auto tests = space_factory{};

    auto Wx = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Wy = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Wz = tests.next<ads::space3>(&mesh, Bx, By, Bz);
    auto Q  = tests.next<ads::space3>(&mesh, Bx, By, Bz);

    auto N = Wx.dof_count() + Wy.dof_count() + Wz.dof_count() + Q.dof_count();
    fmt::print("Test  DoFs: {:10L}\n", N);

    // Trial
    auto p_trial = 4;
    auto c_trial = 0;   // >= 0

    auto bx = ads::make_bspline_basis(xs, p_trial, c_trial);
    auto by = ads::make_bspline_basis(ys, p_trial, c_trial);
    auto bz = ads::make_bspline_basis(zs, p_trial, c_trial);

    auto trials = space_factory{};

    auto Vx = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto Vy = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto Vz = trials.next<ads::space3>(&mesh, bx, by, bz);
    auto P  = trials.next<ads::space3>(&mesh, bx, by, bz);

    auto n = Vx.dof_count() + Vy.dof_count() + Vz.dof_count() + P.dof_count();
    fmt::print("Trial DoFs: {:10L}\n", n);
    fmt::print("Total:      {:10L}\n", N + n);

    auto F       = std::vector<double>(N + n);
    auto problem = mumps::problem{F.data(), N + n};
    auto solver  = mumps::solver{};

    auto G = [&problem](int row, int col, double val) {
        if (val != 0) {
            problem.add(row + 1, col + 1, val);
        }
    };
    auto B = [&problem,N](int row, int col, double val) {
        if (val != 0) {
            problem.add(    row + 1, N + col + 1, val);
            problem.add(N + col + 1,     row + 1, val);
        }
    };
    auto rhs = [&F](int row, double val) { F[row] += val; };

    using ads::dot;
    using ads::grad;

    auto t_before_matrix = std::chrono::steady_clock::now();
    assemble(Wx, quad, G, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Wy, quad, G, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Wz, quad, G, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(Q,  quad, G, [](auto p,  auto q,  auto /*x*/) { return p.val * q.val;           });

    assemble(Vx, Wx, quad, B, [](auto ux, auto vx, auto /*x*/) { return dot(grad(ux), grad(vx)); });
    assemble(Vy, Wy, quad, B, [](auto uy, auto vy, auto /*x*/) { return dot(grad(uy), grad(vy)); });
    assemble(Vz, Wz, quad, B, [](auto uz, auto vz, auto /*x*/) { return dot(grad(uz), grad(vz)); });
    assemble(P,  Wx, quad, B, [](auto p,  auto vx, auto /*x*/) { return - p.val * vx.dx;         });
    assemble(P,  Wy, quad, B, [](auto p,  auto vy, auto /*x*/) { return - p.val * vy.dy;         });
    assemble(P,  Wz, quad, B, [](auto p,  auto vz, auto /*x*/) { return - p.val * vz.dz;         });
    assemble(Vx,  Q, quad, B, [](auto ux, auto  q, auto /*x*/) { return   ux.dx * q.val;         });
    assemble(Vy,  Q, quad, B, [](auto uy, auto  q, auto /*x*/) { return   uy.dy * q.val;         });
    assemble(Vz,  Q, quad, B, [](auto uz, auto  q, auto /*x*/) { return   uz.dz * q.val;         });
    auto t_after_matrix = std::chrono::steady_clock::now();

    auto t_before_boundary = std::chrono::steady_clock::now();
    assemble_facets(mesh.facets(), Wx, quad, G, [](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Wy, quad, G, [](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Wz, quad, G, [](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return 1/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.interior_facets(), Q, quad, G, [](auto p, auto q, auto /*x*/, const auto& face) {
        const auto  h = face.diameter;
        return h * jump(p).val * jump(q).val;
    });

    assemble_facets(mesh.facets(), Vx, Wx, quad, B, [eta](auto ux, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vx)), n) * jump(ux).val
               - dot(grad(avg(ux)), n) * jump(vx).val
               + eta/h * jump(ux).val * jump(vx).val;
    });
    assemble_facets(mesh.facets(), Vy, Wy, quad, B, [eta](auto uy, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vy)), n) * jump(uy).val
               - dot(grad(avg(uy)), n) * jump(vy).val
               + eta/h * jump(uy).val * jump(vy).val;
    });
    assemble_facets(mesh.facets(), Vz, Wz, quad, B, [eta](auto uz, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        return - dot(grad(avg(vz)), n) * jump(uz).val
               - dot(grad(avg(uz)), n) * jump(vz).val
               + eta/h * jump(uz).val * jump(vz).val;
    });
    assemble_facets(mesh.facets(), P, Wx, quad, B, [](auto p, auto vx, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{jump(vx).val, 0, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wy, quad, B, [](auto p, auto vy, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, jump(vy).val, 0};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), P, Wz, quad, B, [](auto p, auto vz, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  v = ads::point3_t{0, 0, jump(vz).val};
        return avg(p).val * dot(v, n);
    });
    assemble_facets(mesh.facets(), Vx, Q, quad, B, [](auto ux, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{jump(ux).val, 0, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, Q, quad, B, [](auto uy, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, jump(uy).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vz, Q, quad, B, [](auto uz, auto q, auto /*x*/, const auto& face) {
        const auto& n = face.normal;
        const auto  u = ads::point3_t{0, 0, jump(uz).val};
        return - dot(u, n) * avg(q).val;
    });
    auto t_after_boundary = std::chrono::steady_clock::now();

    fmt::print("Non-zeros: {}\n", problem.nonzero_entries());
    fmt::print("Computing RHS\n");

    auto t_before_rhs = std::chrono::steady_clock::now();
    assemble_rhs(Wx, quad, rhs, [&stokes](auto vx, auto x) { return vx.val * stokes.fx(x); });
    assemble_rhs(Wy, quad, rhs, [&stokes](auto vy, auto x) { return vy.val * stokes.fy(x); });
    assemble_rhs(Wz, quad, rhs, [&stokes](auto vz, auto x) { return vz.val * stokes.fz(x); });
    auto t_after_rhs = std::chrono::steady_clock::now();

    auto t_before_rhs_bnd = std::chrono::steady_clock::now();
    assemble_rhs(mesh.boundary_facets(), Wx, quad, rhs, [eta,&stokes](auto vx, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vx(x);
        return - dot(grad(vx), n) * g
               + eta/h * g * vx.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wy, quad, rhs, [eta,&stokes](auto vy, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vy(x);
        return - dot(grad(vy), n) * g
               + eta/h * g * vy.val;
    });
    assemble_rhs(mesh.boundary_facets(), Wz, quad, rhs, [eta,&stokes](auto vz, auto x, const auto& face) {
        const auto& n = face.normal;
        const auto  h = face.diameter;
        const auto  g = stokes.vz(x);
        return - dot(grad(vz), n) * g
               + eta/h * g * vz.val;
    });
    auto t_after_rhs_bnd = std::chrono::steady_clock::now();

    fmt::print("Solving\n");
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    auto sol_vx = ads::bspline_function3(&Vx, F.data() + N);
    auto sol_vy = ads::bspline_function3(&Vy, F.data() + N);
    auto sol_vz = ads::bspline_function3(&Vz, F.data() + N);
    auto sol_p  = ads::bspline_function3(&P,  F.data() + N);

    fmt::print("Computing error\n");

    auto t_before_err = std::chrono::steady_clock::now();
    auto mean   = integrate(mesh, quad, sol_p);
    auto err_vx = error(mesh, quad, L2{}, sol_vx, stokes.vx());
    auto err_vy = error(mesh, quad, L2{}, sol_vy, stokes.vy());
    auto err_vz = error(mesh, quad, L2{}, sol_vz, stokes.vz());
    auto err_p  = error(mesh, quad, L2{}, sol_p,  stokes.p(mean));
    auto err    = sum_norms(err_vx, err_vy, err_vz, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("vz error = {:.6}\n", err_vz);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

    auto r_vx = ads::bspline_function3(&Wx, F.data());
    auto r_vy = ads::bspline_function3(&Wy, F.data());
    auto r_vz = ads::bspline_function3(&Wz, F.data());
    auto r_p  = ads::bspline_function3(&Q,  F.data());

    auto norm_r_vx = norm(mesh, quad, L2{}, r_vx);
    auto J_r_vx = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vx(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_vy = norm(mesh, quad, L2{}, r_vy);
    auto J_r_vy = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vy(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_vz = norm(mesh, quad, L2{}, r_vz);
    auto J_r_vz = std::sqrt(integrate_facets(mesh.facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_vz(x, face));
        return 1/h * d * d;
    }));
    auto norm_r_p = norm(mesh, quad, L2{}, r_p);
    auto q_r_p = std::sqrt(integrate_facets(mesh.interior_facets(), mesh, quad, [&](auto x, const auto& face) {
        const auto h = face.diameter;
        auto d = jump(r_p(x, face));
        return h * d * d;
    }));
    auto res_norm = sum_norms(norm_r_vx, J_r_vx, norm_r_vy, J_r_vy, norm_r_vz, J_r_vz, norm_r_p, q_r_p);
    fmt::print("||r_vx|| = {:.6}\n", norm_r_vx);
    fmt::print(" |r_vx|  = {:.6}\n", J_r_vx);
    fmt::print("||r_vy|| = {:.6}\n", norm_r_vy);
    fmt::print(" |r_vy|  = {:.6}\n", J_r_vy);
    fmt::print("||r_vz|| = {:.6}\n", norm_r_vz);
    fmt::print(" |r_vz|  = {:.6}\n", J_r_vz);
    fmt::print("||r_p||  = {:.6}\n", norm_r_p);
    fmt::print(" |r_p|   = {:.6}\n", q_r_p);
    fmt::print("res norm = {:.6}\n", res_norm);

    // auto t_before_output = std::chrono::steady_clock::now();
    // save_to_file("result.data", sol_vx, sol_vy, sol_p);
    // auto t_after_output = std::chrono::steady_clock::now();

    auto as_ms = [](auto d) { return std::chrono::duration_cast<std::chrono::milliseconds>(d); };
    fmt::print("Matrix:  {:>8%Q %q}\n", as_ms(t_after_matrix   - t_before_matrix));
    fmt::print("Bndry:   {:>8%Q %q}\n", as_ms(t_after_boundary - t_before_boundary));
    fmt::print("RHS:     {:>8%Q %q}\n", as_ms(t_after_rhs      - t_before_rhs));
    fmt::print("RHS bd:  {:>8%Q %q}\n", as_ms(t_after_rhs_bnd  - t_before_rhs_bnd));
    fmt::print("Solver:  {:>8%Q %q}\n", as_ms(t_after_solver   - t_before_solver));
    fmt::print("Error:   {:>8%Q %q}\n", as_ms(t_after_err      - t_before_err));
    // fmt::print("Output:  {:>8%Q %q}\n", as_ms(t_after_output   - t_before_output));
}
