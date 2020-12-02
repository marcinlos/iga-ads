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
            point  position;
            double normal;
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
            return {points_[i], normal};
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
            point       normal;
        };

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
                const auto [y, ny] = mesh_y_.facet(iy);
                const auto span    = mesh_x_.subinterval(ix);
                const auto normal  = point{0, ny};
                return {span, y, dir, normal};
            } else {
                assert(dir == orientation::vertical && "Invalid edge orientation");
                const auto [x, nx] = mesh_x_.facet(ix);
                const auto span    = mesh_y_.subinterval(iy);
                const auto normal  = point{nx, 0};
                return {span, x, dir, normal};
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

    struct facet_value {
        value_type avg, jump;
    };

    auto avg(const facet_value& val) noexcept -> value_type {
        return val.avg;
    }

    auto jump(const facet_value& val) noexcept -> value_type {
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

            auto operator ()(dof_index dof, point_index q, point normal) const noexcept -> facet_value {
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

        explicit bspline_function(const space* space, const double* coefficients)
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

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs(e)) {
                const auto u = eval(i, q);
                for (auto j : space.dofs(e)) {
                    const auto v = eval(j, q);
                    const auto iloc = space.local_index(i, e);
                    const auto jloc = space.local_index(j, e);
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

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : trial.dofs(e)) {
                const auto u = eval_trial(i, q);
                for (auto j : test.dofs(e)) {
                    const auto v = eval_test(j, q);
                    const auto iloc = trial.local_index(i, e);
                    const auto jloc = test.local_index(j, e);
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

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : space.dofs_on_facet(f)) {
                const auto u = eval(i, q, facet.normal);
                for (auto j : space.dofs_on_facet(f)) {
                    const auto v = eval(j, q, facet.normal);
                    const auto iloc = space.facet_local_index(i, f);
                    const auto jloc = space.facet_local_index(j, f);
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

        for (auto q : points.indices()) {
            const auto [x, w] = points.data(q);
            for (auto i : trial.dofs_on_facet(f)) {
                const auto u = eval_trial(i, q, facet.normal);
                for (auto j : test.dofs_on_facet(f)) {
                    const auto v = eval_test(j, q, facet.normal);
                    const auto iloc = trial.facet_local_index(i, f);
                    const auto jloc = test.facet_local_index(j, f);
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
    using value_of_f = decltype(f(typename Quad::point{}));
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

/////////////

void DG_laplace();
void DG_stokes();
void DGiGRM_stokes();


int main() {
    // DG_laplace();
    // DG_stokes();
    DGiGRM_stokes();
}

void DG_laplace() {
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
    fmt::print("DoFx: {}\n", n);

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

    auto p(double avg = 0) const noexcept {
        return [this,avg](auto x) { return avg + self()->p(x); };
    }

    auto fx() const noexcept {
        return [this](auto x) { return self()->fx(x); };
    }

    auto fy() const noexcept {
        return [this](auto x) { return self()->fy(x); };
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
    auto elems = 16;
    auto p     = 4;
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

    // auto M = [&problem](int off_row, int off_col) {
    //     return [&problem,off_row,off_col](int row, int col, double val) {
    //         if (val != 0) {
    //             problem.add(row + off_row + 1, col + off_col + 1, val);
    //         }
    //     };
    // };
    // auto rhs = [&F](int off) {
    //     return [&F,off](int J, double val) { F[off + J] += val; };
    // };
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
    auto err_p  = error(mesh, quad, L2{}, sol_p,  stokes.p(mean));
    auto err    = sum_norms(err_vx, err_vy, err_p);
    auto t_after_err = std::chrono::steady_clock::now();

    fmt::print("vx error = {:.6}\n", err_vx);
    fmt::print("vy error = {:.6}\n", err_vy);
    fmt::print("p  error = {:.6}\n", err_p);
    fmt::print("   error = {:.6}\n", err);

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
    fmt::print("Test  DoFs: {:7}\n", N);

    // Trial
    auto p_trial = 4;
    auto c_trial = 0;

    auto bx = ads::make_bspline_basis(xs, p_trial, c_trial);
    auto by = ads::make_bspline_basis(ys, p_trial, c_trial);

    auto trials = space_factory{};

    auto Vx = trials.next<ads::space>(&mesh, bx, by);
    auto Vy = trials.next<ads::space>(&mesh, bx, by);
    auto P  = trials.next<ads::space>(&mesh, bx, by);

    auto n = Vx.dof_count() + Vy.dof_count() + P.dof_count();
    fmt::print("Trial DoFs: {:7}\n", n);

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
    assemble_facets(mesh.facets(), Vx, Q, quad, B, [](auto ux, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{jump(ux).val, 0};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.facets(), Vy, Q, quad, B, [](auto uy, auto q, auto /*x*/, const auto& edge) {
        const auto& n = edge.normal;
        const auto  u = ads::point_t{0, jump(uy).val};
        return - dot(u, n) * avg(q).val;
    });
    assemble_facets(mesh.interior_facets(), P, Q, quad, B, [](auto p, auto q, auto /*x*/, const auto& edge) {
        const auto  h = length(edge.span);
        return h * jump(p).val * jump(q).val;
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

