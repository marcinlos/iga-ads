#include <vector>
#include <tuple>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <type_traits>
#include <chrono>
#include <mutex>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/counting_range.hpp>

#include <ads/util.hpp>
#include <ads/util/iter/product.hpp>
#include <ads/bspline/bspline.hpp>
#include <ads/bspline/eval.hpp>
#include <ads/quad/gauss.hpp>
#include <ads/lin/tensor.hpp>
#include <ads/util/function_value.hpp>
#include <mumps.hpp>


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


    class interval_mesh {
    private:
        partition points_;

    public:
        using point            = double;
        using element_index    = simple_index;
        using element_iterator = simple_index_iterator;
        using element_range    = simple_index_range;

        interval_mesh(partition points) noexcept
        : points_{std::move(points)}
        { }

        auto elements() const noexcept -> element_range {
            const auto elem_count = points_.size() - 1;
            return range(0, elem_count);
        }

        auto subinterval(element_index e) const noexcept -> interval {
            return ads::subinterval(points_, e);
        }
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
        using facet_index      = void;
        using facet_iterator   = void;
        using facet_range      = void;

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

        auto element(element_index e) const noexcept -> element_data {
            const auto [ix, iy] = e;
            const auto sx = mesh_x_.subinterval(ix);
            const auto sy = mesh_y_.subinterval(iy);
            return {sx, sy};
        }
    };


    using local_dof  = simple_index;
    using global_dof = simple_index;


    class bspline_space {
    private:
        bspline::basis   basis_;
        std::vector<int> first_dofs_;

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

        tensor_quadrature_points(interval_quadrature_points ptx, interval_quadrature_points pty)
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


    class quadrature {
    private:
        regular_mesh* mesh_;
        int point_count_;

    public:
        using point         = regular_mesh::point;
        using element_index = regular_mesh::element_index;
        using point_set     = tensor_quadrature_points;

        quadrature(regular_mesh* mesh, int point_count)
        : mesh_{mesh}
        , point_count_{point_count}
        { }

        auto coordinates(element_index e) const -> point_set {
            const auto element = mesh_->element(e);

            auto ptx = data_for_interval(element.span_x);
            auto pty = data_for_interval(element.span_y);

            return {std::move(ptx), std::move(pty)};
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


    class space {
    private:
        regular_mesh* mesh_;
        bspline_space space_x_;
        bspline_space space_y_;

        class evaluator;

    public:
        using point         = regular_mesh::point;
        using element_index = regular_mesh::element_index;
        using facet_index   = regular_mesh::facet_index;
        using dof_index     = index_types::index;
        using dof_iterator  = index_types::index_iterator;
        using dof_range     = index_types::index_range;

        space(regular_mesh* mesh, bspline::basis bx, bspline::basis by)
        : mesh_{mesh}
        , space_x_{std::move(bx)}
        , space_y_{std::move(by)}
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

        auto local_index(dof_index dof, element_index e) const noexcept -> local_dof {
            const auto idx = index_on_element(dof, e);

            const auto ndofs_x = space_x_.dofs_per_element();
            const auto ndofs_y = space_y_.dofs_per_element();

            return linearized(idx, {ndofs_x, ndofs_y});
        }

        auto global_index(dof_index dof) const noexcept -> global_dof {
            const auto ndofs_x = space_x_.dof_count();
            const auto ndofs_y = space_y_.dof_count();

            return linearized(dof, {ndofs_x, ndofs_y});
        }


        auto dof_evaluator(element_index e, const tensor_quadrature_points& points, int ders) const -> evaluator {
            auto data_x = evaluate_basis(points.xs(), space_x_, ders);
            auto data_y = evaluate_basis(points.ys(), space_y_, ders);

            return evaluator{this, e, ders, std::move(data_x), std::move(data_y)};
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

                const auto  Bx = vals_x_(qx, ix, 0);
                const auto dBx = vals_x_(qx, ix, 1);
                const auto  By = vals_y_(qy, iy, 0);
                const auto dBy = vals_y_(qy, iy, 1);

                const auto v  =  Bx *  By;
                const auto dx = dBx *  By;
                const auto dy =  Bx * dBy;

                return {v, dx, dy};
            }
        };

        auto index_on_element(dof_index dof, element_index e) const noexcept -> dof_index {
            const auto [ex, ey] = e;
            const auto [dx, dy] = dof;

            const auto ix = space_x_.local_index(dx, ex);
            const auto iy = space_y_.local_index(dy, ey);

            return {ix, iy};
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
        std::vector<double>       coefficients_;
        mutable bspline::eval_ctx ctx_x_;
        mutable bspline::eval_ctx ctx_y_;
        mutable std::mutex        ctx_lock_;

    public:
        using point  = space::point;

        explicit bspline_function(const space* space)
        : space_{space}
        , coefficients_(space->dof_count())
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

        auto data() noexcept -> double* {
            return coefficients_.data();
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

int main() {
    auto elems = 64;
    auto p     = 3;
    auto c     = 1;

    auto xs = ads::evenly_spaced(0.0, 1.0, elems);
    auto ys = ads::evenly_spaced(0.0, 1.0, elems);

    auto bx = ads::make_bspline_basis(xs, p, c);
    auto by = ads::make_bspline_basis(ys, p, c);

    auto mesh = ads::regular_mesh{xs, ys};
    auto space = ads::space{&mesh, bx, by};
    auto quad = ads::quadrature{&mesh, p + 1};

    auto n = space.dof_count();
    std::cout << "DoFs: " << n << std::endl;

    auto F = ads::bspline_function(&space);
    auto problem = mumps::problem{F.data(), n};
    auto solver  = mumps::solver{};

    auto diff = [](auto before, auto after) {
        return std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
    };

    auto t_before_matrix = std::chrono::steady_clock::now();

    for (auto e : mesh.elements()) {
        auto points = quad.coordinates(e);
        auto eval   = space.dof_evaluator(e, points, 2);

        auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 2>{{n, n}};

        for (auto q : points.indices()) {
            auto [x, w] = points.data(q);

            for (auto i : space.dofs(e)) {
                auto u = eval(i, q);
                for (auto j : space.dofs(e)) {
                    auto v = eval(j, q);

                    auto iloc = space.local_index(i, e);
                    auto jloc = space.local_index(j, e);
                    // M(jloc, iloc) += (u.dx * v.dx + u.dy * v.dy) * w;
                    M(jloc, iloc) += u.val * v.val * w;
                }
            }
        }

        for (auto i : space.dofs(e)) {
            auto iloc = space.local_index(i, e);
            auto I    = space.global_index(i);
            for (auto j : space.dofs(e)) {
                auto jloc = space.local_index(j, e);
                auto J    = space.global_index(j);

                problem.add(J + 1, I + 1, M(jloc, iloc));
            }
        }
    }
    auto t_after_matrix = std::chrono::steady_clock::now();

    std::cout << "Non-zeros: " << problem.nonzero_entries() << std::endl;
    std::cout << "Computing RHS" << std::endl;

    auto func = [](double x, double y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    };

    auto t_before_rhs = std::chrono::steady_clock::now();
    for (auto e : mesh.elements()) {
        auto points = quad.coordinates(e);
        auto eval   = space.dof_evaluator(e, points, 2);

        auto n = space.dof_count(e);
        auto M = ads::lin::tensor<double, 1>{{n}};

        for (auto q : points.indices()) {
            auto [X, w] = points.data(q);
            auto [x, y] = X;

            for (auto j : space.dofs(e)) {
                auto v = eval(j, q);

                auto jloc = space.local_index(j, e);
                M(jloc) += v.val * func(x, y) * w;
            }
        }
        for (auto j : space.dofs(e)) {
            auto jloc = space.local_index(j, e);
            auto J    = space.global_index(j);
            F.data()[J] += M(jloc);
        }
    }
    auto t_after_rhs = std::chrono::steady_clock::now();

    std::cout << "Solving" << std::endl;
    auto t_before_solver = std::chrono::steady_clock::now();
    solver.solve(problem);
    auto t_after_solver = std::chrono::steady_clock::now();

    std::cout << "Evaluating" << std::endl;

    auto t_before_eval = std::chrono::steady_clock::now();

    auto err = 0.0;
    // for (auto e : mesh.elements()) {
    //     auto [sx, sy] = mesh.element(e);
    //     auto x = (sx.left + sx.right) / 2;
    //     auto y = (sy.left + sy.right) / 2;

    //     auto val = F({x, y});
    //     auto exact = func(x, y);
    //     auto d = val - exact;
    //     err = std::max(err, std::abs(d));
    // }
    for (auto e : mesh.elements()) {
        auto points = quad.coordinates(e);
        for (auto q : points.indices()) {
            auto [X, w] = points.data(q);
            auto [x, y] = X;

            auto val   = F(X);
            auto exact = func(x, y);
            auto d = val - exact;
            err += d * d * w;
        }
    }
    err = std::sqrt(err);

    auto t_after_eval = std::chrono::steady_clock::now();

    std::cout << "error = " << err << std::endl;
    std::cout << "Matrix: " << std::setw(6) << diff(t_before_matrix, t_after_matrix) << " ms" << std::endl;
    std::cout << "RHS:    " << std::setw(6) << diff(t_before_rhs, t_after_rhs)       << " ms" << std::endl;
    std::cout << "Solver: " << std::setw(6) << diff(t_before_solver, t_after_solver) << " ms" << std::endl;
    std::cout << "Eval:   " << std::setw(6) << diff(t_before_eval, t_after_eval)     << " ms" << std::endl;


    // auto quad = gauss_quad(space);

    // for (auto e : mesh.elements()) {
    //     auto eval = space.evaluator(e, quad.points(e), ders=2);
    //     auto J = mesh.jacobian(e);
    //     for (auto q : quad.points(e)) {
    //         auto w = quad.weight(e, q);
    //         for (auto i : space.dofs(e)) {
    //             for (auto j : space.dofs(e)) {
    //                 auto u = eval(j, q);
    //                 auto v = eval(i, q);
    //                 auto val = u.dx * v.dx + u.dy * v.dy;

    //                 // possibly...
    //                 auto ii = space.index(i);
    //                 auto jj = space.index(j);
    //                 M(ii, jj) += val * w * J;
    //             }
    //         }
    //     }
    // }

    // for (auto f : mesh.facets()) {
    //     auto eval = space.evaluator(f, quad.points(f), ders=2);
    //     auto J = mesh.jacobian(f);
    //     for (auto q : quad.points(f)) {
    //         auto w = quad.weight(f, q);
    //         for (auto i : space.dofs(f)) {
    //             for (auto j : space.dofs(f)) {
    //                 auto u = eval(j, q);
    //                 auto v = eval(i, q);

    //                 // possibly...
    //                 auto ii = space.index(i);
    //                 auto jj = space.index(j);
    //                 M(ii, jj) = h * jump(u).val * jump(v).val;
    //             }
    //         }
    //     }
    // }
}

