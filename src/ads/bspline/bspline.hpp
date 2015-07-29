#ifndef ADS_SPLINE_HPP_
#define ADS_SPLINE_HPP_

#include <vector>

namespace ads {
namespace bspline {

using knot_vector = std::vector<double>;

struct basis {
    knot_vector knot;
    int degree;

    int knot_size() const {
        return knot.size();
    }

    int dofs() const {
        return knot_size() - degree - 1;
    }

    int elements() const {
        return dofs() - degree;
    }

    int begin_idx() const {
        return degree;
    }

    int end_idx() const {
        return knot_size() - degree - 1;
    }

    double begin() const {
        return knot[begin_idx()];
    }

    double end() const {
        return knot[end_idx()];
    }

    double operator [](int i) const {
        return knot[i];
    }
};

struct eval_array {
private:
    double* data_;
    int cols;

    int linearize_(int i, int j) const {
        return i * cols + j;
    }
public:
    eval_array(double* data, int cols)
    : data_(data)
    , cols(cols)
    { }

    double& operator ()(int i, int j) {
        int idx = linearize_(i, j);
        return data_[idx];
    }
    double operator ()(int i, int j) const {
        int idx = linearize_(i, j);
        return data_[idx];
    }
};

struct eval_ctx {
private:
    std::vector<double> buffer;
    int right_offset;

public:
    eval_array ndu;
    eval_array a;

    eval_ctx(int p)
    : buffer((p + 1) * (2 + (p + 1) + 2))
    , right_offset(p + 1)
    , ndu(buffer.data() + 2 * (p + 1), p + 1)
    , a(buffer.data() + (p + 1) * (2 + p + 1), p + 1)
    { }

    double* left() {
        return buffer.data();
    }

    double* right() {
        return buffer.data() + right_offset;
    }
};

/**
 * Create B-spline basis of specified order with given number of elements on given interval
 *
 * @param a        - left interval endpoint
 * @param b        - right interval endpoint
 * @param p        - order of B-spline basis
 * @param elements - number of elements (spans)
 */
basis create_basis(double a, double b, int degree, int elements);

int find_span(double x, const basis& b);

void eval_basis(int i, double x, const basis& b, double* out, eval_ctx& ctx);

void eval_basis_with_derivatives(int i, double x, const basis& b, double** out, int d, eval_ctx& ctx);

}
}

#endif /* ADS_SPLINE_HPP_ */
