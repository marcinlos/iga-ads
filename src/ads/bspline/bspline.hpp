#ifndef ADS_SPLINE_HPP_
#define ADS_SPLINE_HPP_

#include <vector>
#include "ads/util/multi_array.hpp"

namespace ads {
namespace bspline {

using knot_vector = std::vector<double>;

struct basis {
    knot_vector knot;
    int degree;
    std::vector<double> points;

    basis(knot_vector kn, int degree): knot{ std::move(kn) }, degree{ degree } {
        points.push_back(knot[0]);
        for (auto i = 1u; i < knot.size(); ++ i) {
            if (knot[i] != knot[i - 1]) {
                points.push_back(knot[i]);
            }
        }
    }

    std::size_t knot_size() const {
        return knot.size();
    }

    int dofs() const {
        return knot_size() - degree - 1;
    }

    int dofs_per_element() const {
        return degree + 1;
    }

    int elements() const {
        return points.size() - 1;
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

struct basis_eval_array {
private:
    double* data_;
    int cols;

    int linearize_(int i, int j) const {
        return i * cols + j;
    }
public:
    basis_eval_array(double* data, int cols)
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

struct basis_eval_ctx {
private:
    std::vector<double> buffer;
    int right_offset;

public:
    basis_eval_array ndu;
    basis_eval_array a;

    basis_eval_ctx(int p)
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


struct eval_ctx : public basis_eval_ctx {
private:
    std::vector<double> buffer_;
public:

    eval_ctx(int p)
    : basis_eval_ctx(p)
    , buffer_(p + 1)
    { }

    double* basis_vals() {
        return buffer_.data();
    }
};

struct eval_ders_ctx : public basis_eval_ctx {
private:
    std::vector<double*> buffer_;
    int ders;
public:

    eval_ders_ctx(int p, int ders)
    : basis_eval_ctx(p + 1)
    , buffer_(p + 1)
    , ders(ders)
    {
        for (int i = 0; i <= ders; ++ i) {
            buffer_[i] = new double[p + 1];
        }
    }

    eval_ders_ctx(eval_ders_ctx&& other)
    : basis_eval_ctx(std::move(other))
    , buffer_(std::move(other.buffer_))
    , ders(other.ders)
    {
        other.buffer_.clear();
    }

    double** basis_vals() {
        return buffer_.data();
    }

    ~eval_ders_ctx() {
        for (auto b : buffer_) {
            delete[] b;
        }
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


basis create_basis(double a, double b, int degree, int elements, int repeated_nodes);

basis create_basis_C0(double a, double b, int degree, int elements);

int find_span(double x, const basis& b);

void eval_basis(int i, double x, const basis& b, double* out, basis_eval_ctx& ctx);

void eval_basis_with_derivatives(int i, double x, const basis& b, double** out, int d, basis_eval_ctx& ctx);

double eval(double x, const double* u, const basis& b, eval_ctx& ctx);

std::vector<int> first_nonzero_dofs(const basis& b);

std::vector<std::pair<int, int>> elements_supporting_dofs(const basis& b);


}
}

#endif /* ADS_SPLINE_HPP_ */
