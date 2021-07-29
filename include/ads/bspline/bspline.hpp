// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_BSPLINE_BSPLINE_HPP
#define ADS_BSPLINE_BSPLINE_HPP

#include <vector>

#include "ads/util.hpp"
#include "ads/util/multi_array.hpp"

namespace ads::bspline {

using knot_vector = std::vector<double>;

struct basis {
    knot_vector knot;
    int degree;
    std::vector<double> points;

    basis(knot_vector kn, int degree)
    : knot{std::move(kn)}
    , degree{degree} {
        points.push_back(knot[0]);
        for (auto i = 1; i < knot_size(); ++i) {
            if (knot[i] != knot[i - 1]) {
                points.push_back(knot[i]);
            }
        }
    }

    int knot_size() const { return narrow_cast<int>(knot.size()); }

    int dofs() const { return knot_size() - degree - 1; }

    int dofs_per_element() const { return degree + 1; }

    int elements() const { return static_cast<int>(points.size()) - 1; }

    int begin_idx() const { return degree; }

    int end_idx() const { return knot_size() - degree - 1; }

    double begin() const { return knot[begin_idx()]; }

    double end() const { return knot[end_idx()]; }

    double operator[](int i) const { return knot[i]; }
};

struct basis_eval_array {
private:
    double* data_;
    int cols;

    int linearize_(int i, int j) const { return i * cols + j; }

public:
    basis_eval_array(double* data, int cols)
    : data_(data)
    , cols(cols) { }

    double& operator()(int i, int j) {
        int idx = linearize_(i, j);
        return data_[idx];
    }
    double operator()(int i, int j) const {
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

    explicit basis_eval_ctx(int p)
    : buffer((p + 1) * (2 + (p + 1) + 2))
    , right_offset(p + 1)
    , ndu(buffer.data() + 2 * (p + 1), p + 1)
    , a(buffer.data() + (p + 1) * (2 + p + 1), p + 1) { }

    double* left() { return buffer.data(); }

    double* right() { return buffer.data() + right_offset; }
};

struct eval_ctx : public basis_eval_ctx {
private:
    std::vector<double> buffer_;

public:
    explicit eval_ctx(int p)
    : basis_eval_ctx(p)
    , buffer_(p + 1) { }

    double* basis_vals() { return buffer_.data(); }
};

struct eval_ders_ctx : public basis_eval_ctx {
private:
    std::vector<double*> buffer_;
    int ders;

public:
    eval_ders_ctx(int p, int ders)
    : basis_eval_ctx(p + 1)
    , buffer_(p + 1)
    , ders(ders) {
        for (int i = 0; i <= ders; ++i) {
            buffer_[i] = new double[p + 1];
        }
    }

    eval_ders_ctx(const eval_ders_ctx&) = delete;
    eval_ders_ctx& operator=(const eval_ders_ctx&) = delete;
    eval_ders_ctx(eval_ders_ctx&&) = delete;
    eval_ders_ctx& operator=(eval_ders_ctx&&) = delete;

    double** basis_vals() { return buffer_.data(); }

    ~eval_ders_ctx() {
        for (auto* b : buffer_) {
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
basis create_basis(double a, double b, int p, int elements);

basis create_basis(double a, double b, int p, int elements, int repeated_nodes);

basis create_basis_C0(double a, double b, int p, int elements);

int find_span(double x, const basis& b);

void eval_basis(int i, double x, const basis& b, double* out, basis_eval_ctx& ctx);

void eval_basis_with_derivatives(int i, double x, const basis& b, double** out, int d,
                                 basis_eval_ctx& ctx);

double eval(double x, const double* u, const basis& b, eval_ctx& ctx);

std::vector<int> first_nonzero_dofs(const basis& b);

std::vector<std::pair<int, int>> elements_supporting_dofs(const basis& b);

}  // namespace ads::bspline

#endif  // ADS_BSPLINE_BSPLINE_HPP
