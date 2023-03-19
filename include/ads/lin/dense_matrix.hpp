// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_DENSE_MATRIX_HPP
#define ADS_LIN_DENSE_MATRIX_HPP

#include <iomanip>
#include <iostream>
#include <vector>

#include "ads/lin/lapack.hpp"

namespace ads::lin {

class dense_matrix {
private:
    std::vector<double> data_;
    int rows_;
    int cols_;

public:
    dense_matrix(int rows, int cols)
    : data_(rows * cols)
    , rows_(rows)
    , cols_(cols) { }

    double& operator()(int i, int j) { return data_[j * rows_ + i]; }

    double operator()(int i, int j) const { return data_[j * rows_ + i]; }

    double* data() { return data_.data(); }

    const double* data() const { return data_.data(); }

    int rows() const { return rows_; }

    int cols() const { return cols_; }

    void zero() { std::fill(begin(data_), end(data_), 0); }

    int size() const { return cols_ * rows_; }

    int size(int dim) const { return dim == 0 ? rows_ : cols_; }
};

inline std::ostream& operator<<(std::ostream& os, const dense_matrix& M) {
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            os << std::setw(12) << M(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

template <typename Vec1, typename Vec2>
inline void multiply(const dense_matrix& M, const Vec1& x, Vec2& y, int count = 1,
                     const char* transpose = "N") {
    double alpha = 1;
    double beta = 0;
    int incx = 1;
    int incy = 1;
    int lda = M.rows();
    int rows = M.rows();
    int cols = M.cols();

    for (int i = 0; i < count; ++i) {
        auto in = x.data() + M.rows() * i;
        auto out = y.data() + M.cols() * i;
        dgemv_(transpose, &rows, &cols, &alpha, M.data(), &lda, in, &incx, &beta, out, &incy);
    }
}

inline void multiply(const dense_matrix& A, const dense_matrix& B, dense_matrix& out,
                     const char* transpose = "N") {
    multiply(A, B, out, B.cols(), transpose);
}

}  // namespace ads::lin

#endif  // ADS_LIN_DENSE_MATRIX_HPP
