// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_LIN_BAND_MATRIX_HPP
#define ADS_LIN_BAND_MATRIX_HPP

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "ads/lin/dense_matrix.hpp"
#include "ads/lin/lapack.hpp"

namespace ads::lin {

class band_matrix {
private:
    std::vector<double> data_;

public:
    int kl = 0;
    int ku = 0;
    int rows = 0;
    int cols = 0;
    int row_offset = 0;

    band_matrix() = default;

    band_matrix(int kl, int ku, int n)
    : band_matrix(kl, ku, n, n, kl) { }

    band_matrix(int kl, int ku, int rows, int cols, int row_offset = 0)
    : data_(array_size_(kl, ku, cols, row_offset))
    , kl(kl)
    , ku(ku)
    , rows(rows)
    , cols(cols)
    , row_offset(row_offset) { }

    double& operator()(int i, int j) {
        assert(inside_band(i, j));
        return data_[index_(i, j)];
    }

    double operator()(int i, int j) const {
        if (inside_band(i, j)) {
            return data_[index_(i, j)];
        }
        return 0;
    }

    bool inside_band(int i, int j) const { return i - j <= kl && j - i <= ku; }

    double* full_buffer() { return data_.data(); }

    const double* full_buffer() const { return data_.data(); }

    double* data() { return full_buffer() + row_offset * cols; }

    const double* data() const { return full_buffer() + row_offset * cols; }

    void zero() { std::fill(begin(data_), end(data_), 0); }

    int column_size() const { return row_offset + (ku + 1 + kl); }

private:
    int index_(int i, int j) const {
        int row = row_offset + ku + i - j;
        int col = j;
        return col * column_size() + row;
    }

    static int array_size_(int kl, int ku, int cols, int offset) {
        int rows = offset + (ku + 1 + kl);
        return rows * cols;
    }
};

inline std::ostream& operator<<(std::ostream& os, const band_matrix& M) {
    for (int i = 0; i < M.rows; ++i) {
        for (int j = 0; j < M.cols; ++j) {
            os << std::setw(12) << M(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

template <typename Vec1, typename Vec2>
inline void multiply(const band_matrix& M, const Vec1& x, Vec2& y, int count = 1,
                     const char* transpose = "N") {
    double alpha = 1;
    double beta = 0;
    int incx = 1;
    int incy = 1;
    int lda = M.column_size();

    for (int i = 0; i < count; ++i) {
        auto in = x.data() + M.cols * i;
        auto out = y.data() + M.rows * i;
        dgbmv_(transpose, &M.rows, &M.cols, &M.kl, &M.ku, &alpha, M.data(), &lda, in, &incx, &beta,
               out, &incy);
    }
}

inline void multiply(const band_matrix& A, const dense_matrix& B, dense_matrix& out,
                     const char* transpose = "N") {
    multiply(A, B, out, B.cols(), transpose);
}

inline void to_dense(const band_matrix& M, dense_matrix& out) {
    for (int i = 0; i < M.rows; ++i) {
        for (int j = 0; j < M.cols; ++j) {
            out(i, j) = M(i, j);
        }
    }
}

}  // namespace ads::lin

#endif  // ADS_LIN_BAND_MATRIX_HPP
