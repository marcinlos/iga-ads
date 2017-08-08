#ifndef ADS_LIN_BANDED_MATRIX_HPP_
#define ADS_LIN_BANDED_MATRIX_HPP_

#include <cstddef>
#include <vector>
#include <iomanip>
#include "ads/lin/lapack.hpp"

namespace ads {
namespace lin {

class band_matrix {
private:
    std::vector<double> data_;

public:
    int kl;
    int ku;
    int rows;
    int cols;
    int row_offset;

    band_matrix(int kl, int ku, int n): band_matrix(kl, ku, n, n, kl) { }

    band_matrix(int kl, int ku, int rows, int cols, int row_offset = 0)
    : data_(array_size_(kl, ku, cols, row_offset))
    , kl(kl)
    , ku(ku)
    , rows(rows)
    , cols(cols)
    , row_offset(row_offset)
    { }

    double& operator ()(int i, int j) {
        int row = row_offset + ku + i - j;
        int col = j;
        return data_[col * column_size_() + row];
    }

    double operator ()(int i, int j) const {
        if (inside_band(i, j)) {
            auto& self = const_cast<band_matrix&>(*this);
            return self(i, j);
        }
        return 0;
    }

    bool inside_band(int i, int j) const {
        return i - j <= kl && j - i <= ku;
    }

    double* full_buffer() {
        return data_.data();
    }

    const double* full_buffer() const {
        return data_.data();
    }

    double* data() {
        return full_buffer() + row_offset * cols;
    }

    const double* data() const {
        return full_buffer() + row_offset * cols;
    }

private:

    int column_size_() const {
        return row_offset + (ku + 1 + kl);
    }

    static int array_size_(int kl, int ku, int cols, int offset) {
        int rows = offset + (ku + 1 + kl);
        return rows * cols;
    }
};


inline std::ostream& operator <<(std::ostream& os, const band_matrix& M) {
    for (int i = 0; i < M.rows; ++ i) {
        for (int j = 0; j < M.cols; ++ j) {
            os << std::setw(12) << M(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}


}
}

#endif /* ADS_LIN_BANDED_MATRIX_HPP_ */
