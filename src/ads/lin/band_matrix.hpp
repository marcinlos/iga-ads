#ifndef ADS_LIN_BANDED_MATRIX_HPP_
#define ADS_LIN_BANDED_MATRIX_HPP_

#include <cstddef>
#include <vector>
#include <iostream>
#include <iomanip>

namespace ads {
namespace lin {

class band_matrix {
private:
    std::vector<double> data_;

public:
    int kl;
    int ku;
    int n;

    band_matrix(int kl, int ku, int n)
    : data_(array_size_(kl, ku, n))
    , kl(kl)
    , ku(ku)
    , n(n)
    { }

    double& operator ()(int i, int j) {
        int row = kl + ku + i - j;
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

    int cols() const {
        return n;
    }

    int rows() const {
        return n;
    }

    double* data() {
        return data_.data();
    }

    const double* data() const {
        return data_.data();
    }

private:

    int column_size_() const {
        return kl + (ku + 1 + kl);
    }

    static int array_size_(int kl, int ku, int n) {
        int rows = kl + (ku + 1 + kl);
        int cols = n;
        return rows * cols;
    }
};


inline std::ostream& operator <<(std::ostream& os, const band_matrix& M) {
    for (int i = 0; i < M.rows(); ++ i) {
        for (int j = 0; j < M.cols(); ++ j) {
            os << std::setw(12) << M(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}


}
}

#endif /* ADS_LIN_BANDED_MATRIX_HPP_ */
