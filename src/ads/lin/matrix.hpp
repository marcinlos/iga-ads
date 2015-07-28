#ifndef ADS_LIN_MATRIX_HPP_
#define ADS_LIN_MATRIX_HPP_

#include <vector>
#include <cmath>
#include "ads/lin/print_matrix.hpp"

namespace ads {
namespace lin {

class matrix {
private:
    const int rows_;
    const int cols_;

    std::vector<double> data_;

public:
    matrix(int rows, int cols)
    : rows_(rows)
    , cols_(cols)
    , data_(rows * cols)
    { }

    double& operator ()(int i, int j) {
        return data_[j * rows_ + i];
    }

    double operator ()(int i, int j) const {
        auto& self = const_cast<matrix&>(*this);
        return self(i, j);
    }

    int cols() const {
        return cols_;
    }

    int rows() const {
        return rows_;
    }

    double* data() {
        return data_.data();
    }

    const double* data() const {
        return data_.data();
    }
};

inline std::ostream& operator <<(std::ostream& os, const matrix& m) {
    return print_matrix(os, m);
}

inline bool equal(const matrix& a, const matrix& b, double eps) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        return false;
    }
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = 0; j < a.cols(); ++j) {
            if (std::fabs(a(i, j) - b(i, j)) >= eps) {
                return false;
            }
        }
    }
    return true;
}

}
}

#endif /* ADS_LIN_MATRIX_HPP_ */
