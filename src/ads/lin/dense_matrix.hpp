#ifndef ADS_LIN_DENSE_MATRIX_HPP_
#define ADS_LIN_DENSE_MATRIX_HPP_

#include <vector>
#include <iostream>
#include <iomanip>

namespace ads {
namespace lin {

class dense_matrix {
private:
    std::vector<double> data_;
    int rows_;
    int cols_;

public:
    dense_matrix(int rows, int cols): data_(rows * cols), rows_(rows), cols_(cols) { }

    double& operator ()(int i, int j) {
        return data_[j * rows_ + i];
    }

    double operator ()(int i, int j) const {
        return data_[j * rows_ + i];
    }

    double* data() { return data_.data(); }

    const double* data() const { return data_.data(); }

    int rows() const { return rows_; }

    int cols() const { return cols_; }

    void zero() {
        std::fill(begin(data_), end(data_), 0);
    }

    int size() const { return cols_ * rows_; }

    int size(int dim) const {
        return dim == 0 ? rows_ : cols_;
    }
};

inline std::ostream& operator <<(std::ostream& os, const dense_matrix& M) {
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

#endif /* ADS_LIN_DENSE_MATRIX_HPP_ */
