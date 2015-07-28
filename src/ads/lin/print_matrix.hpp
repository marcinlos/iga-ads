#ifndef ADS_LIN_PRINT_MATRIX_HPP_
#define ADS_LIN_PRINT_MATRIX_HPP_

#include <iostream>
#include <iomanip>

namespace ads {
namespace lin {

template <typename Mat>
inline std::ostream& print_matrix(std::ostream& os, const Mat& m) {

    int precision = os.precision();
    os << std::setprecision(6);

    std::ios::fmtflags original = os.flags();
    os.setf(std::ios::fixed, std:: ios::floatfield);

    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            os << std::setw(11) << m(i, j);
        }
        os << std::endl;
    }

    os << std::setprecision(precision);
    os.flags(original);
    return os;
}

}
}

#endif /* ADS_LIN_PRINT_MATRIX_HPP_ */
