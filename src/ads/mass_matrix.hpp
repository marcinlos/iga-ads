#ifndef ADS_MASS_MATRIX_HPP_
#define ADS_MASS_MATRIX_HPP_

#include "ads/lin/band_matrix.hpp"
#include "ads/basis_data.hpp"


namespace ads {

    void gram_matrix_1d(lin::band_matrix& M, const basis_data& d);

}


#endif /* ADS_MASS_MATRIX_HPP_ */
