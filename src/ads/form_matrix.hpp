#ifndef ADS_FORM_MATRIX_HPP_
#define ADS_FORM_MATRIX_HPP_

#include "ads/lin/band_matrix.hpp"
#include "ads/lin/dense_matrix.hpp"
#include "ads/basis_data.hpp"


namespace ads {

    void gram_matrix_1d(lin::band_matrix& M, const basis_data& d);

    void stiffness_matrix_1d(lin::band_matrix& M, const basis_data& d);

    void advection_matrix_1d(lin::band_matrix& M, const basis_data& d);


    void gram_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

    void stiffness_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

    void advection_matrix_1d(lin::dense_matrix& M, const basis_data& U, const basis_data& V);

}


#endif /* ADS_FORM_MATRIX_HPP_ */
