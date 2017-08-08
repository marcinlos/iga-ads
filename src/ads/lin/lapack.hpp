#ifndef ADS_LIN_LAPACK_HPP_
#define ADS_LIN_LAPACK_HPP_

// LAPACK routines
extern "C" {

using in_int = const int*;
using in_double = const double*;
using in_int_array = const int*;
using out_int = int*;
using out_int_array = int*;

int dgbtrf_(
    in_int m,
    in_int n,
    in_int kl,
    in_int ku,
    double* ab,
    in_int ldab,
    out_int_array ipiv,
    out_int info);

int dgbtrs_(
    const char* trans,
    in_int n,
    in_int kl,
    in_int ku,
    in_int nrhs,
    const double* ab,
    in_int ldab,
    in_int_array ipiv,
    double* b,
    in_int ldb,
    out_int info);

int dgbmv_(
    const char* trans,
    in_int m,
    in_int n,
    in_int kl,
    in_int ku,
    in_double alpha,
    const double* a,
    in_int lda,
    double* x,
    in_int incx,
    in_double beta,
    const double* y,
    in_int incy);

}


#endif /* ADS_LIN_LAPACK_HPP_ */
