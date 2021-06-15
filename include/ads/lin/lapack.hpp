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
    const double* x,
    in_int incx,
    in_double beta,
    double* y,
    in_int incy);

// SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
int dgemv_(
    const char* trans,
    in_int m,
    in_int n,
    in_double alpha,
    const double* a,
    in_int lda,
    const double* x,
    in_int incx,
    in_double beta,
    double* y,
    in_int incy);


int dgetrf_(
    in_int m,
    in_int n,
    double* a,
    in_int lda,
    out_int_array ipiv,
    out_int info);

int dgetrs_(
    const char* trans,
    in_int n,
    in_int nrhs,
    const double* a,
    in_int in_lda,
    in_int_array ipiv,
    double* b,
    in_int ldb,
    in_int info);

}

#endif /* ADS_LIN_LAPACK_HPP_ */
