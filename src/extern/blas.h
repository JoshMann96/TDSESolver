#pragma once

#define dcomplex double __complex__
#define lapack_int int
#define lapack_complex_double dcomplex

#ifdef __cplusplus
extern "C" {
#endif

//lapacke_mangling.h
#ifndef LAPACK_HEADER_INCLUDED
#define LAPACK_HEADER_INCLUDED

#ifndef LAPACK_GLOBAL
#if defined(LAPACK_GLOBAL_PATTERN_LC) || defined(ADD_)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#elif defined(LAPACK_GLOBAL_PATTERN_UC) || defined(UPPER)
#define LAPACK_GLOBAL(lcname,UCNAME)  UCNAME
#elif defined(LAPACK_GLOBAL_PATTERN_MC) || defined(NOCHANGE)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname
#else
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#endif
#endif

#endif

//selected from lapack.h
#define LAPACK_zhpevx_base LAPACK_GLOBAL(zhpevx,ZHPEVX)
void LAPACK_zhpevx_base(
    char const* jobz, char const* range, char const* uplo,
    lapack_int const* n,
    lapack_complex_double* AP,
    double const* vl,
    double const* vu, lapack_int const* il, lapack_int const* iu,
    double const* abstol, lapack_int* m,
    double* W,
    lapack_complex_double* Z, lapack_int const* ldz,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* iwork, lapack_int* IFAIL,
    lapack_int* info
#ifdef LAPACK_FORTRAN_STRLEN_END
    , size_t, size_t, size_t
#endif
);
#ifdef LAPACK_FORTRAN_STRLEN_END
    #define LAPACK_zhpevx(...) LAPACK_zhpevx_base(__VA_ARGS__, 1, 1, 1)
#else
    #define LAPACK_zhpevx(...) LAPACK_zhpevx_base(__VA_ARGS__)
#endif

#define LAPACK_dlamch_base LAPACK_GLOBAL(dlamch,DLAMCH)
double LAPACK_dlamch_base(
    char const* cmach
#ifdef LAPACK_FORTRAN_STRLEN_END
    , size_t
#endif
);
#ifdef LAPACK_FORTRAN_STRLEN_END
    #define LAPACK_dlamch(...) LAPACK_dlamch_base(__VA_ARGS__, 1)
#else
    #define LAPACK_dlamch(...) LAPACK_dlamch_base(__VA_ARGS__)
#endif

#define LAPACK_zgtsv LAPACK_GLOBAL(zgtsv,ZGTSV)
void LAPACK_zgtsv(
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double* DL,
    lapack_complex_double* D,
    lapack_complex_double* DU,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zgtsvx LAPACK_GLOBAL(zgtsvx,ZGTSVX)
void LAPACK_zgtsvx(
    char const* fact, char const* trans,
    lapack_int const* n, lapack_int const* nrhs,
    lapack_complex_double const* DL,
    lapack_complex_double const* D,
    lapack_complex_double const* DU,
    lapack_complex_double* DLF,
    lapack_complex_double* DF,
    lapack_complex_double* DUF,
    lapack_complex_double* DU2, lapack_int* ipiv,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

#define LAPACK_zptsv LAPACK_GLOBAL(zptsv,ZPTSV)
void LAPACK_zptsv(
    lapack_int const* n, lapack_int const* nrhs,
    double* D,
    lapack_complex_double* E,
    lapack_complex_double* B, lapack_int const* ldb,
    lapack_int* info );

#define LAPACK_zptsvx LAPACK_GLOBAL(zptsvx,ZPTSVX)
void LAPACK_zptsvx(
    char const* fact,
    lapack_int const* n, lapack_int const* nrhs,
    double const* D,
    lapack_complex_double const* E,
    double* DF,
    lapack_complex_double* EF,
    lapack_complex_double const* B, lapack_int const* ldb,
    lapack_complex_double* X, lapack_int const* ldx,
    double* rcond,
    double* ferr,
    double* berr,
    lapack_complex_double* work,
    double* rwork,
    lapack_int* info );

//selected from cblas.h
extern void cblas_zscal(const int N, const void *alpha, void *X, const int incX);
extern void cblas_dscal(const int N, const double alpha, double *X, const int incX);
extern void cblas_zdscal(const int N, const double alpha, void *X, const int incX);
extern void cblas_dcopy(const int N, const double *X, const int incX, 
                 double *Y, const int incY);
extern void cblas_zcopy(const int N, const void *X, const int incX, 
                 void *Y, const int incY);
extern void cblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
extern void cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);


enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
extern void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
				 double *Y, const int incY);

#ifdef __cplusplus
}
#endif