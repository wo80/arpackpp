/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE blas1c.h.
   Interface to blas 1 and blas 2 FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arch.h"
#include "blas1f.h"

#ifndef BLAS1C_H
#define BLAS1C_H

// AXPY

inline void axpy(const ARint &n, const float &da, const float dx[],
                 const ARint &incx, float dy[], const ARint &incy) {
  F77NAME(saxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (float)

inline void axpy(const ARint &n, const double &da, const double dx[],
                 const ARint &incx, double dy[], const ARint &incy) {
  F77NAME(daxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (double)

#ifdef ARCOMP_H
inline void axpy(const ARint &n, const arcomplex<float> &da,
                 const arcomplex<float> dx[], const ARint &incx,
                 arcomplex<float> dy[], const ARint &incy) {
  F77NAME(caxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (arcomplex<float>)

inline void axpy(const ARint &n, const arcomplex<double> &da,
                 const arcomplex<double> dx[], const ARint &incx,
                 arcomplex<double> dy[], const ARint &incy) {
  F77NAME(zaxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (arcomplex<double>)
#endif

// COPY

inline void copy(const ARint &n, const float dx[], const ARint &incx,
                 float dy[], const ARint &incy) {
  if (dx != dy) F77NAME(scopy)(&n, dx, &incx, dy, &incy);
} // copy (float)

inline void copy(const ARint &n, const double dx[], const ARint &incx,
                 double dy[], const ARint &incy) {
  if (dx != dy) F77NAME(dcopy)(&n, dx, &incx, dy, &incy);
} // copy (double)

#ifdef ARCOMP_H
inline void copy(const ARint &n, const arcomplex<float> dx[], 
                 const ARint &incx, arcomplex<float> dy[], 
                 const ARint &incy) {
  if (dx != dy) F77NAME(ccopy)(&n, dx, &incx, dy, &incy);
} // copy (arcomplex<float>)

inline void copy(const ARint &n, const arcomplex<double> dx[], 
                 const ARint &incx, arcomplex<double> dy[], 
                 const ARint &incy) {
  if (dx != dy) F77NAME(zcopy)(&n, dx, &incx, dy, &incy);
} // copy (arcomplex<double>)
#endif

// DOT

inline float dot(const ARint &n, const float dx[], const ARint &incx,
                 const float dy[], const ARint &incy) {
  return F77NAME(sdot)(&n, dx, &incx, dy, &incy);
} // dot (float)

inline double dot(const ARint &n, const double dx[], const ARint &incx,
                  const double dy[], const ARint &incy) {
  return F77NAME(ddot)(&n, dx, &incx, dy, &incy);
} // dot (double)

// NRM2

inline float nrm2(const ARint &n, const float dx[], const ARint &incx) {
  return F77NAME(snrm2)(&n, dx, &incx);
} // nrm2 (float)

inline double nrm2(const ARint &n, const double dx[], const ARint &incx) {
  return F77NAME(dnrm2)(&n, dx, &incx);
} // nrm2 (double)

#ifdef ARCOMP_H
inline float nrm2(const ARint &n, const arcomplex<float> dx[], 
                  const ARint &incx) {
  return F77NAME(scnrm2)(&n, dx, &incx);
} // nrm2 (complex <float>)

inline double nrm2(const ARint &n, const arcomplex<double> dx[],
                   const ARint &incx) {
  return F77NAME(dznrm2)(&n, dx, &incx);
} // nrm2 (complex <double>)
#endif

// SCAL

inline void scal(const ARint &n, float &da, float dx[], const ARint &incx) {
  F77NAME(sscal)(&n, &da, dx, &incx);
} // scal (float)

inline void scal(const ARint &n, double &da, double dx[], const ARint &incx) {
  F77NAME(dscal)(&n, &da, dx, &incx);
} // scal (double)

#ifdef ARCOMP_H
inline void scal(const ARint &n, const arcomplex<float> &da,
                 arcomplex<float> dx[], const ARint &incx) {
  F77NAME(cscal)(&n, &da, dx, &incx);
} // scal (arcomplex<float>)

inline void scal(const ARint &n, const arcomplex<double> &da,
                 arcomplex<double> dx[], const ARint &incx) {
  F77NAME(zscal)(&n, &da, dx, &incx);
} // scal (arcomplex<double>)

inline void sscal(const ARint &n, const float &da, arcomplex<float> dx[],
                  const ARint &incx) {
  F77NAME(csscal)(&n, &da, dx, &incx);
} // sscal (arcomplex<float>)

inline void sscal(const ARint &n, const double &da, arcomplex<double> dx[],
                  const ARint &incx) {
  F77NAME(zdscal)(&n, &da, dx, &incx);
} // sscal (arcomplex<double>)
#endif

// GEMV

inline void gemv(const char* trans, const ARint &m, const ARint &n, 
                 const float &alpha, const float a[], const ARint &lda, 
                 const float x[], const ARint &incx, const float &beta, 
                 float y[], const ARint &incy) {
  F77NAME(sgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (float)

inline void gemv(const char* trans, const ARint &m, const ARint &n, 
                 const double &alpha, const double a[], const ARint &lda, 
                 const double x[], const ARint &incx, const double &beta, 
                 double y[], const ARint &incy) {
  F77NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (double)

#ifdef ARCOMP_H
inline void gemv(const char* trans, const ARint &m, 
                 const ARint &n, const arcomplex<float> &alpha,
                 const arcomplex<float> a[], const ARint &lda, 
                 const arcomplex<float> x[], const ARint &incx, 
                 const arcomplex<float> &beta, arcomplex<float> y[],
                 const ARint &incy) {
  F77NAME(cgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (arcomplex<float>)

inline void gemv(const char* trans, const ARint &m, 
                 const ARint &n, const arcomplex<double> &alpha,
                 const arcomplex<double> a[], const ARint &lda, 
                 const arcomplex<double> x[], const ARint &incx, 
                 const arcomplex<double> &beta, arcomplex<double> y[],
                 const ARint &incy) {
  F77NAME(zgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (arcomplex<double>)
#endif

// GBMV

inline void gbmv(const char* trans, const ARint &m, const ARint &n, 
                 const ARint &kl, const ARint &ku, const float &alpha,
                 const float a[], const ARint &lda, const float x[],
                 const ARint &incx, const float &beta, float y[],
                 const ARint &incy) {
  F77NAME(sgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (float)

inline void gbmv(const char* trans, const ARint &m, const ARint &n, 
                 const ARint &kl, const ARint &ku, const double &alpha,
                 const double a[], const ARint &lda, const double x[],
                 const ARint &incx, const double &beta, double y[],
                 const ARint &incy) {
  F77NAME(dgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (double)

#ifdef ARCOMP_H
inline void gbmv(const char* trans, const ARint &m, 
                 const ARint &n, const ARint &kl, 
                 const ARint &ku, const arcomplex<float> &alpha,
                 const arcomplex<float> a[], const ARint &lda, 
                 const arcomplex<float> x[], const ARint &incx, 
                 const arcomplex<float> &beta, arcomplex<float> y[],
                 const ARint &incy) {
  F77NAME(cgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (arcomplex<float>)

inline void gbmv(const char* trans, const ARint &m, 
                 const ARint &n, const ARint &kl, 
                 const ARint &ku, const arcomplex<double> &alpha,
                 const arcomplex<double> a[], const ARint &lda, 
                 const arcomplex<double> x[], const ARint &incx, 
                 const arcomplex<double> &beta, arcomplex<double> y[],
                 const ARint &incy) {
  F77NAME(zgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (arcomplex<double>)
#endif

// SBMV

inline void sbmv(const char* uplo, const ARint &n, const ARint &k, 
                 const float &alpha, const float a[], const ARint &lda, 
                 const float x[], const ARint &incx, const float &beta, 
                 float y[], const ARint &incy) {
  F77NAME(ssbmv)(uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
} // sbmv (float)

inline void sbmv(const char* uplo, const ARint &n, const ARint &k, 
                 const double &alpha, const double a[], const ARint &lda, 
                 const double x[], const ARint &incx, const double &beta, 
                 double y[], const ARint &incy) {
  F77NAME(dsbmv)(uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
} // sbmv (double)


#endif // BLAS1C_H
