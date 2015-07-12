#ifndef _CUBLASCPP_H_
#define _CUBLASCPP_H_

#include <cublas.h>

/**
 *\file CuBlasCpp.h
 *\note all the datas (*A,*B,*C) should pointed to the device memory,
 *and allocated before calling the functions.
*/


namespace LSW_CUDA_MATH{

  bool gemm(char transa, char transb, int m, int 
  			   n, int k, double alpha, double *a, int lda, 
  			   double *b, int ldb, double beta, double *c__, 
  			   int ldc);

  bool gemm(char transa, char transb, int m, int 
  			   n, int k, float alpha, float *a, int lda, 
  			   float *b, int ldb, float beta, float *c__, 
  			   int ldc);

  bool axpy(int n, double da, double *dx, 
			   int incx, double *dy, int incy);

  bool axpy(int n,float da,float *dx, 
			   int incx, float *dy, int incy);

  bool gemv(char trans, int m, int n, double 
			   alpha, double *a, int lda, double *x, int incx, 
			   double beta, double *y, int incy);

  bool gemv(char trans, int m, int n, float 
			   alpha, float *a, int lda, float *x, int incx, 
			   float beta, float *y, int incy);

  bool scal(int n, double sa, double *sx, int incx);

  bool scal(int n, float sa, float *sx, int incx);

}


#endif /* _CUBLASCPP_H_ */
