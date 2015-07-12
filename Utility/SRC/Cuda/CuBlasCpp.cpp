#include "CuBlasCpp.h"

namespace LSW_CUDA_MATH{

  bool gemm(char transa, char transb, int m, int 
  			   n, int k, double alpha, double *a, int lda, 
  			   double *b, int ldb, double beta, double *c__, 
  			   int ldc){
  	 cublasDgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c__,ldc);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool gemm(char transa, char transb, int m, int 
  			   n, int k, float alpha, float *a, int lda, 
  			   float *b, int ldb, float beta, float *c__, 
  			   int ldc){
  	 cublasSgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c__,ldc);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool axpy(int n, double da, double *dx, 
			   int incx, double *dy, int incy){
	 cublasDaxpy(n,da,dx,incx,dy,incy);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool axpy(int n,float da,float *dx, 
			   int incx, float *dy, int incy){
	 cublasSaxpy(n,da,dx,incx,dy,incy);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool gemv(char trans, int m, int n, double 
			   alpha, double *a, int lda, double *x, int incx, 
			   double beta, double *y, int incy){
	 cublasDgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool gemv(char trans, int m, int n, float 
			   alpha, float *a, int lda, float *x, int incx, 
			   float beta, float *y, int incy){
	 cublasSgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool scal(int n, double sa, double *sx, int incx){
	 cublasDscal(n,sa,sx,incx);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

  bool scal(int n, float sa, float *sx, int incx){
	 cublasSscal(n,sa,sx,incx);
	 return (cublasGetError () == CUBLAS_STATUS_SUCCESS);
  }

}
