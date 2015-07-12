#ifndef _WARPEDCUBLAS_H_
#define _WARPEDCUBLAS_H_

#include "CuBlasCpp.h"

/**
 *\file WarpedCuBlas.h
 *\note all the datas (*A,*B,*C) should pointed to the device memory,
 *and allocated before calling the functions.
*/

namespace LSW_CUDA_MATH{
  
  ///C = alpha*T(a)*T(b)+beta*C
  ///\note parameters allocated outside and pointed to the device memory
  template<typename real>
  bool gemm(const real*A,int r_a,int c_a,const real* B,int r_b,int c_b,real *C, char transa='n',char transb='n',real alpha=1.0f,real beta=0.0f){

	bool check_valid = !(r_a<=0||c_a<=0||r_b<=0||c_b<=0||A==NULL||B==NULL||C==NULL);
	if(!check_valid)
	  return false;

	int m = (transa == 'n' ? r_a : c_a);//rows of op(A)
	int n = (transb == 'n' ? c_b :r_b);//columns of op(B)
	int k = (transa == 'n' ? c_a:r_a);//columns of op(A)
	int lda = r_a;//rows of A
	int ldb = r_b;//rows of B
	int ldc = m;//m

	real *a = const_cast<real *>(A);
	real *b = const_cast<real *>(B);

	int op_b_r = (transb == 'n' ? r_b:c_b);//rows of op(B):just for test
	if(op_b_r != k)
	  return false;
	//calling the blas function
	return gemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,C,ldc);  
  }

  ///Y = alpha*T(A)*X + beta*Y
  ///\note parameters allocated outside and pointed to the device memory
  template<typename real>
  bool gemv(const real *A,int r_a,int c_a, const real*X,real *Y, char transa='n',real alpha=1.0f,real beta=0.0f){

	if(A == NULL || X == NULL || Y == NULL || Y == X)
	  return false;

	int m = r_a;
	int n = c_a;
	int lda = m;
	int incx = 1;
	int incy = 1;

	real *a = const_cast<real *>(A);
	real *x = const_cast<real *>(X);
	real *y = Y;

	return gemv(transa,m,n,alpha,a,lda,x,incx,beta,y,incy);		  
  }

  ///B = s*A + B
  ///\note parameters allocated outside and pointed to the device memory
  template<typename real>
  bool axpy(const real *A,real *B,long length,real s = 1.0f){
	
	if (A == NULL || B == NULL)
	  return false;

	int n = length;
	real *dx = const_cast<real *>(A);
	int incx = 1;
	real *dy = (B);
	int incy = 1;

	return axpy(n,s,dx,incx,dy,incy);
  }

  ///A *= s
  ///\note parameters allocated outside and pointed to the device memory
  template<typename real>
  bool scal(real *A,int n,real s){

	if (n == 0)
	  return true;

	if ( n < 0 || A == NULL){
	  return false;
	}
	real *dx = A;
	int incx = 1;
	return scal(n,s,dx,incx);
  }

}

#endif /* _WARPEDCUBLAS_H_ */
