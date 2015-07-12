#ifndef _WARPEDCUBLASFUN_H_
#define _WARPEDCUBLASFUN_H_

#include "CuMatrix.h"
#include "CuVector.h"
#include "WarpedCuBlas.h"
#include "CuFixedVec.h"

namespace LSW_CUDA_MATH{

  ///C := alpha*op( A )*op( B ) + beta*C	
  template <typename real >
  bool M_x_M(const CuFixedMat<real> &A,const CuFixedMat<real> &B,CuMatrix<real> &C,
			 char transa = 'n',char transb = 'n',real alpha = 1.0f,real beta = 0.0f){
	
	//empty matrix
	if (A.cols()*A.rows() == 0){
	  C.resize(0,0);
	  return (B.cols()*B.rows() == 0);
	}
	//resize C
	int m = (transa == 'n' ? A.rows() : A.cols());//rows of op(A)
	int n = (transb == 'n' ? B.cols() :B.rows());//columns of op(B)
	C.resize(m,n);

	return gemm(A.Data(),A.rows(),A.cols(),
				B.Data(),B.rows(),B.cols(),
				C.Data(),transa,transb,alpha,beta);

  }

  template <typename real >
  bool M_x_M(const CuFixedMat<real> &A,const CuFixedMat<real> &B,CuFixedMat<real> &C,
			 char transa = 'n',char transb = 'n',real alpha = 1.0f,real beta = 0.0f){
	
	//empty matrix
	if (A.cols()*A.rows() == 0){
	  return (B.cols()*B.rows() == 0) &&  (C.cols()*C.rows() == 0);
	}

	//check size of C
	int m = (transa == 'n' ? A.rows() : A.cols());//rows of op(A)
	int n = (transb == 'n' ? B.cols() :B.rows());//columns of op(B)
	if(C.rows() != m || C.cols() != n){
	  return false;
	}
	return gemm(A.Data(),A.rows(),A.cols(),
				B.Data(),B.rows(),B.cols(),
				C.Data(),transa,transb,alpha,beta);
  }

  ///y = alpha*op(A)*x + beta*y
  template <typename real>
  bool M_x_V(const CuFixedMat<real> &A,const CuFixedVec<real> &X,CuVector<real> &Y,
			 char transa = 'n',real alpha = 1.0f,real beta = 0.0f){
	
	//resize Y
	int  len_c = ( transa == 'n' ? A.rows() : A.cols() );
	if ( len_c  <=0 ){
	  Y.resize(0);
	  return true;
	}
	if (beta != 0.0f){
	  if ( Y.size() != len_c ){
		return false;
	  }
	}else{
	  Y.resize(len_c);
	}
	return gemv(A.Data(),A.rows(),A.cols(),X.Data(),Y.Data(),transa,alpha,beta);
  }

  template <typename real>
  bool M_x_V(const CuFixedMat<real> &A,const CuFixedVec<real> &X,CuFixedVec<real> &Y,
			 char transa = 'n',real alpha = 1.0f,real beta = 0.0f){
	
	int  len_c = ( transa == 'n' ? A.rows() : A.cols() );
	if ( len_c  <=0 ){
	  return true;
	}
	if ( Y.size() != len_c ){//Y can't be resized
	  return false;
	}else{
	  return gemv(A.Data(),A.rows(),A.cols(),X.Data(),Y.Data(),transa,alpha,beta);
	}
  }
}

#endif /* _WARPEDCUBLASFUN_H_ */
