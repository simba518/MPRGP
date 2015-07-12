#ifndef _CUMATRIX_H_
#define _CUMATRIX_H_

/**
 * @file   CuMatrix.h
 * @author simba <simba@simba>
 * @date   Sat Jul  2 16:43:56 2011
 * 
 * @brief  this is an simple matrix class that allocate memory from device memory.\n
 * used for cuda blas functions.
 * @see CuVector.h
 */

#include <cuda_runtime_api.h>
#include <stdio.h>
#include <CuFixedVec.h>
#include <CuFixedMat.h>

namespace LSW_CUDA_MATH{
  
  template <typename real>
  class CuMatrix:public CuFixedMat<real>{
  public:
	CuMatrix():CuFixedMat<real>(0,0,NULL){}

	CuMatrix(int _r,int _c):CuFixedMat<real>(0,0,NULL){
	  resize(_r,_c);
	}

	CuMatrix(const CuMatrix<real> &other):CuFixedMat<real>(0,0,NULL){
	  (*this) = other;
	}

	bool resize(int i,int j){
	  
	  //check parameters
	  assert (i >=0 && j>=0);
	  if(i < 0 && j < 0)
		return false;
	  
	  //no need to allocate new memory
	  if(r*c == i*j){
		r = i;
		c = j;
		return true;
	  }

	  //release memory allocated before
	  r = i;
	  c = j;
	  bool succ = releaseDeviceMemory(data);

	  //empty matrix
	  if(r*c == 0){
		r = 0;
		c = 0;
		return succ;
	  }

	  //allocate memory
	  succ &= allocateDeviceMemory(r*c,data);
	  if( !succ ){
		r = 0;
		c = 0;
		data = NULL;
	  }
	  return succ;
	}

	bool copyFromHost(const real *h_data,int _r,int _c){
	  if( h_data != NULL && resize(_r,_c) ){
		return ( cudaMemcpy(data,h_data,_r*_c*sizeof(real),cudaMemcpyHostToDevice) == cudaSuccess );
	  }
	  return false;
	}

	~CuMatrix(){
	  
	  if( data != NULL){
		releaseDeviceMemory(data);
	  }
	  data = NULL;
	  r = 0;
	  c = 0;
	}
	
	CuMatrix<real> & operator = (const CuMatrix<real> &other){
	  
	  this->resize(other.rows(),other.cols());
	  int length = this->rows()*this->cols();
	  bool succ = ( cudaMemcpy(this->data,other.Data(),length*sizeof(real),cudaMemcpyDeviceToDevice) == cudaSuccess );
	  assert(succ);
	  return *this;
	}

  protected:
	bool allocateDeviceMemory(const int length,real *&data)const{
	  return (cudaMalloc((void**)&data, length*sizeof(real)) == cudaSuccess);
	}

	bool releaseDeviceMemory(real *&data)const{
	  bool succ = true;
	  if(data != NULL){
		succ = (cudaFree(data) == cudaSuccess);
	  }
	  data = NULL;
	  return succ;
	}

  protected:
	using CuFixedMat<real>:: r;
	using CuFixedMat<real>:: c;
	using CuFixedMat<real>:: data;

  };

  typedef CuMatrix<double> CuMatrixXd;
  typedef CuMatrix<float> CuMatrixXf;
  
}//end of namespace

#endif /*_CUMATRIX_H_*/
