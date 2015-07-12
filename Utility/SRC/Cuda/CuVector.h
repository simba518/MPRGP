#ifndef _CUVECTOR_H_
#define _CUVECTOR_H_

/**
 * @file   CuVector.h
 * @author simba <simba@simba>
 * @date   Sat Jul  2 17:07:06 2011
 * 
 * @brief  an simple vector class that allocate memory from device memory.\n
 * used for cuda blas functions.
 * @see CuMatrix.h
 */

#include <stdio.h>
#include <CuFixedVec.h>

namespace LSW_CUDA_MATH{
  
  template <typename real>
  class CuVector:public CuFixedVec<real>{
	
  public:
	CuVector():CuFixedVec<real>(0,NULL){}

	CuVector(int len):CuFixedVec<real>(0,NULL){
	  real_resize(len);
	}

	CuVector(const CuVector<real> &other):CuFixedVec<real>(0,NULL){
	  
	  (*this) = other;
	}

	bool resize(int len){
	  return real_resize(len);
	}

	bool copyFromHost(const real *h_data,int len){
	  resize(len);
	  return CuFixedVec<real>::copyFromHost(h_data,len);
	}

	~CuVector(){

	  if( data != NULL){
		releaseDeviceMemory(data);
	  }	
	  length = 0;
	  data = NULL;
	}

	CuVector<real> & operator = (const CuVector<real> &other){

	  this->resize(other.size());
	  bool succ = ( cudaMemcpy(this->data,other.Data(),length*sizeof(real),cudaMemcpyDeviceToDevice) == cudaSuccess );
	  assert (succ);
	  return (*this);
	}
	
	CuVector<real> & operator = (const CuFixedVec<real> &other){

	  this->resize(other.size());
	  bool succ = ( cudaMemcpy(this->data,other.Data(),length*sizeof(real),cudaMemcpyDeviceToDevice) == cudaSuccess );
	  assert (succ);
	  return (*this);
	}

  protected:
	bool allocateDeviceMemory(const int length,real *&data)const{
	  bool succ = (cudaMalloc((void**)&data, length*sizeof(real)) == cudaSuccess);
	  if(!succ){
		data = NULL;
	  }
	  return succ;
	}
	bool releaseDeviceMemory(real *&data)const{
	  bool succ = true;
	  if(data != NULL){
		succ = (cudaFree(data) == cudaSuccess);
	  }
	  data = NULL;
	  return succ;
	}
	bool real_resize(int len){
	  
	  //check validation
	  assert (len>=0);
	  if(len < 0)
		return false;

	  //no need to change.
	  if(length == len){
		return true;
	  }
	  
	  //release memory allocated before
	  length = len;
	  bool succ = releaseDeviceMemory(data);
	  
	  //empty vector
	  if(length == 0){
		return succ;
	  }
	  
	  //allocate memory
	  succ = allocateDeviceMemory(length,data);
	  if( !succ ){
		length = 0;
		data = NULL;
	  }
	  return succ;
	}

  protected:
	using  CuFixedVec<real>::length;
	using  CuFixedVec<real>::data;
  };

  typedef CuVector<double> CuVectorXd;
  typedef CuVector<float> CuVectorXf;
  
}//end of namespace

#endif /*_CUVECTOR_H_*/
