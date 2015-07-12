#ifndef _CUDIAMAT_H_
#define _CUDIAMAT_H_

#include <cuda_runtime_api.h>
#include <stdio.h>
#include "WarpedCuBlas.h"
#include <CuVector.h>
#include <EleProduct.h>

namespace LSW_CUDA_MATH{
  
  template <typename real>
  class CuDiaMat{
	
  public:
	CuDiaMat(){
	  dim = 0;
	  data = NULL;
	}
	CuDiaMat(int _dim){
	  dim = 0;
	  data = NULL;
	  resize(_dim);
	}
	bool resize(int _dim){

	  //check validation
	  assert (_dim>=0);
	  if(_dim < 0)
		return false;

	  //no need to change.
	  if(dim == _dim){
		return true;
	  }
	  
	  //release memory allocated before
	  dim = _dim;
	  bool succ = releaseDeviceMemory(data);
	  
	  //empty vector
	  if(dim == 0){
		return succ;
	  }
	  
	  //allocate memory
	  succ = allocateDeviceMemory(dim,data);
	  if( !succ ){
		dim = 0;
		data = NULL;
	  }
	  return succ;
	}
	const real *Data()const{
	  return data;
	}
	real *Data(){
	  return data;
	}
	
	int Dim()const{
	  return dim;
	}

	bool copyFromHost(const real *h_data,int _dim){
	  if( h_data != NULL && resize(_dim) ){
		return ( cudaMemcpy(data,h_data,_dim*sizeof(real),cudaMemcpyHostToDevice) == cudaSuccess );
	  }
	  return false;
	}
	
	bool copyToHost(real *h_data)const{
	  if( h_data != NULL && this->Data() != NULL){
		return ( cudaMemcpy(h_data,data,dim*sizeof(real),cudaMemcpyDeviceToHost) == cudaSuccess );
	  }
	  return false;
	}

	bool dot(const CuFixedVec<real> v, CuFixedVec<real> result)const{
	  assert (v.size() == result.size() && v.size()==Dim() );
	  return EleProduct::getInstance()->product(v.Data(),Data(),result.Data(),Dim());
	}

	bool dot(const CuFixedVec<real> v, CuVector<real> &result)const{
	  if(result.resize(dim)){
		return this->dot(v,(CuFixedVec<real> &)result);
	  }
	  return false;
	}

	~CuDiaMat(){

	  if( data != NULL){
		releaseDeviceMemory(data);
	  }	
	  dim = 0;
	  data = NULL;
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
	int dim;
	real *data;
  };
  
}//end of namespace

#endif /*_CUDIAMAT_H_*/
