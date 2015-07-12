#ifndef _CUVECTOREXT_H_
#define _CUVECTOREXT_H_

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense>
using Eigen::Matrix;
using Eigen::Dynamic;

#include <CuVector.h>

namespace LSW_CUDA_MATH{
  
  /**
   * @class CuVectorExt extension of the CuVector class.
   * 
   */
  template <typename real>
  class CuVectorExt: public CuVector<real>{
	
  public:
	CuVectorExt(){
	  
	}

	CuVectorExt(int len):CuVector<real>(len){
	  
	}

	CuVectorExt(const Matrix< real, Dynamic, 1> &host_v):CuVector<real>(0){
	  (*this) = host_v;
	}

	CuVectorExt<real>& operator = (const Matrix< real, Dynamic, 1> &host_v){

	  if(host_v.size() > 0){
		CuVector<real>::copyFromHost(&host_v(0),host_v.size());
	  }else{
		CuVector<real>::resize(0);
	  }	
	  return *this;
	}

	bool copyTo(Matrix< real, Dynamic, 1> &host_v)const{

	  host_v.resize(this->size());
	  if(host_v.size() > 0){
		return copyToHost(&host_v(0));
	  }
	  return true;
	}	
  };

  typedef CuVectorExt<double> CuVectorXdE;
  typedef CuVectorExt<float> CuVectorXfE;

  
}//end of namespace

#endif /*CuVectorExt*/
