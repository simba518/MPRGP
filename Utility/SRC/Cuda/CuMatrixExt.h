#ifndef _CUMATRIXEXT_H_
#define _CUMATRIXEXT_H_

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense>
using Eigen::Matrix;
using Eigen::Dynamic;

#include <CuMatrix.h>

namespace LSW_CUDA_MATH{
  
  /**
   * @class CuMatrixExt extension of the CuMatrix class.
   * 
   */
  template <typename real>
  class CuMatrixExt: public CuMatrix<real>{
	
  public:
	CuMatrixExt(){}

	CuMatrixExt(int _r,int _c):CuMatrix<real>(_r,_c){
	  
	}

	CuMatrixExt(const Matrix<real, Dynamic, Dynamic> &host_m):CuMatrix<real>(0,0){
	  (*this) = host_m;
	}

	CuMatrixExt<real>& operator = (const Matrix<real, Dynamic, Dynamic> &host_m){
	  
	  if(host_m.rows() >0 && host_m.cols() > 0){

		CuMatrix<real>::copyFromHost(&host_m(0,0),host_m.rows(),host_m.cols());
	  }else{
		CuMatrix<real>::resize(0,0);
	  }
	  return (*this);
	}

	bool copyTo(Matrix<real, Dynamic, Dynamic> &host_m)const{

	  host_m.resize(this->rows(),this->cols());
	  if(host_m.size() > 0){
		return copyToHost(&host_m(0,0));
	  }
	  return true;
	}

  };
  
  typedef CuMatrixExt<double> CuMatrixXdE;
  typedef CuMatrixExt<float> CuMatrixXfE;
  
}//end of namespace

#endif /*CuMatrixExt*/
