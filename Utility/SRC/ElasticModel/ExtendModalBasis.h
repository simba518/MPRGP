#ifndef _EXTENDMODALBASIS_H_
#define _EXTENDMODALBASIS_H_

#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <MatrixTools.h>
using namespace Eigen;

namespace SIMULATOR{
  
  /**
   * @class ExtendModalBasis generate the extended modal basis.
   * @see An efficient construction of reduced deformable objects, siggraph asia, 2013.
   */
  class ExtendModalBasis{
	
  public:
	// construct basis for constrained object
	template<class T>
	static void construct(const Matrix<T,-1,-1> &linearBasis, Matrix<T,-1,-1> &extendedBasis){
	  extendNonRigidBasis(linearBasis, extendedBasis);
	  EIGEN3EXT::GramSchmidt(extendedBasis);
	}

	// construct basis for unconstrained object
	template<class T>
	static void construct(const Matrix<T,-1,-1> &linearBasis,const Matrix<T,-1,1> &rest_shape,Matrix<T,-1,-1> &extendedBasis){
	  
	  assert_gt(linearBasis.cols(), 6);
	  Matrix<T,-1,-1> non_rigid_B, rigid_B;
	  extendRigidBasis(rest_shape, rigid_B);
	  const Matrix<T,-1,-1> B1 = linearBasis.rightCols(linearBasis.cols()-6);
	  extendNonRigidBasis(B1, non_rigid_B);
	  extendedBasis.resize(non_rigid_B.rows(), non_rigid_B.cols()+rigid_B.cols());
	  extendedBasis.leftCols(rigid_B.cols()) = rigid_B;
	  extendedBasis.rightCols(non_rigid_B.cols()) = non_rigid_B;
	  EIGEN3EXT::GramSchmidt(extendedBasis);
	}

	template<class T>
	static void extendRigidBasis(const Matrix<T,-1,1> &rest_shape, Matrix<T,-1,-1> &extendedBasis){

	  extendedBasis.resize(rest_shape.size(),12);
	  extendedBasis.setZero();
	  for (int i = 0; i < 3; ++i)
		for (int j = 0;  j < rest_shape.size()/3; ++j)
		  extendedBasis(j*3+i,i) = 1.0f;

	  Matrix<T,-1,-1> B;
	  extendNonRigidBasis(rest_shape,B);
	  extendedBasis.rightCols(9) = B;
	}

	template<class T, int cols>
	static void extendNonRigidBasis(const Matrix<T,-1,cols> &linearBasis, Matrix<T,-1,-1> &extendedBasis){

	  const int r = linearBasis.cols();
	  const int n = linearBasis.rows()/3;
	  assert_eq(n*3, linearBasis.rows());

	  extendedBasis.resize(n*3, 9*r);
	  extendedBasis.setZero();

	  for (int i = 0; i < r; ++i)
		for (int j = 0; j < n; ++j)
		  for (int k = 0; k < 3; ++k)
			for (int s = 0; s < 3; ++s)
			  extendedBasis(j*3+k,i*9+k*3+s) = linearBasis(j*3+s,i);
	}

  };
  
}//end of namespace

#endif /*_EXTENDMODALBASIS_H_*/
