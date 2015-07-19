#ifndef _MPRGPPRECONDITION_H_
#define _MPRGPPRECONDITION_H_

#include <iostream>
#include <stdio.h>
#include <vector>
#include <set>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "MPRGPUtility.h"
using namespace Eigen;
using namespace std;

namespace MATH{

  // no preconditioning
  template <typename T, typename MAT >
  class InFaceNoPreconSolver{

	typedef Eigen::Matrix<T,-1,1> Vec;

  public:
	InFaceNoPreconSolver(const vector<char>& face):_face(face){}
	int solve(const Vec&rhs, Vec&result, const Vec &phi){
	  MASK_FACE(rhs,result,_face);
	  return 0;
	}
	void setMatrix(const MAT& matrix){}

  protected:
	const vector<char>& _face;
  };

  // Jacobian preconditioner for boundary constraints
  template <typename T, typename MAT, bool NO_PRECOND=false >
  class DiagonalInFacePreconSolver{

	typedef Eigen::Matrix<T,-1,1> Vec;

  public:
	DiagonalInFacePreconSolver(const MAT &M,const vector<char> &face):face(face){

	  invert_diag.resize(M.rows());
	  if (!NO_PRECOND){
		for (int i = 0; i < invert_diag.size(); ++i){
		  assert_ge( M.diag(i), ScalarUtil<T>::scalar_eps );
		  invert_diag[i] = 1.0/M.diag(i);
		}
	  }
	}
	int solve(const Vec &g, Vec &z, const Vec &phi)const{

	  if (!NO_PRECOND){
		z.resize(phi.size());
		for(int i = 0; i < z.size(); i++){
		  assert_eq(invert_diag[i], invert_diag[i]);
		  z[i] = phi[i]*invert_diag[i];
		}
	  }else{
		z = phi;
	  }
	  return 0;
	}

  protected:
	Vec invert_diag;
	const vector<char> &face;
  };

  // Jacobian preconditioner for decoupled constraints
  template <typename T, typename MAT, bool NO_PRECOND=false >
  class DiagonalDecouplePreconSolver:public DiagonalInFacePreconSolver<T,MAT,NO_PRECOND>{

	typedef Eigen::Matrix<T,-1,1> Vec;

  public:
	DiagonalDecouplePreconSolver(const MAT &M,const vector<char> &face,const SparseMatrix<T> &J):
	  DiagonalInFacePreconSolver<T,MAT,NO_PRECOND>(M,face){

	  if (!NO_PRECOND && J.rows() > 0){
		JMinv = J*this->invert_diag.asDiagonal();
		const SparseMatrix<T> JMinvJt = JMinv*J.transpose();
		getDiagonal(JMinvJt, invert_JMinvJt);
		for (int i = 0; i < invert_JMinvJt.size(); ++i){
		  assert_ge( invert_JMinvJt[i], ScalarUtil<T>::scalar_eps );
		  invert_JMinvJt[i] = 1.0/invert_JMinvJt[i];
		}
	  }
	}
	int solve(const Vec &g, Vec &z, const Vec &phi)const{

	  if (!NO_PRECOND && JMinv.rows() > 0){
		Vec lambda = invert_JMinvJt.asDiagonal()*(JMinv*g);
		assert_eq(lambda.size(), this->face.size());
		for (size_t i = 0; i < this->face.size(); ++i){
		  if(0 == this->face[i])
			lambda[i] = 0;
		}
		z = this->invert_diag.asDiagonal()*g-JMinv.transpose()*lambda;
	  }else{
		z = phi;
	  }
	  return 0;
	}
	
  private:
	SparseMatrix<T> JMinv; // J*M^{-1}
	Vec invert_JMinvJt;
  };

}//end of namespace

#endif /* _MPRGPPRECONDITION_H_ */
