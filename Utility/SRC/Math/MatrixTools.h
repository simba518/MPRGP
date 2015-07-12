#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <assertext.h>
#include <boost/foreach.hpp>

namespace EIGEN3EXT{

  /************************************ create ********************************/
  template <class T>
  inline const Eigen::Matrix<T,3,3> &createFromRowMajor (Eigen::Matrix<T,3,3> &M, const T *A){

	const int r = 3;
	const int c = 3;
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		M(i,j) = A[i*c + j];
	  }
	}
	return M;
  }


  template <class T>
  inline const Eigen::Matrix<T,-1,-1> &createFromRowMajor(Eigen::Matrix<T,-1,-1> &M,
														  const T *A, const int r, const int c){

	assert_ge (r,0);
	assert_ge (c,0);
	M.resize(r,c);
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		M(i,j) = A[i*c + j];
	  }
	}
	return M;
  }

  template <class T, class MATRIX>
  inline void createRowMajor(const MATRIX &M,T *A){

	const int r = M.rows();
	const int c = M.cols();
	for (int i = 0; i < r; ++i) {
	  for (int j = 0; j < c; ++j) {
		A[i*c + j] = M(i,j);
	  }
	}
  }

  template <class T>
  inline const Eigen::Matrix<T,-1,-1> &createFromColMajor(Eigen::Matrix<T,-1,-1> &M, 
														  const T *A, const int r, const int c){
  	assert_ge (r,0);
  	assert_ge (c,0);
  	M.resize(r,c);
  	for (int i = 0; i < r; ++i) {
  	  for (int j = 0; j < c; ++j) {
  		M(i,j) = A[j*r + i];
  	  }
  	}
  	return M;
  }

  template <class T, class MATRIX> 
  inline void convert(const MATRIX &M, std::vector<Eigen::Matrix<T,-1,1> > &vv){
	vv.clear();
	vv.reserve(M.cols());
	for (int c = 0; c < M.cols(); ++c){
	  vv.push_back(M.col(c));
	}
  }

  template <class T> 
  inline void convert(const std::vector<Eigen::Matrix<T,-1,1> > &vv,Eigen::Matrix<T,-1,-1> &M){

	M.resize(0,0);
	if(vv.size() > 0){
	  const int cols = (int)vv.size();
	  const int rows = vv[0].size();
	  if(cols > 0){
		M.resize(rows,cols);
		for (int c = 0; c < M.cols(); ++c){
		  assert_eq(vv[c].size(),rows);
		  M.col(c) = vv[c];
		}
	  }
	}
  }

  /************************************ decomposition *************************/
  
  /*
   * Modified Sigular Value Decomposition.
   * F = U*D*Vt, where U and Vt are pure rotation matrices, e.g. det(U)=det(V)=1.
   * @see study record: Polar Decomposition.
   */
  template <class T>
  inline void ModifiedSVD3x3(const Eigen::Matrix<T,3,3> &F,
							 Eigen::Matrix<T,3,3> &U,
							 Eigen::Matrix<T,3,3> &Vt,
							 Eigen::Matrix<T,3,3> &D){

	Eigen::JacobiSVD<Eigen::Matrix<T,3,3> > svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = svd.matrixU();
	Vt = svd.matrixV().transpose();
	assert_eq (svd.singularValues().size(),3);
	D.setZero();
	D(0,0) = svd.singularValues()[0];
	D(1,1) = svd.singularValues()[1];
	D(2,2) = svd.singularValues()[2];

	if (U.determinant() < 0) {
	  U.col(2) *= (-1.0f);
	  D(2,2) *= -1.0f;
	}
	if (Vt.determinant() < 0){
	  Vt.row(2) *= (-1.0f);
	  D(2,2) *= -1.0f;
	}
  }

  /*
   * Modified Polar Decomposition.
   * F = RS, where R is a pure rotation matrix, e.g. det(R)=1, and S is symetric
   * but maybe not positive-definit.
   * @see study record: Polar Decomposition.
   */
  template <class T>
  inline void ModifiedPD3x3(const Eigen::Matrix<T,3,3> &F,
							Eigen::Matrix<T,3,3> &R,
							Eigen::Matrix<T,3,3> &S){

  	Eigen::Matrix<T,3,3> &U = R;
	Eigen::Matrix<T,3,3> &Vt = S;
	Eigen::Matrix<T,3,3> D;
  	ModifiedSVD3x3(F,U,Vt,D);
  	R = U*Vt;
  	S = Vt.transpose()*D*Vt;
  }


  // comput the pseudoinverse of a dense matrix using SVD.
  // https://inst.eecs.berkeley.edu/~ee127a/book/login/def_pseudo_inv.html
  template <class T>
  inline void PseudoInverse(const Eigen::Matrix<T,-1,-1> &A,Eigen::Matrix<T,-1,-1> &invA,const T tol=0.0){
	Eigen::JacobiSVD <Eigen::Matrix<T,-1,-1> > svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);
	const Eigen::Matrix<T,-1,-1> &U = svd.matrixU();
	const Eigen::Matrix<T,-1,-1> &V = svd.matrixV();
	const Eigen::Matrix<T,-1,1>  &s = svd.singularValues();
	Eigen::Matrix<T,-1,1> invS = s;
	for (int i = 0; i < invS.size(); ++i){
	  if(abs(invS[i]) > tol)
		invS[i] = 1.0f/invS[i];
	}
	invA = V*invS.asDiagonal()*U.transpose();
  }

  template <class T>
  inline Eigen::Matrix<T,-1,-1> PseudoInverse(const Eigen::Matrix<T,-1,-1> &A, const T tol=0.0){
	Eigen::Matrix<T,-1,-1> invA;
	PseudoInverse(A,invA,tol);
	return invA;
  }

  template<class MATRIX_TMP, class T>
  inline void MGramSchmidt(const MATRIX_TMP &M,Eigen::Matrix<T,-1,-1> &U){

	assert_eq(M.rows(),M.cols());
	const int n = U.cols();
	for(int i=0;i<n;i++){
	  for(int j=0;j<i;j++){
		const T a = (U.col(j).dot(M*U.col(j)));
		assert_ne(a,0.0f);
		assert_eq(a,a);
		const T alpha=(U.col(i).dot(M*U.col(j)))/a;
		U.col(i)-=alpha*U.col(j);
	  }
	  U.col(i)/=sqrt((U.col(i).dot(M*U.col(i))));
	}
  }

  template<class T>
  inline void GramSchmidt(Eigen::Matrix<T,-1,-1> &U){

	const int n = U.cols();
	for(int i=0;i<n;i++){
	  for(int j=0;j<i;j++){
		const T a = (U.col(j).dot(U.col(j)));
		assert_ne(a,0.0f);
		assert_eq(a,a);
		const T alpha=(U.col(i).dot(U.col(j)))/a;
		U.col(i)-=alpha*U.col(j);
	  }
	  U.col(i)/=sqrt((U.col(i).dot(U.col(i))));
	}
  }
  
}

#endif /* _MATRIXTOOLS_H_ */
