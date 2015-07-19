#ifndef _TEST_UTILITY_H_
#define _TEST_UTILITY_H_

#include "MPRGPUtility.h"
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace MATH;

typedef Matrix<double,4,1> Vec4d;
typedef vector<Vec4d,aligned_allocator<Vec4d> > VVec4d;
typedef vector<VVec4d > VVVec4d;

template <class T>
const SparseMatrix<T> &createFromDense(const Matrix<T,-1,-1> &M, SparseMatrix<T> &S, const T tol=1e-16){
  
  typedef Triplet<T> E_Triplet;
  std::vector<E_Triplet> striplet;
  striplet.reserve(M.size());
  for (int i = 0; i < M.rows(); ++i) {
	for (int j = 0; j < M.cols(); ++j) {
	  if ( fabs(M(i,j)) >= tol )
		striplet.push_back( E_Triplet(i,j,M(i,j)) );
	}
  }
  S.resize(M.rows(), M.cols());
  S.setFromTriplets(striplet.begin(), striplet.end());
  return S;
}

template <class T> 
const SparseMatrix<T> createFromDense(const Matrix<T,-1,-1> &M, const T tol=1e-16){
  SparseMatrix<T> S;
  return createFromDense(M,S,tol);
}

#endif /* _TEST_UTILITY_H_ */
