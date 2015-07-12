#include "SparseGenEigenSolver.h"
#include <memory.h>
#include <stdlib.h>
#include <arpack++/arlsmat.h>
#include <arpack++/arlgsym.h>

template <typename T>
bool EIGEN3EXT::eigensym( const T *valA,int nnzA,  const int *irowA,  const int *pcolA,
						 const T *valB,int nnzB,  const int *irowB,  const int *pcolB,
						  int n,int n_eig,T *vals,T *vecs,
						  const char mode,
						  const std::string which){
  
  //check parameters
  const bool valid = (valA!=NULL&&irowA!=NULL&&pcolA!=NULL&&
				valB!=NULL&&irowB!=NULL&&pcolB!=NULL&&
				n_eig>0&&n_eig<n&&vals!=NULL&&vecs!=NULL);
  if(!valid){
	return false;
  }

  ARluSymMatrix<T> A(n, nnzA, const_cast<T*>(valA), const_cast<int*>(irowA), const_cast<int*>(pcolA));
  ARluSymMatrix<T> B(n, nnzB, const_cast<T*>(valB), const_cast<int*>(irowB), const_cast<int*>(pcolB));

  if ('R' == mode){
	ARluSymGenEig<T> dprob(n_eig, A, B , const_cast<char*>(which.c_str()));
	int coveraged = dprob.EigenValVectors(vecs,vals);
	return (coveraged==n_eig);
  }else{
	ARluSymGenEig<T> dprob(mode,n_eig, A, B,0.0f,const_cast<char*>(which.c_str()));
	int coveraged = dprob.EigenValVectors(vecs,vals);
	return (coveraged==n_eig);
  }
}

template bool EIGEN3EXT::eigensym( const double *valA,int nnzA,  const int *irowA,  const int *pcolA,
								   const double *valB,int nnzB,  const int *irowB,  const int *pcolB,
								   int n,int n_eig,double *vals,double *vecs,
								   const char mode,
								   const std::string which);

template bool EIGEN3EXT::eigensym( const float *valA,int nnzA,  const int *irowA,  const int *pcolA,
								   const float *valB,int nnzB,  const int *irowB,  const int *pcolB,
								   int n,int n_eig,float *vals,float *vecs,
								   const char mode,
								   const std::string which);
