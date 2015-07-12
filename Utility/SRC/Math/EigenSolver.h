#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <assertext.h>
#include <Log.h>

namespace EIGEN3EXT{
  
  // solve the largest eigenvalue in absolute value using power method.
  // @see http://ergodic.ugr.es/cphys/LECCIONES/FORTRAN/power_method.pdf
  // @note this mehtod maynot convergent unless A is diagonalizable and have a dominant eigenvalue.
  template<class VECTOR, class SCALOR, class SOLVER>
  int powerMethod(SOLVER &solver, VECTOR &eig_vec, const int dim,const int max_it=0,const SCALOR tol=1e-8){

	assert_ge(max_it,0);
	assert_ge(tol,0);
	assert_gt(dim,0);

	VECTOR x(dim);
	for (int i = 0; i<dim; ++i) x[i] = 1.0;
	x *= (1.0/x.norm());

	int it = 0;
	for (; it<max_it || max_it<=0; ++it){
	  eig_vec = solver.solve(x);
	  eig_vec *= (1.0/eig_vec.norm());
	  if ( (eig_vec-x).norm() <= tol )
		break;
	  x = eig_vec;
	}
	return it;
  }

  template<class MATRIX>
  class MatrixMultplyVector{
	
  public:
    MatrixMultplyVector(){ A = NULL; }
	MatrixMultplyVector(const MATRIX &A){ compute(A); }
	void compute(const MATRIX &A){
	  this->A = &A;
	}
	template<class VECTOR>
	VECTOR solve(const VECTOR &b){
	  assert(A);
	  return (*A)*b;
	}
	
  private:
	const MATRIX *A;
  };

  // linear solver for general eigen value problem K*x = lambda*M*x.
  template<class SOLVER_K, class SOLVER_M>
  class GeneralProblemSolver{
	
  public:
	GeneralProblemSolver(SOLVER_K&Ks,SOLVER_M&Ms):K_solver(Ks),M_solver(Ms){}
	template<class VECTOR>
	VECTOR solve(const VECTOR &b){
	  return M_solver.solve(K_solver.solve(b));
	}
	
  private:
	SOLVER_K &K_solver;
	SOLVER_M &M_solver;
  };
  
  template<class MATRIX, class VECTOR>
  double eigenValue(const MATRIX &A, const VECTOR &eig_vec){
	assert_eq(A.rows(), A.cols());
	assert_eq(A.rows(), eig_vec.size());
	return eig_vec.dot(A*eig_vec)/(eig_vec.dot(eig_vec));
  }

  // K*x = lambda*M*x.
  template<class MATRIX_K, class MATRIX_M, class VECTOR>
  double eigenValue(const MATRIX_K &K, const MATRIX_M &M, const VECTOR &eig_vec){

	assert_eq(K.rows(), K.cols());
	assert_eq(M.rows(), M.cols());
	assert_eq(M.rows(), K.cols());
	assert_eq(K.rows(), eig_vec.size());

	return eig_vec.dot(K*eig_vec)/(eig_vec.dot(M*eig_vec));
  }

}//end of namespace

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/UmfPackSupport>
using namespace Eigen;

namespace EIGEN3EXT{

  // solve the largest/smallest |lambda| for K*x = lambda*x.
  template<class T>
  int largestEigen(const SparseMatrix<T> &K,
					  Matrix<T,-1,1>& eig_vec, T &eig_value,
					  const int max_it=0,const T tol=1e-8){

	assert_eq(K.cols(), K.rows());
	MatrixMultplyVector<SparseMatrix<double> > K_s(K);
	const int it = powerMethod(K_s,eig_vec,K.rows(),max_it,tol);
	eig_value = eigenValue(K,eig_vec);
	return it;
  }

  template<class T>
  int smallestEigen(const SparseMatrix<T> &K,
					Matrix<T,-1,1>& eig_vec, T &eig_value,
					const int max_it=0,const T tol=1e-8){

	assert_eq(K.cols(), K.rows());
	UmfPackLU<SparseMatrix<double> > K_s(K);
	if(K_s.info()!=Eigen::Success) {
	  ERROR_LOG("UmfPackLU decomposition failed");
	  return -1;
	}

	const int it = powerMethod(K_s,eig_vec,K.rows(),max_it,tol);
	eig_value = eigenValue(K,eig_vec);
	return it;
  }

  template<class T>
  int smallestEigenSPD(const SparseMatrix<T> &K,
					   Matrix<T,-1,1>& eig_vec, T &eig_value,
					   const int max_it=0,const T tol=1e-8){

	assert_eq(K.cols(), K.rows());
	CholmodSupernodalLLT<SparseMatrix<double> > K_s(K);
	if(K_s.info()!=Eigen::Success) {
	  ERROR_LOG("CholmodSupernodalLLT decomposition failed");
	  return -1;
	}

	const int it = powerMethod(K_s,eig_vec,K.rows(),max_it,tol);
	eig_value = eigenValue(K,eig_vec);
	return it;
  }

  // solve the largest/smallest |lambda| for K*x = lambda*M*x, where K is sym, and M is SPD.
  template<class T>
  int largestGenEigenSym(const SparseMatrix<T> &K,const SparseMatrix<T> &M,
						 Matrix<T,-1,1>& eig_vec, T &eig_value,
						 const int max_it=0,const T tol=1e-8){

	assert_eq(K.cols(), K.rows());
	assert_eq(M.cols(), M.rows());
	assert_eq(K.rows(), M.rows());

	typedef MatrixMultplyVector<SparseMatrix<double> > MV_Solver;
	// typedef CholmodSupernodalLLT<SparseMatrix<double> > LLT_Solver;
	typedef UmfPackLU<SparseMatrix<double> > LU_Solver;

	MV_Solver K_s(K);
	LU_Solver M_s(M);
	if(M_s.info()!=Eigen::Success) {
	  ERROR_LOG("UmfPackLU decomposition failed");
	  return -1;
	}

	GeneralProblemSolver<MV_Solver,LU_Solver> solver(K_s,M_s);
	const int it = powerMethod(solver,eig_vec,K.rows(),max_it,tol);
	eig_value = eigenValue(K,M,eig_vec);

	return it;
  }

  template<class T>
  int smallestGenEigenSym(const SparseMatrix<T> &K,const SparseMatrix<T> &M,
						  Matrix<T,-1,1>& eig_vec, T &eig_value, 
						  const int max_it=0,const T tol=1e-8){

	assert_eq(K.cols(), K.rows());
	assert_eq(M.cols(), M.rows());
	assert_eq(K.rows(), M.rows());

	typedef MatrixMultplyVector<SparseMatrix<double> > MV_Solver;
	// typedef CholmodSupernodalLLT<SparseMatrix<double> > LLT_Solver;
	typedef UmfPackLU<SparseMatrix<double> > LU_Solver;

	MV_Solver K_s(M);
	LU_Solver M_s(K);
	if(M_s.info()!=Eigen::Success) {
	  ERROR_LOG("UmfPackLU decomposition failed");
	  return -1;
	}

	GeneralProblemSolver<MV_Solver,LU_Solver> solver(K_s,M_s);
	const int it = powerMethod(solver,eig_vec,K.rows(),max_it,tol);
	eig_value = eigenValue(K,M,eig_vec);

	return it;
  }


}//end of namespace

#endif /*_EIGENSOLVER_H_*/
