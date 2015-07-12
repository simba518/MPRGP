#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <SparseGenEigenSolver.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(SparseGenEigenSolverTest)

/**
 * K:
 * 1.0  0.1  0.0 
 * 0.1  2.0  0.0
 * 0.0  0.0  3.0
 *
 * M:
 * 2.3  0.0  0.0 
 * 0.0  1.2  0.0
 * 0.0  0.0  3.9
 */
BOOST_AUTO_TEST_CASE(TestEigenSolver){
  
  const int n = 3;
  Eigen::SparseMatrix<double> K(n,n);
  K.insert (0,0) = 1.0f;
  K.insert (1,0) = 0.10f;
  K.insert (1,1) = 2.0f;
  K.insert (2,2) = 3.0f;

  Eigen::DiagonalMatrix<double,-1> M(n);
  M.diagonal()[0] = 2.3;
  M.diagonal()[1] = 1.2;
  M.diagonal()[2] = 3.9;

  MatrixXd eig_vec;
  VectorXd eig_val;
  const int max_eig_num = 2;
  const MatrixXd Km = K;
  const MatrixXd Mm = M;
  Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd > es(Km,Mm);
  const MatrixXd eigen_m = es.eigenvectors().topLeftCorner(n,max_eig_num);

  K.makeCompressed();
  TEST_ASSERT( EigenSparseGenEigenSolver::solve(K,M,eig_vec,eig_val,max_eig_num) );
  ASSERT_EQ_SMALL_VEC_TOL( es.eigenvalues(), eig_val, eig_val.size(), 1e-8 );
  for (int i = 0; i < max_eig_num; ++i){
	TEST_ASSERT( ((eigen_m.col(i)-eig_vec.col(i)).norm() < 1e-8) || 
					((eigen_m.col(i)+eig_vec.col(i)).norm() < 1e-8) );
  }

  TEST_ASSERT( EigenSparseGenEigenSolver::solve(K,M,eig_vec,eig_val,max_eig_num,'R',"LM") );
  ASSERT_EQ_SMALL_VEC_TOL( es.eigenvalues().tail(max_eig_num), eig_val, eig_val.size(),1e-8);

  TEST_ASSERT( EigenSparseGenEigenSolver::solve(K,M,eig_vec,eig_val,2,'R',"BE") );
  ASSERT_EQ_TOL(es.eigenvalues()[0], eig_val[0], 1e-8);
  ASSERT_EQ_TOL(es.eigenvalues()[K.rows()-1], eig_val[1], 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
