#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <EigenSolver.h>
#include <MatrixIO.h>
#include <SparseMatrixTools.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(EigenSolverTest)

BOOST_AUTO_TEST_CASE(testLargestEigenValue){
  
  MatrixXd A(3,3);
  A <<  1, 2, 0,
	   -2, 1, 2,
	    1, 3, 1;

  VectorXd eig;
  MatrixMultplyVector<MatrixXd> solver(A);
  powerMethod(solver,eig,A.rows(),1000,1e-8);
  const double lambda = eigenValue(A,eig);

  VectorXd correct_eig(3);
  correct_eig << 0.408248289878, 0.408248289878, 0.816496581513;
  ASSERT_EQ_SMALL_VEC_TOL(correct_eig, eig, eig.size(), 1e-10);
  ASSERT_EQ_TOL(lambda, 2.99999999857, 1e-10);

}

BOOST_AUTO_TEST_CASE(testSmallestEigenValue){

  MatrixXd A(2,2);
  A <<  13, -4,
	    -4, 7;

  VectorXd eig;
  FullPivLU<MatrixXd> solver(A);
  powerMethod(solver,eig,A.rows(),1000,1e-8);
  const double lambda = eigenValue(A,eig);

  VectorXd correct_eig(2);
  correct_eig << 0.447213597809, 0.894427189846;
  ASSERT_EQ_SMALL_VEC_TOL(correct_eig, eig, eig.size(), 1e-10);
  ASSERT_EQ_TOL(lambda, 5, 1e-10);
}

BOOST_AUTO_TEST_CASE(testLargestGeneralEigenValue){

  MatrixXd K(3,3);
  K <<  1, 2, 0,
	   -2, 1, 2,
	    1, 3, 1;

  MatrixXd M(3,3);
  M <<  2, 0, 0,
	    0, 2, 0,
        0, 0, 2;

  VectorXd eig;
  MatrixMultplyVector<MatrixXd> K_solver(K);
  LLT<MatrixXd> M_solver(M);
  GeneralProblemSolver<MatrixMultplyVector<MatrixXd>,LLT<MatrixXd> >solver(K_solver,M_solver);
  powerMethod(solver,eig,K.rows(),1000,1e-8);
  const double lambda = eigenValue(K,M,eig);

  VectorXd correct_eig(3);
  correct_eig << 0.408248289878, 0.408248289878, 0.816496581513;
  ASSERT_EQ_SMALL_VEC_TOL(correct_eig, eig, eig.size(), 1e-10);
  ASSERT_EQ_TOL(lambda, 1.5, 1e-7);

  const SparseMatrix<double> Ks = createFromDense(K);
  const SparseMatrix<double> Ms = createFromDense(M);
  VectorXd eig_vec;
  double eig_value = 1.0f;
  largestGenEigenSym(Ks,Ms,eig_vec,eig_value,1000,1e-8);
  ASSERT_EQ_SMALL_VEC_TOL(eig_vec, correct_eig, eig_vec.size(), 1e-10);
  ASSERT_EQ_TOL(eig_value, 1.5, 1e-7);

}

BOOST_AUTO_TEST_CASE(testSmallestGeneralEigenValue){

  MatrixXd K(2,2);
  K <<  13, -4,
	    -4, 7;

  MatrixXd M(2,2);
  M <<  2, 0,
	    0, 2;

  VectorXd eig;
  MatrixMultplyVector<MatrixXd> K_solver(M);
  LLT<MatrixXd> M_solver(K);
  GeneralProblemSolver<MatrixMultplyVector<MatrixXd>,LLT<MatrixXd> >solver(K_solver,M_solver);
  powerMethod(solver,eig,K.rows(),1000,1e-8);
  const double lambda = eigenValue(K,M,eig);

  VectorXd correct_eig(2);
  correct_eig << 0.447213597809, 0.894427189846;
  ASSERT_EQ_SMALL_VEC_TOL(correct_eig, eig, eig.size(), 1e-10);
  ASSERT_EQ_TOL(lambda, 2.5, 1e-7);

  const SparseMatrix<double> Ks = createFromDense(K);
  const SparseMatrix<double> Ms = createFromDense(M);
  VectorXd eig_vec;
  double eig_value = 1.0f;
  smallestGenEigenSym(Ks,Ms,eig_vec,eig_value,1000,1e-8);
  ASSERT_EQ_SMALL_VEC_TOL(eig_vec, correct_eig, eig_vec.size(), 1e-10);
  ASSERT_EQ_TOL(eig_value, 2.5, 1e-7);

}

BOOST_AUTO_TEST_CASE(testGeneralEigenValueLargeProblem){
  
  SparseMatrix<double> K,M;
  TEST_ASSERT( load(K,std::string(TEST_DATA_DIR)+"beam_sparse_K.b") );
  TEST_ASSERT( load(M,std::string(TEST_DATA_DIR)+"beam_sparse_M.b") );

  VectorXd eig_vec;
  double eig_value = 1.0f;
  const int it = largestGenEigenSym(K,M,eig_vec,eig_value,1000,1e-8);
  ASSERT_GT(it,0);
  const int it2 = smallestGenEigenSym(K,M,eig_vec,eig_value,1000,1e-8);
  ASSERT_GT(it2,0);

}

BOOST_AUTO_TEST_SUITE_END()
