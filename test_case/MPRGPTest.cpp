#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MPRGPSolver.h>
#include "TestUtility.h"
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(MPRGPTest)

BOOST_AUTO_TEST_CASE(test_solver){

  const int n = 5;
  const MatrixXd M = MatrixXd::Random(n,n) + MatrixXd::Identity(n,n)*3.0f;
  const MatrixXd MtM = M.transpose()*M;
  const SparseMatrix<double> A = createFromDense(MtM);
  const MatrixXd JM = MatrixXd::Identity(std::min<int>(A.rows(), A.rows()), A.rows());
  const SparseMatrix<double> J = createFromDense(JM);

  const VectorXd B = VectorXd::Random(n);
  const VectorXd c = VectorXd::Random(n);
  VectorXd x = VectorXd::Random(n);

  typedef FixedSparseMatrix<double> MAT;
  MAT FA(A);
  const int rlst_code = MPRGPGeneralCon<double>::solve<MAT,false>(FA,B,J,c,x);
  ASSERT_EQ(rlst_code, 0);
 }

BOOST_AUTO_TEST_SUITE_END()
