#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <SparseMatrixTools.h>
#include <vector>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(SparseMatrixToolsTest)

BOOST_AUTO_TEST_CASE(testCreateFromDense){
	  
  const MatrixXd denseM = MatrixXd::Random(10,7);
  const MatrixXd sparseM = EIGEN3EXT::createFromDense(denseM);
  ASSERT_EQ_SMALL_MAT_TOL (denseM, sparseM, 1e-10);
}

BOOST_AUTO_TEST_CASE(testEye){
	  
  const MatrixXd S = eye(5,2.5f);
  const MatrixXd D = MatrixXd::Identity(5,5)*2.5f;
  ASSERT_LT( (D - S).norm(), 1e-10 );

  const MatrixXd m1 = MatrixXd::Random(2,3);
  const MatrixXd m2 = MatrixXd::Random(7,4);
  MatrixXd m(m1.rows()+m2.rows(), m1.cols()+m2.cols());
  m.setZero();
  m.topLeftCorner(m1.rows(),m1.cols()) = m1;
  m.bottomRightCorner(m2.rows(),m2.cols()) = m2;
  
  vector<MatrixXd> ms;
  ms.push_back(m1);
  ms.push_back(m2);
  SparseMatrix<double> m3;
  const MatrixXd m1m2 = eye(m3,ms);
  ASSERT_LT( (m1m2 - m).norm(), 1e-10 );
}

BOOST_AUTO_TEST_CASE(testBlock){
  
  const SparseMatrix<double> S = EIGEN3EXT::random(10,11,1.0);
  const MatrixXd A = S;

  const int r0 = 2;
  const int c0 = 3;
  const int rows = 4;
  const int cols = 5;
  const MatrixXd sub1 = EIGEN3EXT::block(S,r0,c0,rows,cols);
  ASSERT_EQ (sub1.rows(), rows);
  ASSERT_EQ (sub1.cols(), cols);

  ASSERT_LT ( (sub1-A.block(r0,c0,rows,cols)).norm(), 1e-10 );
}

BOOST_AUTO_TEST_CASE(testInverse){
	  
  const int n = 10;
  const MatrixXd A = MatrixXd::Random(n,n) + 10.0f*MatrixXd::Identity(n,n);
  const SparseMatrix<double> I = EIGEN3EXT::eye(n,(double)1.0f);
  const SparseMatrix<double> P = EIGEN3EXT::random(n,n,1.0) + I*10.0;

  SparseMatrix<double> inv_P;
  inverse(P, inv_P);
  ASSERT_EQ (P.rows(), inv_P.rows());
  ASSERT_EQ (P.cols(), inv_P.cols());
  ASSERT_LT ((I - (P*inv_P)).norm(),1e-8);
}

BOOST_AUTO_TEST_CASE(addToBlockTest){
  
  SparseMatrix<double> oldM = random(7,8,1.0);
  SparseMatrix<double> M = oldM;
  const SparseMatrix<double> sub = random(3,2,1.0);
  int r0 = 1;
  int c0 = 2;
  const int r = sub.rows();
  const int c = sub.cols();
  addToBlock(M,sub,r0,c0);
  ASSERT_EQ((M.block(r0,c0,r,c)-oldM.block(r0,c0,r,c)-sub).norm(),0.0f);
}

BOOST_AUTO_TEST_CASE(getLowerTest){
 
  const SparseMatrix<double> M = random(8,8,1.0);
  const SparseMatrix<double> L = getLower(M);
  const SparseMatrix<double> M2 = L.selfadjointView<Lower>();
  const SparseMatrix<double> M3 = M.selfadjointView<Lower>();
  ASSERT_EQ((M2-M3).norm(),0.0);
}

BOOST_AUTO_TEST_SUITE_END()
