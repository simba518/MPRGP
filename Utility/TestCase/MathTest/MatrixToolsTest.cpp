#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MatrixTools.h>
#include <iostream>
using namespace Eigen;
using namespace EIGEN3EXT;
using namespace std;

BOOST_AUTO_TEST_SUITE(MatrixTools)

BOOST_AUTO_TEST_CASE(createTest){

  const double A[] = {1,2,3,4,5,6};

  MatrixXd m;
  createFromRowMajor(m, &A[0], 2, 3);
  MatrixXd m2(2,3);
  m2 << 1,2,3,4,5,6;
  ASSERT_EQ (m,m2);

  createFromColMajor(m, &A[0], 2, 3);
  m2 << 1,3,5,2,4,6;
  ASSERT_EQ (m,m2);
}

BOOST_AUTO_TEST_CASE(SVDTest){
  
  const Matrix3d M = Matrix3d::Random();
  Matrix3d U, Vt, D;
  ModifiedSVD3x3(M,U,Vt,D);
  ASSERT_EQ_TOL (U.determinant(), 1, 1e-10);
  ASSERT_EQ_TOL (Vt.determinant(), 1, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((U.transpose()*U), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((Vt.transpose()*Vt), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((U*D*Vt), M, 1e-10);
}

BOOST_AUTO_TEST_CASE(PolarDecompositionTest){

  const Matrix3d M = Matrix3d::Random();
  Matrix3d R,S;
  ModifiedPD3x3(M,R,S);
  ASSERT_EQ_TOL (R.determinant(), 1, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((R.transpose()*R), Matrix3d::Identity(), 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL (S.transpose(), S, 1e-10);
  ASSERT_EQ_SMALL_MAT_TOL ((R*S), M, 1e-10);
}

BOOST_AUTO_TEST_CASE(convertTest){
  
  MatrixXd U(2,3);
  std::vector<VectorXd> u;
  convert(U,u);
  ASSERT_EQ((int)u.size(),U.cols());
  ASSERT_EQ(u[0].size(),U.rows());
  ASSERT_EQ(u[1].size(),U.rows());
  ASSERT_EQ(u[2].size(),U.rows());
  
  MatrixXd M;
  convert(u,M);
  ASSERT_EQ(M.rows(),U.rows());
  ASSERT_EQ(M.cols(),U.cols());

  ASSERT_EQ(M,U);
}

BOOST_AUTO_TEST_CASE(pseudoinverseTest){
  
  MatrixXd U2x2(2,2);
  U2x2 << 3,0.5,1,4;
  const MatrixXd invU2x2 = PseudoInverse(U2x2);
  const MatrixXd I2x2 = MatrixXd::Identity(2,2);
  const MatrixXd A = (U2x2*invU2x2);
  const MatrixXd B = (invU2x2*U2x2);
  ASSERT_EQ_SMALL_MAT_TOL(A,I2x2,1e-10);
  ASSERT_EQ_SMALL_MAT_TOL(B,I2x2,1e-10);

  // this example is obatined from
  // https://inst.eecs.berkeley.edu/~ee127a/book/login/exa_pinv_4by5.html
  MatrixXd U4x5(4,5);
  U4x5.setZero();
  U4x5(0,0) = 1;
  U4x5(0,4) = 2;
  U4x5(1,2) = 3;
  U4x5(3,1) = 4;
  const MatrixXd invU4x5 = PseudoInverse(U4x5);
  MatrixXd cU(5,4);
  cU.setZero();
  cU(0,0) = 0.2;
  cU(1,3) = 0.25;
  cU(2,1) = 1.0f/3.0f;
  cU(4,0) = 0.4;
  ASSERT_EQ_SMALL_MAT_TOL(invU4x5,cU,1e-8);
}

BOOST_AUTO_TEST_CASE(MGramSchmidtTest){

  const int n = 10;
  const int r = 3;
  MatrixXd M1 = MatrixXd::Random(n,n)+MatrixXd::Identity(n,n)*10.0f;
  const MatrixXd M = M1.transpose() + M1;
  MatrixXd U = MatrixXd::Random(n,r);
  
  EIGEN3EXT::MGramSchmidt(M,U);
  const MatrixXd UtMU = U.transpose()*M*U;
  const MatrixXd I = MatrixXd::Identity(r,r);
  ASSERT_EQ_SMALL_MAT_TOL(UtMU,I,1e-12);
}

BOOST_AUTO_TEST_CASE(GramSchmidtTest){

  const int n = 10;
  const int r = 3;
  MatrixXd U = MatrixXd::Random(n,r);
  
  EIGEN3EXT::GramSchmidt(U);
  const MatrixXd UtU = U.transpose()*U;
  const MatrixXd I = MatrixXd::Identity(r,r);
  ASSERT_EQ_SMALL_MAT_TOL(UtU,I,1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
