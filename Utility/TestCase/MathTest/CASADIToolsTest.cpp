#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CASADITools.h>
using namespace Eigen;
using namespace CASADI;
using casadi::SXElement;
using casadi::SX;

BOOST_AUTO_TEST_SUITE(CASADITools)

BOOST_AUTO_TEST_CASE(convert2VecTest){

  VectorXd v(3);
  v << 1,2,5;
  const casadi::SX SV = convert(v);
  ASSERT_EQ(SV.size1(),v.size());
  ASSERT_EQ(SV.size2(),1);
  const VectorXd a = convert2Vec<double>(SV);
  ASSERT_EQ(a,v);
}

BOOST_AUTO_TEST_CASE(makeEyeMatrixTest){
  
  VectorXd a(3);
  a << 1,2,3;
  const casadi::SX A = makeEyeMatrix(a);
  const MatrixXd Ms = convert<double>(A);
  MatrixXd M(3,3);
  M.setZero();
  M(0,0)=1;   M(1,1)=2;   M(2,2)=3;
  ASSERT_EQ(Ms,M);
}

BOOST_AUTO_TEST_CASE(convertTest){

  VectorXd v(6);
  v << 1,2,3,4,5,6;
  VSX vec;
  for (int i = 0; i < v.size(); ++i){
    vec.push_back(v[i]);
  }

  const SX M1 = convert(vec);
  ASSERT_EQ(M1.size1(),6);
  ASSERT_EQ(M1.size2(),1);
  const VectorXd m1 = convert2Vec<double>(M1);
  ASSERT_EQ(m1,v);

  const SX M2 = convert(vec,2);
  ASSERT_EQ(M2.size1(),3);
  ASSERT_EQ(M2.size2(),2);
  const MatrixXd m2 = convert<double>(M2);
  const Matrix<double,3,2> &em2 = Map<Matrix<double,3,2> >(&(v[0]));
  ASSERT_EQ(m2,em2);

  const SX M3 = convert(vec,2,false);
  ASSERT_EQ(M3.size1(),3);
  ASSERT_EQ(M3.size2(),2);
  const MatrixXd m3 = convert<double>(M3);
  Matrix<double,3,2> em3(3,2);
  em3 << 1,2,3,4,5,6;
  ASSERT_EQ(m3,em3);
  
}

BOOST_AUTO_TEST_SUITE_END()
