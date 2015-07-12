#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ExtendModalBasis.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace SIMULATOR;

BOOST_AUTO_TEST_SUITE(ExtendModalBasisTest)

BOOST_AUTO_TEST_CASE(testExtendNonRigidBasis){
  
  MatrixXd linearBasis, extendedBasis;
  linearBasis.resize(6,2);
  linearBasis.col(0) << 1,2,3,4,5,6;
  linearBasis.col(1) = linearBasis.col(0)*10;

  ExtendModalBasis::extendNonRigidBasis(linearBasis, extendedBasis);
  // cout << linearBasis << endl << endl;
  // cout << extendedBasis << endl;
}

BOOST_AUTO_TEST_CASE(testExtendRigidBasis){
  
  MatrixXd extendedBasis;
  const VectorXd rest_shape = VectorXd::Random(4*3);
  ExtendModalBasis::extendRigidBasis(rest_shape, extendedBasis);
  // cout << extendedBasis << endl;
  EIGEN3EXT::GramSchmidt(extendedBasis);
  // cout << (extendedBasis.transpose()*extendedBasis) << endl;
  // cout << extendedBasis << endl;
}

BOOST_AUTO_TEST_CASE(testConstruct){
  
  MatrixXd linearBasis, extendedBasis;
  linearBasis = MatrixXd::Random(6*3,2);
  ExtendModalBasis::construct(linearBasis, extendedBasis);
}

BOOST_AUTO_TEST_SUITE_END()
