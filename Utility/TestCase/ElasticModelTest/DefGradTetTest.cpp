#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DefGradTet.h>
#include "../CommonTest/TetMeshFactoryForTest.h"
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(DefGradTetTest)

BOOST_AUTO_TEST_CASE(testfun){
  
  DefGradTet def_grad;
  pTetMesh_const vol_mesh = TetMeshFactoryForTest::tet2();
  def_grad.prepare(vol_mesh);
  ASSERT_EQ( def_grad.invDm().size(), 2);
  ASSERT_EQ( def_grad.dF().size(), 2);
  
  VectorM3x3 F;
  const VectorXd X = VectorXd::Random(vol_mesh->nodes().size()*3);
  def_grad.evalF(F, X);
  ASSERT_EQ(F.size(),2);
  ASSERT_GT((F[0]-Matrix3d::Identity()).norm(),0);

  VectorXd X0;
  vol_mesh->nodes(X0);
  def_grad.evalF(F, X0);
  ASSERT_EQ(F.size(),2);
  ASSERT_EQ(F[0],Matrix3d::Identity());
  ASSERT_EQ(F[1],Matrix3d::Identity());
}

BOOST_AUTO_TEST_SUITE_END()
