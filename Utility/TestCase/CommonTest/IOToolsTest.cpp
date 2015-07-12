#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <IOTools.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(IOTools)

BOOST_AUTO_TEST_CASE(vectorIO)
{
  std::vector<double> v1,v2;
  v1.push_back(1);
  v1.push_back(2);
  v1.push_back(3);
  const std::string fname = std::string(TEST_DATA_DIR)+"/temptV.b";
  TEST_ASSERT(writeVec(fname,v1));
  TEST_ASSERT(loadVec(fname,v2));
  ASSERT_EQ_SMALL_VEC(v1,v2,(int)v1.size());
}

BOOST_AUTO_TEST_CASE(matrixIO)
{
  MatrixXd m1 = MatrixXd::Random(3,2)*2;
  MatrixXd m2 = MatrixXd::Random(3,2);
  ASSERT_NE(m1,m2);
  const std::string fname = std::string(TEST_DATA_DIR)+"temptM.b";
  TEST_ASSERT(writeMat(fname,m1));
  TEST_ASSERT(loadMat(fname,m2));
  ASSERT_EQ_SMALL_MAT(m1,m2);
}

BOOST_AUTO_TEST_SUITE_END()
