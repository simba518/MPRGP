#include <AuxTools.h>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(AuxToolsTest)

BOOST_AUTO_TEST_CASE(testConvert){
  
  const string a = TOSTR(20);
  ASSERT_NE(a,string("30"));
  ASSERT_EQ(a,string("20"));

  const int b = TOINT(a);
  ASSERT_EQ(b,20);

  const double c = TOINT(a);
  ASSERT_EQ(c,20.0f);
}

BOOST_AUTO_TEST_CASE(testMath){
  
  std::vector<double> a;
  ASSERT_EQ(UTILITY::norm2<vector<double> >(a),0);
  a.push_back(1);
  ASSERT_EQ(UTILITY::norm2<vector<double> >(a),1);
  a.push_back(-2);
  ASSERT_EQ_TOL(UTILITY::norm2<vector<double> >(a),sqrt(1+2*2),1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
