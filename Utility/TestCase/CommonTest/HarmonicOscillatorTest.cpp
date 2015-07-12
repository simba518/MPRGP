#include <vector>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <HarmonicOscillator.h>
using namespace UTILITY;
using namespace std;

BOOST_AUTO_TEST_SUITE(HarmonicOscillatorTest)

BOOST_AUTO_TEST_CASE(generateSequenceTest){

  const double lambda = 0.5;
  const double alpha_k = 0.1;
  const double alpha_m = 0.2;
  const double z0 = 0.5;
  const double v0 = 0.8;

  HarmonicOscillator<double> z(lambda,alpha_k,alpha_m,z0,v0);
  ASSERT_EQ(z(0),0.5);
  ASSERT_EQ_TOL(z(1.83),1.059137778977706,1e-10);

  vector<double> Lambda(2,lambda), Z0(2,z0), V0(2,v0);
  HarmonicOscillatorSet<double> Z(Lambda,alpha_k,alpha_m,Z0,V0);

  const vector<double> q = Z.at<vector<double> >(0);
  ASSERT_EQ(q.size(),2);
  ASSERT_EQ(q[0],0.5);
  ASSERT_EQ(q[1],0.5);

  const vector<double> qq = Z.at<vector<double> >(1.83);
  ASSERT_EQ(qq.size(),2);
  ASSERT_EQ_TOL(qq[0],1.059137778977706,1e-10);
  ASSERT_EQ_TOL(qq[1],1.059137778977706,1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
