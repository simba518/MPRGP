#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <TextForRender.h>
using namespace Eigen;
using namespace QGLVEXT;

BOOST_AUTO_TEST_SUITE(TextForRenderTest)

BOOST_AUTO_TEST_CASE(testAll){

  TextForRender text;
  text.insert(TextWithPosition("h",1,2));
  text.update(TextWithPosition("c",1,2));
  ASSERT_EQ (text.size(),1);
  ASSERT_EQ (text.begin()->getContent(),"c");
}

BOOST_AUTO_TEST_SUITE_END()
