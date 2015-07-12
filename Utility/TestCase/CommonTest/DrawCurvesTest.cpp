#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DrawCurves.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(DrawCurvesTest)

BOOST_AUTO_TEST_CASE(testfun){

  const VectorXd v = VectorXd::Random(10);
  TEST_ASSERT( PythonScriptDraw2DCurves<VectorXd>::write(string(TEST_DATA_DIR)+"/tempt.py",v) );

  const VectorXd v2 = VectorXd::Random(10);
  PythonScriptDraw2DCurves<VectorXd> curves;
  curves.add("v",v,1,0,"--");
  curves.add("v2",v2);
  curves.write(string(TEST_DATA_DIR)+"t.py");
}

BOOST_AUTO_TEST_SUITE_END()
