#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <fstream>
#include <FEMMeshFormat.h>
using namespace std;
USE_PRJ_NAMESPACE;

BOOST_AUTO_TEST_SUITE(VtkToAbq)

BOOST_AUTO_TEST_CASE(testfun){

  ifstream in("/home/simba/Workspace/CollisionHandle/data/dino/tempt_cubes_test/frame_0.vtk");
  ofstream out("/home/simba/Workspace/CollisionHandle/data/dino/tempt_cubes_test/mesh8.abq");
  FEMMeshFormat::VTKToABQ(in,out);
}

BOOST_AUTO_TEST_SUITE_END()
