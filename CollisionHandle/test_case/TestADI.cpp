#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ADIGroup.h>
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestADI)

BOOST_AUTO_TEST_CASE(testfun){
  
  // ADIGroup adi;
  // TEST_ASSERT( adi.loadMesh("./test_case/test_data/beam/mesh.abq") );
  // adi.split(17,9,9);
  // TEST_ASSERT( adi.save("./tempt/beam_adi_groups.txt") );
  // TEST_ASSERT( adi.checkGroups() );
  // VVVI groups;
  // TEST_ASSERT( ADIGroup::loadAdiGroups("./tempt/beam_adi_groups.txt", groups) );
  // ASSERT_EQ (groups.size(), 3);

  // ASSERT_EQ (groups[0].size(), 17);
  // ASSERT_EQ (groups[1].size(), 9);
  // ASSERT_EQ (groups[2].size(), 9);

  // ASSERT_EQ (groups[0][0].size(), 81);
  // ASSERT_EQ (groups[1][0].size(), 153);
  // ASSERT_EQ (groups[2][0].size(), 153);
  
  // ASSERT_EQ (groups[0][16].size(), 81*2);
  // ASSERT_EQ (groups[1][8].size(), 153);
  // ASSERT_EQ (groups[2][8].size(), 153);
}

BOOST_AUTO_TEST_SUITE_END()
