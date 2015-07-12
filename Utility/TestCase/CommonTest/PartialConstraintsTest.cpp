#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <PartialConstraints.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(PartialConstraintsTest)

BOOST_AUTO_TEST_CASE(testIO){
  
  PartialConstraintsSet parcons;
  TEST_ASSERT(parcons.load(string(TEST_DATA_DIR)+"PartialConTest.txt"));
  ASSERT_EQ( parcons.getPartialConSet().size(),2);
  for (int f = 0; f < 2; ++f){
	pPartialConstraints_const p = parcons.getPartialCon(f);
	TEST_ASSERT( p!= NULL);
	if(p){
	  ASSERT_EQ( p->numConNodes(),61 );
	  ASSERT_EQ( *(p->getConNodesSet()[0].begin()), 170 );
	}
  }
  TEST_ASSERT(parcons.write(string(TEST_DATA_DIR)+"tempt_PartialConTest.txt"));
}

BOOST_AUTO_TEST_SUITE_END()
