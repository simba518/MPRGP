#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <BBox.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(BBoxTest)

BOOST_AUTO_TEST_CASE(allTest){
  
  double v[9] = {0.0f,0.0f,0.0f,  1.0f,0.1f,2.0f,  0.5f,3.0f,0.2f  };
  Vector3d min,max;
  min << 0.0f,0.0f,0.0f;
  max << 1.0f,3.0f,2.0f;

  BBox<double,Vector3d> box(v,3);

  ASSERT_EQ(box.getMaxAxis(),3.0);
  ASSERT_EQ(box.getWidth() , 1.0);
  ASSERT_EQ(box.getDeepth() , 2.0);
  ASSERT_EQ(box.getMaxConner() , max);
  ASSERT_EQ(box.getMinConner() , min);

  double vmin[3],vmax[3];
  box.getMinConner(vmin);
  box.getMaxConner(vmax);

  ASSERT_EQ(vmin[0],min[0]);
  ASSERT_EQ(vmin[1],min[1]);
  ASSERT_EQ(vmin[2],min[2]);

  ASSERT_EQ(vmax[0],max[0]);
  ASSERT_EQ(vmax[1],max[1]);
  ASSERT_EQ(vmax[2],max[2]);

}

BOOST_AUTO_TEST_SUITE_END()
