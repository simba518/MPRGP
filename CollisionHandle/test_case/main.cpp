#define BOOST_TEST_MODULE "C++ Unit Tests for CollisionHandle Lib"
#include <boost/test/included/unit_test.hpp>

#include <float.h>
#include <limits>
#include <MPRGPSolver.h>
using namespace MATH;

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;
