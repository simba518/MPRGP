#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <Simulator.h>
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(SimulatorTest)

// BOOST_AUTO_TEST_CASE(test_one_tet){
  
//   Simulator simulator;
//   simulator.init("./test_case/test_data/one_tet/collision_ball.ini");
//   simulator.print();
//   simulator.run();

//   boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();

//   double ek, ep;
//   fem_solver->getSystemEnergy(ek, ep);
//   // cout<<"ek, ep: " << setprecision(16) << ek << ", " << ep << endl;

//   ASSERT_EQ_TOL(ek, 18.0642544048533, 1e-10);
//   ASSERT_EQ_TOL(ep, -28.73275355226063, 1e-10);

//   const VVVec4d &linear_con = fem_solver->getLinearCon();

//   ASSERT_EQ(linear_con.size(), 4);
//   ASSERT_EQ(linear_con[0].size(), 0);
//   ASSERT_EQ(linear_con[1].size(), 1);
//   ASSERT_EQ(linear_con[2].size(), 1);
//   ASSERT_EQ(linear_con[3].size(), 0);

//   // cout<<"p1: " << setprecision(16) << linear_con[1][0].transpose() << endl;
//   // cout<<"p1: " << setprecision(16) << linear_con[2][0].transpose() << endl;

//   Vector4d p1,p2;
//   p1 << -0.9161056049828989, -0.2890042487666153,  0.2778975795392272,  2.270165175095369;
//   p2 << -0.9527177556000175, -0.2890044237539933,  0.0938366730819447,  2.55179605456626;
//   ASSERT_EQ_SMALL_VEC_TOL(linear_con[1][0], p1, 4, 1e-10);
//   ASSERT_EQ_SMALL_VEC_TOL(linear_con[2][0], p2, 4, 1e-10);

//   // VectorXd pos, vel, accel;
//   // fem_solver->getPos(pos);
//   // fem_solver->getPos(vel);
//   // fem_solver->getPos(accel);
  
// }

// BOOST_AUTO_TEST_CASE(test_coarse_beam_in_ball){
 
//   Simulator simulator;
//   simulator.init("./test_case/test_data/beam-coarse/collision_in_ball.ini");
//   simulator.print();
//   simulator.run();

//   boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();

//   double ek, ep;
//   fem_solver->getSystemEnergy(ek, ep);
//   // cout<<"ek, ep: " << setprecision(16) << ek << ", " << ep << endl;

//   ASSERT_EQ_TOL(ek, 3.192540512505974, 1e-10);
//   ASSERT_EQ_TOL(ep, 879.6552520662923, 1e-10);

// }

// BOOST_AUTO_TEST_CASE(test_coarse_beam_out_ball){
 
//   Simulator simulator;
//   simulator.init("./test_case/test_data/beam-coarse/collision_out_ball.ini");
//   simulator.print();
//   simulator.run();

//   boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();
//   double ek, ep;
//   fem_solver->getSystemEnergy(ek, ep);
//   // cout<<"ek, ep: " << setprecision(16) << ek << ", " << ep << endl;

//   ASSERT_EQ_TOL(ek, 44.55807589628412, 1e-10);
//   ASSERT_EQ_TOL(ep, 68.832263223794, 1e-10);
// }

BOOST_AUTO_TEST_SUITE_END()
