#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MPRGPSolver.h>
#include "TestUtility.h"
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(MPRGPTest)

BOOST_AUTO_TEST_CASE(test_projector){

  const int m = 3;
  const int n = 4;

  const MatrixXd JM = MatrixXd::Identity(m,n);
  const SparseMatrix<double> J = createFromDense(JM);
  const SparseMatrix<double> JJt_mat=J*J.transpose();
  VectorXd JJt;
  getDiagonal(JJt_mat, JJt);

  const VectorXd c = VectorXd::Random(m);

  DecoupledConProjector<double> dp(J, JJt, c);
  GeneralConProjector<double> gp(J, c);
  const VectorXd x = VectorXd::Random(n);
  const VectorXd g = VectorXd::Random(n);

  {// project, decide face
	VectorXd y_d, y_g;
	gp.project(x, y_g);
	dp.project(x, y_d);
	ASSERT_EQ(y_d.transpose(), y_g.transpose());

	gp.DECIDE_FACE(y_g);
	dp.DECIDE_FACE(y_d);
	// cout << "x: " << y_g.transpose() << endl;
	// cout << "Jx-c: "<< (J*y_g-c).transpose() << endl;
	// cout << "f: ";
	for (int i = 0; i < m; ++i){
	  // cout << (int)gp.getFace()[i] << " ";
	  ASSERT_EQ( gp.getFace()[i], dp.getFace()[i]);
	}
	// cout << endl;
  }

  {// phi
	VectorXd phi_g, phi_d;
	gp.PHI(g, phi_g);
	dp.PHI(g, phi_d);
	// cout<< "J: " << gp.getConMatrix() << endl;
	// cout<< "Jh: " << gp.getActiveConMatrix() << endl;
	// cout<< "J*phi: " << (gp.getActiveConMatrix()*phi_d).transpose() << endl;
	ASSERT_EQ(phi_d.transpose(), phi_g.transpose());
  }

  {// beta
	const VectorXd phi = VectorXd::Random(n);
	VectorXd beta_p, beta_d;
	gp.BETA(g, beta_p, phi);
	dp.BETA(g, beta_d, phi);
	cout << "g-bd: " << (g-beta_d).transpose() << endl;
	cout << "g-bp: " << (g-beta_p).transpose() << endl;
	ASSERT_EQ(beta_d.transpose(), beta_p.transpose());
  }
}

BOOST_AUTO_TEST_CASE(test_solver){

  const int n = 5;
  const MatrixXd M = MatrixXd::Random(n,n) + MatrixXd::Identity(n,n)*3.0f;
  const MatrixXd MtM = M.transpose()*M;
  const SparseMatrix<double> A = createFromDense(MtM);

  const MatrixXd JM = MatrixXd::Identity(n,n);
  const SparseMatrix<double> J = createFromDense(JM);

  const VectorXd B = VectorXd::Random(n);
  const VectorXd c = VectorXd::Random(n);
  VectorXd x = VectorXd::Random(n);

  typedef FixedSparseMatrix<double> MAT;
  MAT FA(A);
  {
	cout << "general: " << endl;
	VectorXd y = x;
	const int rlst_code = MPRGPGeneralCon<double>::solve<MAT,true>(FA,B,J,c,y);
	ASSERT_EQ(rlst_code, 0);
  }
  {
	cout << "decoupled: " << endl;
	VectorXd y = x;
	const int rlst_code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA,B,J,c,y);
	ASSERT_EQ(rlst_code, 0);
  }
 }

BOOST_AUTO_TEST_CASE(test_solver_large){
 
  SparseMatrix<double> A, J_all;
  VectorXd B, c_all, x;
  const string qp_file = "./test_case/test_data/qp33.b";
  const bool load_qp_succ = loadQP(A, B, J_all, c_all, x, qp_file);
  const int m = std::min(20,J_all.rows());
  const SparseMatrix<double> J = J_all.topRows(m);
  const VectorXd c = c_all.head(m);

  TEST_ASSERT(load_qp_succ);
  typedef FixedSparseMatrix<double> MAT;
  MAT FA(A);

  cout << "dimension: " << A.rows() << endl;
  cout << "constraints: " << J.rows() << endl;

  // solve
  VectorXd y2 = x;
  const int rlst_code2 = MPRGPGeneralCon<double>::solve<MAT,true>(FA, B, J, c, y2);

  VectorXd y1 = x;
  const int rlst_code1 = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, y1);

  cout << "diff: " << (y1-y2).norm() << endl;

  // check results
  ASSERT_EQ(rlst_code1, 0);
  ASSERT_EQ(rlst_code2, 0);
  const VectorXd cy1 = J*y1-c;
  const VectorXd cy2 = J*y2-c;
  for (int i = 0; i < cy1.size(); ++i){
	ASSERT_GE(cy1[i],0);
	ASSERT_GE(cy2[i],0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
