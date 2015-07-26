#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MPRGPSolver.h>
#include <ICASolver.h>
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(ICAConvergencyCompare)

void compICA(const string qp){

  cout << qp << endl;

  SparseMatrix<double> A, J;
  VectorXd B, c, init_x;
  UTILITY::Timer timer;

  const bool load_succ = loadQP(A, B, J, c, init_x, qp);
  assert(load_succ);

  // print information
  if(true){
	cout << "A.rows() = " << A.rows() << endl;
	cout << "J.rows() = " << J.rows() << endl;
	cout << "A.nonZeros() = " << A.nonZeros() << endl;
	cout << "J.nonZeros() = " << J.nonZeros() << endl;
  }

  // analysis matrix A
  if(false){
	// VectorXd eig_vec;
	// double eig_val = 0.0f;
	// EIGEN3EXT::largestEigen(A,eig_vec,eig_val);
	// INFO_LOG("largest eigen value: " << eig_val);
	// EIGEN3EXT::smallestEigen(A,eig_vec,eig_val,(int)(A.rows()/10),1e-3*A.norm());
	// INFO_LOG("smallest eigen value: " << eig_val);
  }

  const double tol = 8e-5; // *B.norm()
  const int max_it = 100000;
  VectorXd uncon_x = init_x;

  // mprgp with decoupled general con
  if(true){

  	const SparseMatrix<double> JJt_mat = J*J.transpose();
  	assert_eq_ext(JJt_mat.nonZeros(), J.rows(), "Matrix J is not decoupled.\n" << J);
  	VectorXd JJt;
  	MATH::getDiagonal(JJt_mat, JJt);
  	DecoupledConProjector<double> projector(J, JJt, c);
	  
  	VectorXd x0;
  	projector.project(init_x, x0);

  	typedef FixedSparseMatrix<double> MAT;
  	MAT FA(A);

	// decoupled mprgp with preconditioning
	VectorXd x = x0;
  	timer.start();
	const int code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, x, tol, max_it, "decoupled mprgp with precond");
  	timer.stop("decoupled mprgp with con solving time: ");
	ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
  }

  // ica with con
  if(true){

  	timer.start();
  	ICASolver ica_solver(max_it, tol);
  	ica_solver.setName("ica with con");
  	ica_solver.reset(A,B);
  	VectorXd x = init_x;
  	const bool succ = ica_solver.solve(J, c, x);
  	timer.stop("ica with con solving time: ");
  	ica_solver.printSolveInfo(A, J, c, x);
  	ERROR_LOG_COND("ICA is not convergent, (iterations, residual) = " << 
  				   ica_solver.getIterations() << ", " << ica_solver.getResidual(), succ);
  }
}

BOOST_AUTO_TEST_CASE(compare){
  
  const string base_name = "/home/simba/Workspace/MPRGP/data/dino/tempt_two_mprgp_h005_";
  vector<string> qp;
  qp.push_back(base_name+"1303/QP/qpdec11.b");
  qp.push_back(base_name+"1493/QP/qpdec11.b");
  qp.push_back(base_name+"5910/QP/qpdec11.b");
  qp.push_back(base_name+"9381/QP/qpdec11.b");
  qp.push_back(base_name+"20044/QP/qpdec11.b");

  for (int i = 0; i < (int)qp.size(); ++i){
    compICA(qp[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
