#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MPRGPSolver.h>
#include <ICASolver.h>
#include <FlexiblePCG.h>
#include <EigenSolver.h>
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(ConvergencyCompare)

BOOST_AUTO_TEST_CASE(test_iterative_solvers){
  
  const string project_dir = "/home/simba/Workspace/CollisionHandle/";
  const string qp = project_dir+"/data/bunny/tempt_one2/QP/qp11.b";
  // const string qp = project_dir+"/data/dino/tempt_cubes/QP/qp201.b";
  // const string qp = project_dir+"/data/dragon/tempt_selfcon_mprgp/QP/qp60.b";
  // const string qp = project_dir+"/data/longcube/tempt_mprgp/QP/qp20.b";
  
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
	return;
  }

  // analysis matrix A
  if(false){
	VectorXd eig_vec;
	double eig_val = 0.0f;
	EIGEN3EXT::largestEigen(A,eig_vec,eig_val);
	INFO_LOG("largest eigen value: " << eig_val);
	EIGEN3EXT::smallestEigen(A,eig_vec,eig_val,(int)(A.rows()/10),1e-3*A.norm());
	INFO_LOG("smallest eigen value: " << eig_val);
  }

  const double tol = 1e-3; // *B.norm()
  const int max_it = 10000;
  VectorXd uncon_x = init_x;

  // un-constrained solving with eigen cg
  if (false){
	ConjugateGradient<SparseMatrix<double> > cg_solver;
	timer.start();
	uncon_x = cg_solver.compute(A).solve(B);
	timer.stop("eigen-cg solving time: ");
  }

  // mprgp with decoupled general con
  if(false){

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
	if(true){
	  VectorXd x = x0;
	  const int code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, x, tol, max_it, "decoupled mprgp with precond");
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
	}

	// decoupled mprgp without preconditioning
	if(false){
	  VectorXd x = x0;
	  const int code = MPRGPDecoupledCon<double>::solve<MAT,false>(FA, B, J, c, x, tol, max_it, "decoupled mprgp no precond");
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
	}

	// without constraints
	if(true){
	  SparseMatrix<double> J(0, init_x.size());
	  VectorXd c;
	  VectorXd x = x0;
	  const int code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, x, tol, max_it, "decoupled mprgp without con");
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
	}
  }

  // ica with con
  if(false){

  	timer.start();
  	ICASolver ica_solver(max_it, tol);
  	ica_solver.setName("ica with con");
  	ica_solver.reset(A,B);
  	// VectorXd x = init_x;
  	VectorXd x = uncon_x;
  	const bool succ = ica_solver.solve(J, c, x);
  	timer.stop("ica with con solving time: ");
  	ica_solver.printSolveInfo(A, J, c, x);
  	ERROR_LOG_COND("ICA is not convergent, (iterations, residual) = " << 
  				   ica_solver.getIterations() << ", " << ica_solver.getResidual(), succ);
  }

  // gs without con
  if(false){
  	GaussSeidel<false> GS(max_it, tol);
  	timer.start();
  	GS.setName("gs without con");
  	GS.reset(A);
  	VectorXd x = init_x;
  	const bool succ = GS.solve(B, x);
  	timer.stop("GS solving time: ");
  	GS.printSolveInfo(A, B, x);
  	ERROR_LOG_COND("GS is not convergent\n", succ);
  }

  // block-gs without con
  if(false){
  	BlockGaussSeidel BlockGS(max_it, tol);
  	BlockGS.setName("block without con");
  	BlockGS.reset(A);
  	VectorXd x = init_x;
  	const bool succ = BlockGS.solve(B, x);
  	// BlockGS.printSolveInfo(A, B, x);
  	ERROR_LOG_COND("Block GS is not convergent\n", succ);
  }

}

BOOST_AUTO_TEST_SUITE_END()
