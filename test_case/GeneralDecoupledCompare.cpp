#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MPRGPSolver.h>
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(gdcompare)

void loadQPFile(const string qp, SparseMatrix<double> &A, 
				SparseMatrix<double> &J,VectorXd &B,
				VectorXd &c, VectorXd &init_x, const string ext){
  
  const bool load_succ = loadQP(A, B, J, c, init_x, qp);
  assert(load_succ);
  cout<< ext << "A.rows() = " << A.rows() << endl;
  cout<< ext << "J.rows() = " << J.rows() << endl;
  cout<< ext << "A.nonZeros() = " << A.nonZeros() << endl;
  cout<< ext << "J.nonZeros() = " << J.nonZeros() << endl;
}

BOOST_AUTO_TEST_CASE(compare){
  
  UTILITY::Timer timer;
  const string project_dir = "/home/simba/Workspace/MPRGP/";
  const string qp_base = project_dir+"/data/dino/tempt_cubes_decoupled/QP/";

  const double tol = 8e-5; // *B.norm()
  const int max_it = 1000;

  typedef FixedSparseMatrix<double> MAT;
  const int T = 140;
  SparseMatrix<double> A, J;
  VectorXd B, c, init_x, dec_ev, gen_ev;
  for (int f = 100; f < T ; ++f){

	cout << "frame: " << f << endl;
	{// decoupled mprgp
	  ostringstream qp;
	  qp << qp_base << "qpdec"<< f << ".b";
	  loadQPFile(qp.str(),A,J,B,c,init_x,"dec ");
	  MAT FA(A);
	  VectorXd x = init_x;
	  if(dec_ev.size() != init_x.size()){
		dec_ev.setRandom();
		dec_ev.normalize();
	  }
	  timer.start();
	  const int code=MPRGPDecoupledCon<double>::solve<MAT,true>(FA,B,J,c,x,tol,max_it,"decoupled",&dec_ev);
	  timer.stop("decoupled mprgp solving time: ");
	}

	{// general mprgp
	  ostringstream qp;
	  qp << qp_base << "qpfull"<< f << ".b";
	  loadQPFile(qp.str(),A,J,B,c,init_x,"gen ");
	  MAT FA(A);
	  VectorXd x = init_x;
	  if(gen_ev.size() != init_x.size()){
		gen_ev.setRandom();
		gen_ev.normalize();
	  }
	  timer.start();
	  const int code = MPRGPGeneralCon<double>::solve<MAT,false>(FA,B,J,c,x,tol,max_it,"general",&gen_ev);
	  timer.stop("general mprgp solving time: ");
	}
  }
}

BOOST_AUTO_TEST_SUITE_END()
