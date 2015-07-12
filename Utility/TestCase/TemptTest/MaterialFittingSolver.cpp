#include <alglib/optimization.h>
#include "nnls.h"
#include "MaterialFitting.h"
#include "QPSolver.h"

using namespace alglib;
using namespace UTILITY;
using namespace ELASTIC_OPT;

void MaterialFitting::solveByNNLS(){

  const SXMatrix M1 = assembleObjMatrix();
  vector<SX> smooth_funs;
  computeSmoothObjFunctions(smooth_funs);
  
  SXMatrix obj_fun_m(M1.size1()*M1.size2()+smooth_funs.size(),1);
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  obj_fun_m.elem(i*M1.size2()+j,0) = M1.elem(i,j);
  }
  for (int i = 0; i < smooth_funs.size(); ++i){
    obj_fun_m.elem(M1.size2()*M1.size2()+i,0) = smooth_funs[i];
  }

  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction g_fun = CasADi::SXFunction(variable_x,obj_fun_m);
  g_fun.init();
  
  VectorXd x0(variable_x.size()), b;
  x0.setZero();
  CASADI::evaluate(g_fun, x0, b);
  b = -b;

  MatrixXd J = CASADI::convert<double>(g_fun.jac());
  assert_eq(J.rows(), b.size());
  assert_eq(J.cols(), x0.size());
  
  VectorXd x(b.size()),  w(J.cols()), zz(J.rows());
  VectorXi index(J.cols()*2);
  getInitValue(x);
  
  int exit_code = 0;
  double residual = 1e-18;
  nnls(&J(0,0), J.rows(), J.rows(), J.cols(), 
  	   &b[0], &x[0], &residual,
  	   &w[0], &zz[0], &index[0], &exit_code);
  INFO_LOG("residual: " << residual);
  print_NNLS_exit_code(exit_code);

  rlst.resize(x.size());
  for (int i = 0; i < x.size(); ++i){
    rlst[i] = x[i];
  }

  const MatrixXd H = J.transpose()*J;
  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(H);
  cout << "The eigenvalues of H are: " << es.eigenvalues().transpose() << endl;
}

void MaterialFitting::solveByIpopt(){

  TRACE_FUN();
  FUNC_TIMER();

  // init solver
  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  CasADi::IpoptSolver solver = CasADi::IpoptSolver(fun);
  solver.setOption("generate_hessian",use_hessian);
  solver.setOption("tol",1e-9);
  solver.setOption("max_iter",3000);
  Timer timer;
  timer.start();
  solver.init();
  timer.stop("solver init: ");

  // set init value
  VectorXd init_x;
  getInitValue(init_x);
  solver.setInput(&init_x[0],CasADi::NLP_X_INIT);

  // set bounds
  vector<double> lower(variable_x.size()), upper(variable_x.size());
  const double lower_x = lower_bound >= 0 ? lower_bound : 0;
  const double upper_x = upper_bound > lower_x ? upper_bound : 1e20;
  for (int i = 0; i < lower.size(); ++i){
    lower[i] = lower_x;
    upper[i] = upper_x;
  }

  if (lower_bound >= 0.0f){
	solver.setInput(lower,CasADi::NLP_LBX);
  }
  if (upper_bound > 0.0f){
	assert_gt(upper_bound, lower_bound);
	solver.setInput(upper,CasADi::NLP_UBX);
  }

  // solving
  solver.solve();
  rlst.resize(variable_x.size());
  solver.getOutput(rlst,CasADi::NLP_X_OPT);

}

void MaterialFitting::solveByLinearSolver(){
  
  MatrixXd H;
  VectorXd g;
  hessGrad(H,g);

  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(H);
  cout << "H.norm() = " << H.norm() << endl;
  cout << "The eigenvalues of H are: \n" << es.eigenvalues().transpose() << endl;

  const VectorXd x = H.lu().solve(-g);
  rlst.resize(x.size());
  for (int i = 0; i < x.size(); ++i)
    rlst[i] = x[i];
  cout<< "\n(H*x+g).norm() = " << (H*x+g).norm() << "\n\n";
}

void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr){

  CasADi::SXFunction &fun = ((CasADi::SXFunction*)ptr)[0];
  CasADi::SXFunction &grad_fun = ((CasADi::SXFunction*)ptr)[1];

  VectorXd eig_x(x.length());
  for (int i = 0; i < eig_x.size(); ++i)
	eig_x[i] = x[i];
  
  static VectorXd f,g;
  CASADI::evaluate(fun, eig_x, f);
  CASADI::evaluate(grad_fun, eig_x, g);
  assert_eq(f.size(),1);
  assert_eq(g.size(), grad.length());
  assert_eq(g.norm(), g.norm());
  assert_eq(f[0], f[0]);

  func = f[0];
  for (int i = 0; i < g.size(); ++i)
    grad[i] = g[i];

  DEBUG_LOG("f = "<<func);
}

void MaterialFitting::solveByAlglib(){

  // create fun, grad
  VSX variable_x;
  CasADi::SXFunction fun_grad_fun[2];
  {
	initAllVariables(variable_x);
	fun_grad_fun[0] = CasADi::SXFunction(variable_x,objfun);
	fun_grad_fun[0].init();

	fun_grad_fun[1] = CasADi::SXFunction(variable_x,fun_grad_fun[0].jac());
	fun_grad_fun[1].init();
  }

  // init value, boundary
  real_1d_array x,bndl, bndu;
  {
	VectorXd init_x;
	getInitValue(init_x);
	x.setlength(init_x.size());
	for (int i = 0; i < init_x.size(); ++i)
	  x[i] = init_x[i];

	bndl.setlength(init_x.size());
	bndu.setlength(init_x.size());
	const double lower = lower_bound >= 0 ? lower_bound : 0;
	const double upper = upper_bound > lower ? upper_bound : 1e20;
	DEBUG_LOG("lower_bound: " << lower);
	DEBUG_LOG("upper_bound: " << upper);
	for (int i = 0; i < bndu.length(); ++i){
	  bndl[i] = lower;
	  bndu[i] = upper;
	}
  }

  minbleicstate state;
  minbleicreport rep;
  rep.terminationtype = 100;

  // stop conditions
  double epsg = 1e-12;
  double epsf = 1e-12;
  double epsx = 1e-12;
  ae_int_t maxits = 1000;

  // solve
  minbleiccreate(x, state);
  minbleicsetbc(state, bndl, bndu);
  minbleicsetcond(state, epsg, epsf, epsx, maxits);
  alglib::minbleicoptimize(state, function1_grad, NULL, fun_grad_fun);

  DEBUG_LOG("alglib exit code: "<<(int)rep.terminationtype); // EXPECTED: 4

  // get results
  {
	minbleicresults(state, x, rep);
	rlst.resize(x.length());
	for (int i = 0; i < x.length(); ++i)
	  rlst[i] = x[i];
  }

}

#include <float.h>
float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

template <typename T>
class MatrixA{
public:
  MatrixA(const SX &obj, const VSX &x){

	TRACE_FUN();
	fun = CasADi::SXFunction(x,obj);
	fun.init();

	CasADi::SXMatrix g = fun.jac();
	grad_fun = CasADi::SXFunction(x,g);
	grad_fun.init();

	VectorXd x0(x.size());
	x0.setZero();
	b.resize(x.size());
	CASADI::evaluate(grad_fun, x0, b);
	b *= -1.0f;

	// H_diag.resize(x.size());
	// assert_eq(g.size2(),x.size());
	// for (int i = 0; i < x.size(); ++i){
	//   CasADi::SXFunction gi = CasADi::SXFunction(x[i],g(0,i));
	//   gi.init();
	//   H_diag[i] = (gi.jac()).elem(0,0).getValue();
	//   assert_eq(H_diag[i],H_diag[i]);
	//   assert_gt(H_diag[i],0);
	// }
  }
  template <typename VEC,typename VEC_OUT>
  void multiply(const VEC& x,VEC_OUT& result)const{
	CASADI::evaluate(const_cast<CasADi::SXFunction &>(grad_fun), x, result);
	result += b;
  }
  double diag(const int i)const{
	assert_in(i,0,H_diag.size()-1);
	return H_diag[i];
  }
  const VectorXd &diag()const{
	H_diag;
  }
  int rows()const{
	return b.size();
  }
  const VectorXd &B()const{
	return b;
  }
  double funValue(const VectorXd &x)const{
	VectorXd rlst;
	CASADI::evaluate(const_cast<CasADi::SXFunction &>(fun), x, rlst);  
	assert_eq(rlst.size(),1);
	return rlst[0];
  }
  
private:
  CasADi::SXFunction fun;
  CasADi::SXFunction grad_fun;
  VectorXd b;
  VectorXd H_diag;
};

void MaterialFitting::solveByMPRGP(const string init_x){

  // init A
  VSX variable_x;
  initAllVariables(variable_x);
  MatrixA<double> A(objfun, variable_x);

  // set bounds
  const int n = variable_x.size();
  VectorXd lower(n),upper(n);
  const double lower_x = lower_bound >= 0 ? lower_bound : 0;
  const double upper_x = upper_bound > lower_x ? upper_bound : 1e20;
  for (int i = 0; i < lower.size(); ++i){
    lower[i] = lower_x;
    upper[i] = upper_x;
  }

  // solve
  MPRGPQPSolver<double, MatrixA<double>, Kernel<double>, MatrixA<double> > solver(A,A.B(),lower,upper,false);
  solver.setCallback(boost::shared_ptr<Callback<double, MatrixA<double> > >(new Callback<double, MatrixA<double> >(&A)));
  solver.setSolverParameters(1e-8,100);
  VectorXd x;
  getInitValue(x);
  // if ( init_x.size() <= 0 || !load(init_x,x)){
  // }

  for (int i = 0; i < x.size(); ++i){
    x[i] = x[i] > lower_x? x[i]:lower_x;
    x[i] = x[i] < upper_x? x[i]:upper_x;
  	assert_in(x[i],lower_x,upper_x);
  }

  Timer timer;
  timer.start();
  solver.solve(x);
  timer.stop("x: ");

  rlst.resize(x.size());
  for (int i = 0; i < x.size(); ++i)
    rlst[i] = x[i];
}
