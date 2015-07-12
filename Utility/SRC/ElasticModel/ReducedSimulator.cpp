#include "ReducedSimulator.h"
#include <JsonFilePaser.h>
#include <Log.h>
#include <ConMatrixTools.h>
#include <Timer.h>
using namespace UTILITY;
using namespace SIMULATOR;

bool ReducedSimulator::init(const string filename){

  JsonFilePaser jsonf;
  bool succ = false;
  if (jsonf.open(filename)){
	succ = true;
	jsonf.read("alpha_k",alpha_k,0.0);
	jsonf.read("alpha_m",alpha_m,0.0);
	jsonf.read("h",h,1.0);
  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

void ReducedImpLogConSimulator::setConGroups(const vector<int> &con_nodes){

  if (con_nodes.size() <= 0){
	removeAllCon();
	return ;
  }

  assert_ge(fullDim(), con_nodes.size()*3);
  SparseMatrix<double> sC;
  computeConM(con_nodes, sC, fullDim()/3);
  assert(model);
  assert_eq(sC.cols(), model->getModalBasis().rows());
  C = sC*model->getModalBasis();

  const int r = reducedDim();
  const int rc = C.rows();
  const MatrixXd hC = h*C;

  A.resize(r+rc,r+rc);
  A.setZero();
  A.topRightCorner(C.cols(),C.rows()) = hC.transpose();
  A.bottomLeftCorner(C.rows(),C.cols()) = hC;
  
  b.resize(r+rc);
}

void ReducedImpLogConSimulator::removeAllCon(){
  
  uc.resize(0);
  C.resize(0,0);
  const int r = reducedDim();
  A.resize(r,r);
  b.resize(r);
}

bool ReducedImpLogConSimulator::forward(){

  assert(model);
  assert_eq(C.rows(), uc.size());
  const int r = reducedDim();
  assert_eq(q.size(), r);
  assert_eq(v.size(), r);
  assert_eq(red_fext.size(), r);
  assert_eq(model->getReducedMassM().rows(), r);

  static MatrixXd L;
  static VectorXd s;
  model->evaluateFK(q,s,L);

  L *= (h*h+h*alpha_k);
  L += (1.0f+h*alpha_m)*model->getReducedMassM();

  s *= (-1.0f*h);
  s += h*red_fext;
  s += model->getReducedMassM()*v;

  if (C.rows() <= 0){
	A = L;
	b = s;
  }else{
	assert_eq(b.size(), r+C.rows());
	assert_eq(A.rows(), b.size());
	assert_eq(A.cols(), b.size());
	A.topLeftCorner(r,r) = L;
	b.head(r) = s;
  }

  assert_eq(b.size(), r+C.rows());
  if (C.rows() > 0){
	b.tail(C.rows()) = uc-C*q;
  }

  v = A.lu().solve(b).head(r);
  q += h*v;

  return true;
}

bool ReducedStaticPenConSimulator::init(const string init_filename){

  bool succ = ReducedSimulator::init(init_filename);
  if (succ){
	UTILITY::JsonFilePaser jsonf;
	succ = jsonf.open(init_filename);
	if (!succ){
	  ERROR_LOG("failed to open" << init_filename);
	}else{
	  jsonf.read("con_penalty", lambda, 100.0);
	  jsonf.read("max_iter", max_it, 5);
	  jsonf.read("tolerance", tolerance, 0.1);
	}
  }
  return succ;
}

void ReducedStaticPenConSimulator::setConGroups(const vector<int> &con_nodes){

  if (con_nodes.size() <= 0){
	removeAllCon();
	return ;
  }

  assert_ge(fullDim(), con_nodes.size()*3);
  SparseMatrix<double> sC;
  computeConM(con_nodes, sC, fullDim()/3);
  assert(model);
  assert_eq(sC.cols(), model->getModalBasis().rows());
  C = sC*model->getModalBasis();
  lambda_CtC = lambda*(C.transpose()*C);
}

void ReducedStaticPenConSimulator::setUc(const VectorXd &uc){

  lambda_CtUc = lambda*(C.transpose()*uc);
}

void ReducedStaticPenConSimulator::removeAllCon(){

  C.resize(0,0);
}

bool ReducedStaticPenConSimulator::forward(){

  // solve the nonlinear equation using newton method.
  bool succ = true;
  for (int i = 0; i < max_it; ++i){

	const VectorXd &f = grad(q);
	const MatrixXd &G = jac(q);
	const VectorXd p = G.lu().solve(f);
	q -= p;
	if(p.norm() <= tolerance){
	  break;
	}
  }
  return succ;
}

const VectorXd &ReducedStaticPenConSimulator::grad(const VectorXd &q){

  model->evaluateF(q, g);
  g -= red_fext;
  if (C.size() > 0){
	assert_eq(lambda_CtUc.size(), g.size());
	assert_eq(lambda_CtC.rows(), g.size());
	g += lambda_CtC*q - lambda_CtUc;
  }
  return g;
}

const MatrixXd &ReducedStaticPenConSimulator::jac(const VectorXd &q){

  model->evaluateK(q, J);
  if (C.size() > 0){
	assert_eq(J.rows(), lambda_CtC.rows());
	assert_eq(J.cols(), lambda_CtC.cols());
	J += lambda_CtC;
  }
  return J;
}
