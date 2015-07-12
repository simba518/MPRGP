#ifndef _MATERIALFITTING_H_
#define _MATERIALFITTING_H_

#include <iostream>
#include <MatrixTools.h>
#include <MatrixIO.h>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
#include <MassMatrix.h>
#include <ElasticForceTetFullStVK.h>
#include <Timer.h>
#include "ComputeStiffnessMat.h"
#include "ComputeMassMat.h"
using namespace std;
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace ELASTIC_OPT;
using namespace CASADI;
using CasADi::SX;

namespace ELASTIC_OPT{
  
  // fit G,L,and density by using MA,
  // K(G,L)*W = M(\rho)*W*\Lambda
  class MaterialFitting{
	
  public:
	MaterialFitting(){
	  tetmesh = pTetMesh(new TetMesh());
	  use_hessian = false;
	  mu_neigh_G = 0.0f;
	  mu_neigh_L = 0.0f;
	  mu_neigh_rho = 0.0f;
	  mu_average_rho = 0.0f;
	  mu_neigh_poission = 0.0f;
	  mu_neigh_E = 0.0f;
	  lower_bound = 0.0f;
	  upper_bound = -1.0f;
	  scalor = 1.0f;
	  scaled_lambda = 1.0f;
	  scaled_mass = 1.0f;
	  fullspace_con_pen = 0.0f;
	}
	void loadTetMesh(const string filename){
	  const bool succ = tetmesh->load(filename); assert(succ);
	}
	void loadMtl(const string filename){
	  bool succ = tetmesh->loadElasticMtl(filename); assert(succ);
	}
	void loadFixednodes(const string filename){
	  vector<int> fixednodes_vec;
	  const bool succ = loadVec(filename,fixednodes_vec, TEXT); assert(succ);
	  fixednodes.clear();
	  for (int i = 0; i < fixednodes_vec.size(); ++i){
		fixednodes.insert(fixednodes_vec[i]);
	  }
	}
	virtual void setWLambda(const MatrixXd &eigen_W,const VectorXd &eigen_lambda){

	  CASADI::convert(eigen_W, W);
	  lambda.resize(eigen_lambda.size(), eigen_lambda.size());
	  ScaleMat = MatrixXd::Ones(eigen_lambda.size(), eigen_lambda.size());
	  for (int i = 0; i < eigen_lambda.size(); ++i)
		lambda.elem(i,i) = eigen_lambda[i];
	  assert_lt(W.size2(), W.size1());
	  assert_eq(lambda.size1(), W.size2());
	  assert_eq(lambda.size2(), W.size2());
	}
	void setScaleMatrix(const MatrixXd &Scale){
	  ScaleMat = Scale;
	}
	void useHessian(const bool use){
	  use_hessian = use;
	}
	void setBounds(const double lower, const double upper){
	  lower_bound = lower;
	  upper_bound = upper;
	}
	void setFullSpaceConPen(const double pen){
	  assert_ge(pen,0);
	  fullspace_con_pen = pen;
	}
	void setMuSmoothGL(const double mu_G,const double mu_L){
	  assert_ge(mu_G,0.0f);
	  assert_ge(mu_L,0.0f);
	  mu_neigh_G = mu_G;
	  mu_neigh_L = mu_L;
	}	
	void setMuSmoothEv(const double mu_E,const double mu_v){
	  assert_ge(mu_E,0.0f);
	  assert_ge(mu_v,0.0f);
	  mu_neigh_E = mu_E;
	  mu_neigh_poission = mu_v;
	}
	void setMuSmoothDensity(const double mu_rho){
	  assert_ge(mu_rho,0.0f);
	  mu_neigh_rho = mu_rho;
	}
	void setMuAverageDensity(const double mu){
	  assert_ge(mu,0.0f);
	  mu_average_rho = mu;
	}
	void computeK();
	void computeM();
	void computeK(SparseMatrix<double> &K)const;
	void computeK(SparseMatrix<double> &K,const vector<double>&G,const vector<double>&L)const;
	void computeM(DiagonalMatrix<double,-1> &M)const;
	void computeM(DiagonalMatrix<double,-1> &M,const vector<double> &rho)const;
	void computeM(SparseMatrix<double> &M)const;
	void computeM(SparseMatrix<double> &M,const vector<double> &rho)const;
	void removeFixedDOFs();
	SparseMatrix<double> getMatrixForRemovingFixedDOFs()const;

	void hessGrad(MatrixXd &H, VectorXd &g)const;
	void hess(MatrixXd &H)const;
	void grad(VectorXd &g)const;

	void scale();
	void unScale();
	virtual void assembleObjfun();
	void solveByIpopt();
	void solveByNNLS();
	void solveByLinearSolver();
	void solveByAlglib();
	void solveByMPRGP(const string init_x="");

	virtual void saveResults(const string filename)const;
	virtual vector<double> getShearGResult()const;
	virtual vector<double> getLameResult()const;
	virtual vector<double> getDensityResult()const;
	void testFixedNodes();
	void printResult()const;
	int fixedNodeBefore(const int i)const{
	  int b = 0;
	  BOOST_FOREACH(int ele, fixednodes){
		if(i < ele)
		  break;
		b ++;
	  }
	  return b;
	}
	bool isFixed(const int i)const{
	  return fixednodes.find(i) != fixednodes.end();
	}

	bool saveAllInputs(const string mesh_name)const;
	
  protected:
	virtual void initShearG(const int num_tet, VSX &G){
	  G = CASADI::makeSymbolic(num_tet, "G");
	}
	virtual void initLame(const int num_tet, VSX &Lame){
	  Lame = CASADI::makeSymbolic(num_tet, "lame");
	}
	virtual void initDensity(const int num_tet, VSX &rho){
	  TRACE_FUN();
	  rho = CASADI::makeSymbolic(num_tet, "rho");
	}
	virtual void initAllVariables(VSX &x)const{
	  x = CASADI::connect(CASADI::connect(G,Lame),rho);
	}
	virtual void getInitValue(VectorXd &init_x)const;
	virtual SXMatrix assembleObjMatrix();
	void addFullSpaceConstraints(SX &objfun)const;

	virtual VSX getG()const{return G;}
	virtual VSX getLame()const{return Lame;}
	virtual VSX getRho()const{return rho;}
	virtual VSX getE()const{return getYoungE(G,Lame);}
	virtual VSX getV()const{return getPoissonV(G,Lame);}

	void addSmoothObjfun(SX &objfun)const;
	void computeSmoothObjFunctions(vector<SX> &funs)const;
	void addAverageDensityObjfun(SX &objfun)const;
	void addFixedNodesObjfun(SX &objfun)const;
	void print_NNLS_exit_code(const int exit_code)const;

	template<class T>
	T getPoissonV(const T&G, const T&L)const{return L/(2.0*(L+G));}
	template<class T>
	T getYoungE(const T&G, const T&L)const{return G*(3.0*L+2.0*G)/(L+G);}
	template<class T>
	vector<T> getPoissonV(const vector<T>&G, const vector<T>&L)const{
	  assert_eq(G.size(),L.size());
	  vector<T> v(G.size());
	  for (size_t i = 0; i < G.size(); ++i)
		v[i] = getPoissonV(G[i],L[i]);
	  return v;
	}
	template<class T>
	vector<T> getYoungE(const vector<T>&G, const vector<T>&L)const{
	  assert_eq(G.size(),L.size());
	  vector<T> E(L.size());
	  for (size_t i = 0; i < L.size(); ++i)
		E[i] = getYoungE(G[i],L[i]);
	  return E;
	}
	
  protected:
	SXMatrix W;
	SXMatrix lambda;
	MatrixXd ScaleMat;
	pTetMesh tetmesh;
	set<int> fixednodes;
	bool use_hessian;
	double mu_neigh_G, mu_neigh_L, mu_neigh_rho,mu_average_rho,mu_neigh_poission,mu_neigh_E;
	SX objfun;
	CASADI::VSX G;
	CASADI::VSX Lame;
	CASADI::VSX rho;
	SXMatrix K;
	SXMatrix M;
	SX scalor;
	vector<double> rlst;
	double lower_bound, upper_bound, scaled_mass, scaled_lambda, fullspace_con_pen;
  };
  typedef boost::shared_ptr<MaterialFitting> pMaterialFitting;

  // fit only G, L using MA.
  // K(G,L)*W = M*W*\Lambda
  class MaterialFitting_MA_K: public MaterialFitting{

  public:
	MaterialFitting_MA_K():MaterialFitting(){
	  mu_mass = 1.0f;
	  mu_stiff = 1.0f;
	}
	vector<double> getDensityResult()const{
	  return tetmesh->material()._rho;
	}

  protected:
	void initDensity(const int num_tet, VSX &rho);
	void initAllVariables(VSX &x)const{
	  x = CASADI::connect(G,Lame);
	}
	void getInitValue(VectorXd &init_x)const;

  private:
	double mu_mass;
	double mu_stiff;
  };

  // fit G,L,and density by diagonalizing K, M.
  // W^t*K(G,L)*W = \Lambda
  // W^t*M(\rho)*W = I
  class MaterialFitting_Diag_KM: public MaterialFitting{

  public:
	MaterialFitting_Diag_KM():MaterialFitting(){
	  mu_mass = 1.0f;
	  mu_stiff = 1.0f;
	}
	void setMuMass(const double mu){
	  assert_ge(mu,0.0f);
	  mu_mass = mu;
	}
	void setMuStiff(const double mu){
	  assert_ge(mu,0.0f);
	  mu_mass = mu;
	}
	void assembleObjfun();

  private:
	double mu_mass;
	double mu_stiff;
  };

  // fit only density by diagonalizing M.
  // W^t*M(\rho)*W = I
  class MaterialFitting_Diag_M: public MaterialFitting{

  public:
	vector<double> getShearGResult()const{return tetmesh->material()._G;}
	vector<double> getLameResult()const{return tetmesh->material()._lambda;}
	vector<double> getDensityResult()const;

  protected:
	virtual SXMatrix assembleObjMatrix();
	void initShearG(const int num_tet, VSX &G);
	void initLame(const int num_tet, VSX &Lame);
	virtual void initAllVariables(VSX &x)const{x = rho;}
	virtual void getInitValue(VectorXd &init_x)const;
  };

  // fit only G, L by diagonalizing K.
  // W^t*K(G,L)*W = \Lambda
  class MaterialFitting_Diag_K: public MaterialFitting_MA_K{
	
  protected:
	SXMatrix assembleObjMatrix();
  };

  // use Young's parameter E and Poisson's ratio v as parameters.
  // W^t*K(E,v)*W = \Lambda
  class MaterialFitting_EV:public MaterialFitting{
	
  public:
	vector<double> getShearGResult()const{return getShareG(getYoungE(),getPoissonV());}
	vector<double> getLameResult()const{return getLameL(getYoungE(),getPoissonV());}
	
  protected:
	void initShearG(const int num_tet, VSX &G){
	  init_Ev(num_tet,E,v);
	  G = getShareG(E,v);
	}
	void initLame(const int num_tet, VSX &Lame){
	  init_Ev(num_tet,E,v);
	  Lame = getLameL(E,v);
	}
	virtual void initDensity(const int num_tet, VSX &rho){
	  TRACE_FUN();
	  const vector<double> &r = tetmesh->material()._rho;
	  assert_eq(num_tet, r.size());
	  rho.resize(num_tet);
	  for (int i = 0; i < num_tet; ++i)	rho[i] = r[i];
	}
	virtual vector<double> getDensityResult()const{
	  return tetmesh->material()._rho;
	}
	virtual void initAllVariables(VSX &x)const{x = CASADI::connect(E,v);}
	virtual void getInitValue(VectorXd &init_x)const;
	virtual SXMatrix assembleObjMatrix();

	virtual VSX getE()const{return E;}
	virtual VSX getV()const{return v;}

	virtual void init_Ev(const int num_tet, VSX &E, VSX &v){
	  assert_eq(E.size(),v.size());
	  if (E.size() != num_tet){
		E = CASADI::makeSymbolic(num_tet, "E");
		v = CASADI::makeSymbolic(num_tet, "v");
	  }
	}
	virtual vector<double> getYoungE()const{
	  const int num_tet = tetmesh->tets().size();
	  assert_ge(rlst.size(),num_tet);
	  return vector<double>(rlst.begin(),rlst.begin()+num_tet);
	}
	virtual vector<double> getPoissonV()const{
	  const int num_tet = tetmesh->tets().size();
	  assert_ge(rlst.size(),num_tet*2);
	  return vector<double>(rlst.begin()+num_tet,rlst.begin()+2*num_tet);
	}

	template<class T>
	T getShareG(const T&E, const T&v)const{return E/(2.0*(1.0+v));}
	template<class T>
	T getLameL(const T&E, const T&v)const{return (E*v)/((1.0+v)*(1.0-2.0*v));}
	template<class T>
	vector<T> getShareG(const vector<T>&E, const vector<T>&v)const{
	  assert_eq(E.size(),v.size());
	  vector<T> G(E.size());
	  for (size_t i = 0; i < E.size(); ++i)
		G[i] = getShareG(E[i],v[i]);
	  return G;
	}
	template<class T>
	vector<T> getLameL(const vector<T>&E, const vector<T>&v)const{
	  assert_eq(E.size(),v.size());
	  vector<T> L(E.size());
	  for (size_t i = 0; i < E.size(); ++i)
		L[i] = getLameL(E[i],v[i]);
	  return L;
	}
	
  protected:
	VSX E,v;
  };

  // fix v, and fit only E by by diagonalizing K.
  // W^t*K(E)*W = \Lambda
  class MaterialFitting_EV_Diag_K:public MaterialFitting_EV{
	
  protected:
	void init_Ev(const int num_tet, VSX &E, VSX &v){
	  assert_eq(E.size(),v.size());
	  if (E.size() != num_tet){
		E = CASADI::makeSymbolic(num_tet, "E");
		const vector<double> vv = getPoissonV();
		assert_eq(vv.size(),E.size());
		v.resize(vv.size());
		for (int i = 0; i < vv.size(); ++i)
		  v[i] = vv[i];
	  }
	}
	void initAllVariables(VSX &x)const{x = E;}
	void getInitValue(VectorXd &init_x)const{

	  const vector<double> &g = tetmesh->material()._G;
	  const vector<double> &la = tetmesh->material()._lambda;
	  const vector<double> E = MaterialFitting::getYoungE(g,la);
	  const int num_tet = tetmesh->tets().size();
	  init_x.resize(num_tet);
	  for (int i = 0; i < num_tet; ++i)
		init_x[i] = E[i];
	}
	vector<double> getPoissonV()const{
	  const vector<double> &g = tetmesh->material()._G;
	  const vector<double> &la = tetmesh->material()._lambda;
	  return MaterialFitting::getPoissonV(g,la);
	}
  };

  // fix v, and fit only E by by using modal analysis.
  // K(E)*W = M0*W*\Lambda
  class MaterialFitting_EV_MA_K:public MaterialFitting_EV_Diag_K{
  protected:
	SXMatrix assembleObjMatrix();
  };

  // fit density and scale of the basis by diagonalizing M.
  // ||W^t*M(\rho)*W-S||+||S-I|| where S is a diagonal matrix.
  class MaterialFitting_Diag_M_ScaleW: public MaterialFitting_Diag_M{

  public:
	void setWLambda(const MatrixXd &eigen_W,const VectorXd &eigen_lambda){
	  MaterialFitting_Diag_M::setWLambda(eigen_W, eigen_lambda);
	  scale_W = CASADI::makeSymbolic(eigen_lambda.size(),"s");
	}

  protected:
	SXMatrix assembleObjMatrix();
	void initAllVariables(VSX &x)const{
	  assert(tetmesh);
	  assert_eq(rho.size(), tetmesh->tets().size());
	  assert_eq(W.size2(), scale_W.size());
	  x = CASADI::connect(rho,scale_W);
	}
	void getInitValue(VectorXd &init_x)const;

  private:
	VSX scale_W;
  };

}//end of namespace

#endif /*_MATERIALFITTING_H_*/
