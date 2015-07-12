#ifndef _FULLSIMULATOR_H_
#define _FULLSIMULATOR_H_

#include <string>
#include <assertext.h>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <FullElasticModel.h>
#include <JsonFilePaser.h>
using namespace std;
using namespace Eigen;

namespace SIMULATOR{
  
  /**
   * @class BaseFullSim base class for making simulation in full space, child
   * class would implement different integration approaches (implicit,
   * semi-implicit, and explicit integration method), and different constraint
   * approaches (Lagrangian or Penalty method.).
   * 
   */
  class BaseFullSim{

  public: 
	BaseFullSim(){
	  h = 0.01;
	  alpha_k = 0.01;
	  alpha_m = 0.01;
	  correct_initialized = false;
	}
	BaseFullSim(pBaseFullModel def_model){
	  this->def_model = def_model;
	  h = 0.01;
	  alpha_k = 0.01;
	  alpha_m = 0.01;
	  correct_initialized = false;
	}
	virtual void setElasticModel(pBaseFullModel def_model){
	  this->def_model = def_model;
	}

	virtual bool init(const string init_filename){

	  UTILITY::JsonFilePaser jsonf;
	  if (!jsonf.open(init_filename)){
		ERROR_LOG("failed to open" << init_filename);
		return false;
	  }
	  
	  correct_initialized = jsonf.read("h",h);               assert_gt(h,0);
	  correct_initialized &= jsonf.read("alpha_m",alpha_m);  assert_ge(alpha_m,0);
	  correct_initialized &= jsonf.read("alpha_k",alpha_k);  assert_ge(alpha_k,0);
	  correct_initialized &= def_model->init(init_filename);
	  if (correct_initialized) BaseFullSim::setInitValue();
	  ERROR_LOG_COND("failed to initialze. ", correct_initialized);
	  return correct_initialized;
	}
	virtual bool prepare(){
	  const bool succ = def_model->prepare();
	  this->reset();
	  return succ;
	}
	bool initialized()const{
	  return correct_initialized;
	}
	void setU0(const VectorXd &u0){u = u0;}
	void setV0(const VectorXd &v0){v = v0;}
	virtual void setTimeStep(const double h){
	  assert_gt (h ,0.0f);
	  this->h = h;  
	}
	virtual void setDampings(const double alpha_k,const double alpha_m){
	  assert_ge (alpha_k , 0.0f);
	  assert_ge (alpha_m , 0.0f);
	  this->alpha_k = alpha_k;
	  this->alpha_m = alpha_m;  
	}
	virtual void reset(){
	  const int r = getDim();
	  v.resize(r);
	  v.setZero();
	  u.resize(r);
	  u.setZero();
	  fext.resize(r);
	  fext.setZero();
	}

	/**
	 * set constraints,where C_triplet is the triplet format of the constrain
	 * matrix, uc is displacements of the the constrained nodes (or the
	 * barycenters of constrained groups), where Cu = uc.
	 * C_rows, C_cols: the rows and cols of C.
	 */
	virtual void setConM(const VecT &C_triplet,const int C_rows,const int C_cols) = 0;
	virtual void setUc(const VectorXd &uc){
	  this->uc = uc;
	}
	virtual void removeAllCon(){
	  uc.resize(0);
	}
	const VectorXd &getPosCon()const{
	  return uc;
	}

	// set external forces in fullspace.
	void setExtForce(const VectorXd &full_fext){
	  this->fext = full_fext;
	}
	virtual void setExtForceOfNode(const double force[3],int vertex_id){
	  assert_ge(vertex_id,0);
	  assert_gt(fext.size(),(vertex_id*3+2));
	  fext[vertex_id*3+0] = force[0];
	  fext[vertex_id*3+1] = force[1];
	  fext[vertex_id*3+2] = force[2];
	}
	void setExtForceForAllNodes(const double force_x,const double force_y,const double force_z){
	  const int n = getDim()/3;
	  fext.resize(3*n);
	  for (int vertex_id = 0; vertex_id < n; ++vertex_id){
		fext[vertex_id*3+0] = force_x;
		fext[vertex_id*3+1] = force_y;
		fext[vertex_id*3+2] = force_z;    
	  }
	}

	// simulation
	virtual bool forward() = 0;
	
	// get results
	const VectorXd &getU()const{
	  return u;
	}
	const VectorXd &getV()const{
	  return v;
	}
	VectorXd &getU(){
	  return u;
	}
	VectorXd &getV(){
	  return v;
	}
	int getDim()const{
	  if (def_model){
		return def_model->dimension();
	  }
	  return 0;
	}
	double getTimestep()const{
	  return h;
	}
	double getAlphaK()const{
	  return alpha_k;
	}
	double getAlphaM()const{
	  return alpha_m;
	}

  protected:
	void setInitValue(){
	  const int n3 = this->getDim();
	  v.resize(n3);
	  v.setZero();
	  u.resize(n3);
	  u.setZero();
	  fext.resize(n3);
	  fext.setZero();
	}
	
  protected:
	pBaseFullModel def_model; // deformation data model, calculate F(u), K(u).
	bool correct_initialized; // the forward() should only be called when this
							  // class is correctly initialized.
	double h;
	double alpha_k;
	double alpha_m;

	VectorXd v; // velocity.
	VectorXd u; // displacements.
	VectorXd fext;// external forces.

	VectorXd uc;// displacements of the the constrained nodes (or the
				// barycenters of the constrained groups) with respect to the
				// rest shape.
  };
  typedef boost::shared_ptr<BaseFullSim> pBaseFullSim;

  /**
   * @class LagImpFullSim simulator in full space using Lagrangian constraints
   * and implicit integration scheme.
   * @see Doc/IntegrationMethods.pdf
   */
  class LagImpFullSim: public BaseFullSim{
	
  public:
	LagImpFullSim():BaseFullSim(){}
	LagImpFullSim(pBaseFullModel def_model):BaseFullSim(def_model){}
	bool prepare(){
	  bool succ = false;
	  if (BaseFullSim::prepare()){
		succ = def_model->evaluateM(M);
	  }
	  C_Ct_triplet.clear();
	  resetM_triplet();
	  return succ;
	}
	void setTimeStep(const double h){
	  BaseFullSim::setTimeStep(h);
	  resetM_triplet();
	}
	void setDampings(const double alpha_k,const double alpha_m){
	  BaseFullSim::setDampings(alpha_k, alpha_m);
	  resetM_triplet();
	}
	void setConM(const VecT &C_triplet,const int C_rows,const int C_cols);
	void removeAllCon(){
	  BaseFullSim::removeAllCon();
	  C_Ct_triplet.clear();
	  C.resize(0,0);
	}
	bool forward();

  protected:
	bool assembleA();
	bool assembleB();
	bool resetM_triplet();
	bool resetK_triplet();
	void resetA_triplet();
	void copyTriplet(VecT &Full_triplet, const VecT &sub_triplet, const int start)const;
	
  protected:
	VectorXd f; // internal forces.
	SparseMatrix<double> M; // mass matrix.
	SparseMatrix<double> C; // constraint matrix.
	
	/**
	 * right hand side and left hand side for the linear system at each time
	 * step to sovle for the integration, where
	 *     |H C^t|
	 * A = |C  O |, b = M*v + h*(fext-f),
	 * and H = (1+h*alpha_m)*M + h*(h+alpha_k)*K(u),
	 * C is the constraint matrix.
	 */
	SparseMatrix<double> A; 
	VectorXd b;  

	/**
	 * Inorder to fastly assemble sparse matrix A, we record the triplets for
	 * each submatrix in A_triplet, then generate A from A_triplet each step.
	 * A_triplet is consisted of three parts: 
	 * 1. scaled mass matrix: update when intialization or time step is changed.
	 * 2. constraint matrix and its transpose, update when constraints is changed.
	 * 3. scaled stiffness matrix, update at each step.
	 */
	VecT C_Ct_triplet;     // C and C^t.
	VecT Scaled_M_triplet; // (1+alpha_m)*M.
	VecT Scaled_K_triplet; // h*(h+alpha_k)*K(u).
	VecT A_triplet; // [scaled_M, C, C^T, scaled_K].
  };
  
  /**
   * @class PenStaticFullSim quasi-static simulator with penalty constraints, 
   * and we solve the nonlinear equation using newton method.
   * @see Doc/StaticDeformation.pdf
   */
  class PenStaticFullSim: public BaseFullSim{
	
  public:
	PenStaticFullSim():BaseFullSim(){
	  lambda = 100.0f;
	  max_it = 5;
	  tolerance = 0.1f;
	}
	PenStaticFullSim(pBaseFullModel def_model):BaseFullSim(def_model){
	  lambda = 100.0f;
	  max_it = 5;
	  tolerance = 0.1f;
	}
	bool init(const string init_filename);
	void setConM(const VecT &C_triplet,const int C_rows,const int C_cols);
	void setUc(const VectorXd &uc);
	void removeAllCon();
	virtual bool forward();

  protected:
	const VectorXd &grad(const VectorXd &u);
	const SparseMatrix<double> &jac(const VectorXd &u);
	
  protected:
	double lambda; // penalty for con.
	SparseMatrix<double> C; // constraint matrix.
	SparseMatrix<double> lambda_CtC; // C^t*C
	VectorXd lambda_CtUc; // C^t*Uc

	VectorXd g;
	SparseMatrix<double> J;
	int max_it;
	double tolerance;
  };


  /**
   * @class PenStaticFullSim semi-implicit simulator with penalty constraints.
   * @see Doc/IntegrationMethods.pdf
   * @note currently, we only use the mass-damping, and ignore the stiffness damping.
   */
  class PenSemiImpFullSim: public PenStaticFullSim{
	
  public:
    PenSemiImpFullSim(): PenStaticFullSim(){}
	PenSemiImpFullSim(pBaseFullModel def_model):PenStaticFullSim(def_model){}
	bool prepare();
	bool forward();

  private:
	DiagonalMatrix<double,-1> h_M_inv; // h*(M^{-1})
  };
  
}//end of namespace

#endif /* _FULLSIMULATOR_H_ */
