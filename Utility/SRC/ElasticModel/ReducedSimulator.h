#ifndef _REDUCEDSIMULATOR_H_
#define _REDUCEDSIMULATOR_H_

#include <set>
#include <ReducedElasticModel.h>
using namespace std;

namespace SIMULATOR{
  
  /**
   * @class ReducedSimulator base class for subspace simulation.
   * 
   */
  class ReducedSimulator{
	
  public:
	// set
	ReducedSimulator(pReducedElasticModel m):model(m){
	  h = 1.0f;
	  alpha_k = 0.0f;
	  alpha_m = 0.0f;
	}
	virtual bool init(const string init_filename);
	virtual bool prepare(){
	  bool succ = false;
	  if (model) succ = model->prepare();
	  reset();
	  return succ;
	}
	void setQ0(const VectorXd &q0){q = q0;}
	void setV0(const VectorXd &v0){v = v0;}
	void setTimeStep(const double h){
	  assert_gt (h, 0.0f);
	  this->h = h;  
	}
	void setDampings(const double alpha_k,const double alpha_m){
	  assert_ge (alpha_k, 0.0f);
	  assert_ge (alpha_m, 0.0f);
	  this->alpha_k = alpha_k;
	  this->alpha_m = alpha_m;  
	}
	virtual void reset(){
	  const int r = reducedDim();
	  v.resize(r);
	  v.setZero();
	  q.resize(r);
	  q.setZero();
	  red_fext.resize(r);
	  red_fext.setZero();
	}

	// set constraints
	virtual void setConGroups(const vector<int> &con_nodes) = 0;
	virtual void setUc(const VectorXd &uc) = 0;
	virtual void removeAllCon() = 0;

	// set external forces in full space
	virtual void setExtForce(const VectorXd &full_ext){
	  assert(model);
	  const MatrixXd &B = model->getModalBasis();
	  assert_eq(B.rows(), full_ext.size());
	  red_fext = B.transpose()*full_ext;
	}

	// simulate
	virtual bool forward() = 0;
	
	// get
	const VectorXd &getQ()const{return q;}
	const VectorXd &getV()const{return v;}
	VectorXd &getQ(){return q;}
	VectorXd &getV(){return v;}
	double getTimestep()const{
	  return h;
	}
	double getAlphaK()const{
	  return alpha_k;
	}
	double getAlphaM()const{
	  return alpha_m;
	}
	int reducedDim()const{
	  return (NULL==model) ? 0:model->reducedDim();
	}
	int fullDim()const{
	  return (NULL==model) ? 0:model->fullDim();
	}

  protected:
	pReducedElasticModel model;
	double h;
	double alpha_k;
	double alpha_m;

	VectorXd v; // velocity.
	VectorXd q; // position in sub space.
	VectorXd red_fext; // external forces in subspace.
  };
  typedef boost::shared_ptr<ReducedSimulator> pReducedSimulator;

  /**
   * @class ReducedImpLogConSimulator implicit integration with Lagrangian constraints.
   * 
   */  
  class ReducedImpLogConSimulator:public ReducedSimulator{
	
  public:
	ReducedImpLogConSimulator(pReducedElasticModel m):ReducedSimulator(m){}
	void setConGroups(const vector<int> &con_nodes);
	void setUc(const VectorXd &uc){
	  this->uc = uc;
	}
	void removeAllCon();
	bool forward();

  private:
	MatrixXd A;
	VectorXd b;
	MatrixXd C;
	VectorXd uc;
  };

  /**
   * @class ReducedStaticPenConSimulator quasi-static simulator with penalty constraints, 
   * and we solve the nonlinear equation using newton method.
   * @see Doc/Latex/StaticDeformation.pdf
   */
  class ReducedStaticPenConSimulator:public ReducedSimulator{
	
  public:
	ReducedStaticPenConSimulator(pReducedElasticModel m):ReducedSimulator(m){
	  lambda = 100.0f;
	  max_it = 5;
	  tolerance = 0.1f;
	}
	bool init(const string init_filename);
	void setConGroups(const vector<int> &con_nodes);
	void setUc(const VectorXd &uc);
	void removeAllCon();
	bool forward();

  protected:
	const VectorXd &grad(const VectorXd &u);
	const MatrixXd &jac(const VectorXd &u);	
	
  private:
	double lambda; // penalty for con.
	MatrixXd C; // constraint matrix.
	MatrixXd lambda_CtC; // C^t*C
	VectorXd lambda_CtUc; // C^t*Uc

	VectorXd g;
	MatrixXd J;
	int max_it;
	double tolerance;
  };
  
}//end of namespace

#endif /*_REDUCEDSIMULATOR_H_*/
