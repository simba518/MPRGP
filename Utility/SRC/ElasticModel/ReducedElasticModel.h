#ifndef _REDUCEDELASTICMODEL_H_
#define _REDUCEDELASTICMODEL_H_

#include <string>
#include <vector>
#include <assertext.h>
#include <ElasticForceTetFullStVK.h>
using namespace std;
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class ReducedElasticModel base class for reduced elastic data model.
   * 
   */
  class ReducedElasticModel{
	
  public:
	// initialize from the ini file
	virtual bool init(const string init_filename);

	virtual void setModalBasis(const MatrixXd &B){
	  this->B = B;
	}

	virtual void setVolMesh(pTetMesh_const vol_mesh){}

	virtual bool prepare(){}

	// compute the internal forces in subspace
	virtual bool evaluateF(const VectorXd &reduced_u,VectorXd &f) = 0;

	// compute the stiffness matrix in subspace
	virtual bool evaluateK(const VectorXd &reduced_u,MatrixXd &K) = 0;

	// compute the internal forces in subspace
	virtual bool evaluateFK(const VectorXd &reduced_u,VectorXd &f,MatrixXd &K){
	  bool succ = evaluateF(reduced_u, f);
	  succ &= evaluateK(reduced_u, K);
	  return succ;
	}

	// return the reduced mass matrix: U^t*M*U
	virtual const MatrixXd &getReducedMassM()const{
	  assert_eq(M.rows(), M.cols());
	  assert_eq(M.rows(), reducedDim());
	  return M;
	}

	// convert the displacements in subsapce to fullspace
	virtual bool computeFullDisp(const VectorXd &reduced_u,VectorXd &full_u){
	  assert_eq(reduced_u.size(), B.cols());
	  full_u = B*reduced_u;
	  return true;
	}

	// dimension of the reduced subspace
	int reducedDim()const{
	  return B.cols();
	}

	// dimension of the full space.
	int fullDim()const{
	  return B.rows();
	}

	// get basis of the subspace
	const MatrixXd &getModalBasis()const{
	  return B;
	}
	
  protected:
	MatrixXd B;
	MatrixXd M;
  };
  typedef boost::shared_ptr<ReducedElasticModel> pReducedElasticModel;

  /**
   * @class DirectReductionElasticModel Given B, reduce the F_full and K_full directly, where,
   * f = B^T*F_full(x0+B*q)
   * K = B^T*K_full(x0+B*q)*B
   */
  class DirectReductionElasticModel:public ReducedElasticModel,public ElasticForceTetFullStVK{
	
  public:
	DirectReductionElasticModel():ElasticForceTetFullStVK(){}
	DirectReductionElasticModel(pTetMesh_const tet):ElasticForceTetFullStVK(tet){}
	void setVolMesh(pTetMesh_const vol_mesh){
	  ElasticForceTetFullStVK::setVolMesh(vol_mesh);
	}
	bool prepare();
	bool evaluateF(const VectorXd &reduced_u,VectorXd &f){
	  assert_eq(rest_shape.size(), B.rows());
	  assert_eq(reduced_u.size(), B.cols());
	  const VectorXd x = rest_shape+B*reduced_u;
	  static VectorXd f_full;
	  ElasticForceTetFullStVK::force(x,f_full);
	  f = B.transpose()*f_full;
	  f *= -1.0f;
	  return true;
	}
	bool evaluateK(const VectorXd &reduced_u,MatrixXd &K){
	  const VectorXd x = rest_shape+B*reduced_u;
	  K = B.transpose()*(ElasticForceTetFullStVK::K(x)*B);
	  K *= -1.0f;
	  return true;
	}
	bool evaluateFK(const VectorXd &reduced_u,VectorXd &f,MatrixXd &K){
	  assert_eq(rest_shape.size(), B.rows());
	  assert_eq(reduced_u.size(), B.cols());
	  const VectorXd x = rest_shape+B*reduced_u;
	  static VectorXd f_full;
	  ElasticForceTetFullStVK::force(x,f_full);
	  f = B.transpose()*f_full;
	  K = B.transpose()*(ElasticForceTetFullStVK::K(x)*B);
	  f *= -1.0f;
	  K *= -1.0f;
	  return true;
	}

  private:
	VectorXd rest_shape;
  };

  /**
   * @class CubaturedElasticModel reduce the F_full and K_full using cubatured samples, where,
   * f = sum wi*Bi^t*fi(q), for i in S.
   * K = sum wi*Bi^t*Ki(q)*Bi, for i in S.
   * @see 1. An efficient construction of reduced deformable objects, asia 2013.
   * 2. Optimizing Cubature for Efficient Integration of Subspace Deformations. asia 2008.
   */
  class CubaturedElasticModel:public ReducedElasticModel,public ElasticForceTetFullStVK{
	
  public:
	CubaturedElasticModel():ElasticForceTetFullStVK(){}
	CubaturedElasticModel(pTetMesh_const tet):ElasticForceTetFullStVK(tet){}
	void setVolMesh(pTetMesh_const vol_mesh){
	  ElasticForceTetFullStVK::setVolMesh(vol_mesh);
	}
	bool prepare();
	void setCubature(const vector<double> &w, const vector<int> &S){
	  assert_eq(w.size(), S.size());
	  weights = w;
	  sampledTets = S;
	}
	bool init(const string init_filename);
	bool evaluateF(const VectorXd &reduced_u,VectorXd &f);
	bool evaluateK(const VectorXd &reduced_u,MatrixXd &K);

  private:
	vector<double> weights;
	vector<int> sampledTets;
	VectorXd rest_shape;
  };
  typedef boost::shared_ptr<CubaturedElasticModel> pCubaturedElasticModel;
  
}//end of namespace

#endif /*_REDUCEDELASTICMODEL_H_*/
