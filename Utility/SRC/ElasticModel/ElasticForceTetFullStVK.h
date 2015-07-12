#ifndef _ELASTICFORCETETFULLSTVK_H_
#define _ELASTICFORCETETFULLSTVK_H_

#include <boost/shared_ptr.hpp>
#include <ElasticForceTetFull.h>

namespace UTILITY{

  struct HessianOpHcIJ{
  public:
	HessianOpHcIJ(VectorXd& out,const VectorXd& I,const VectorXd& J):_out(out),_I(I),_J(J){out.setZero();}
	virtual void operator()(const int& depth,const int& i,const int& j,const Matrix3d& m){
	  _out.block<3,1>(i*3,0) += (m*_I[depth])*_J.block<3,1>(j*3,0);
	}
	VectorXd& _out;
	const VectorXd& _I;
	const VectorXd& _J;
  };
  
  /**
   * @class ElasticForceTetFullStVK full force using StVK constitutive model.
   * @bug curently, return -f(x) and -K(x)
   * @bug the eigenvalues of K(x) are all minus.
   */
  class ElasticForceTetFullStVK: public ElasticForceTetFull{
	
  public:
	ElasticForceTetFullStVK():ElasticForceTetFull(){}
	ElasticForceTetFullStVK(pTetMesh_const vol_mesh):ElasticForceTetFull(vol_mesh){}

	double energy(const VectorXd &X);
	void forceHessian(HessianOpHcIJ& oper,const VectorXd& X,const VectorXd& V);
    void forceHessianHcIJ(VectorXd& HIJ,const VectorXd& I,const VectorXd& J,const VectorXd& X,const VectorXd& V){
	  HessianOpHcIJ op(HIJ,I,J);
	  forceHessian(op,X,V);
	}

  public:
	void force_tet(mat3x4& f,const int& i,const VectorXd& X);
	void forceDerivX_tet(TetDF &df, const int& i, const VectorXd& X);
	void forceDerivXdX_tet(mat3x4& dfdX,const int& i,const VectorXd& dx,const VectorXd& X);

  protected:
	void computeTetForces(VecMat3x4 &tet_forces, const VectorXd &X){
	  tet_forces.resize(_vol_mesh->tets().size()); 
	  for (size_t i = 0; i < _vol_mesh->tets().size(); ++i){
		force_tet(tet_forces[i],i,X);
	  }
	}
	void computeTetForceDerivX(std::vector<TetDF>&df,const VectorXd &X){
	  df.resize(_vol_mesh->tets().size());
	  for (size_t i = 0; i < _vol_mesh->tets().size(); ++i){
		forceDerivX_tet(df[i],i,X);
	  }
	}
	void computeTetForceDerivXdX(VecMat3x4 &tet_kdx,const VectorXd &dx,const VectorXd &X){
	  tet_kdx.resize(_vol_mesh->tets().size());
	  for (size_t i = 0; i < _vol_mesh->tets().size(); ++i){
		forceDerivXdX_tet(tet_kdx[i],i,dx,X);
	  }
	}

  protected: 
	void W(double& W,const Matrix3d& F,const double& G,const double& lambda,const double& V) const;
	void PF(Matrix3d& P,const Matrix3d& F,const double& G,const double& lambda) const;
	void dPF(Matrix3d& deriv,const Matrix3d& dF,const Matrix3d& F,const double& G,const double& lambda) const;
	void ddPF(Matrix3d& deriv,const Matrix3d& dF,const Matrix3d& dF2,const Matrix3d& F,const double&G,const double&lambda,const int&i)const;
  };
  
  typedef boost::shared_ptr<ElasticForceTetFullStVK> pElasticForceTetFullStVK;
}//end of namespace

#endif /*_ELASTICFORCETETFULLSTVK_H_*/
