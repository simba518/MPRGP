#ifndef _COMPUTESTIFFNESSMAT_H_
#define _COMPUTESTIFFNESSMAT_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Sparse>
#include <DefGradTet.h>
#include <CASADITools.h>
using namespace UTILITY;
using CasADi::SX;
using CasADi::SXMatrix;

namespace ELASTIC_OPT{

  struct TetDF{ // dF/dX for one tetrahedron.
    SXMatrix df[4][4];
  };
  
  /**
   * @class compute stiffness matrix using casadi.
   * 
   */
  class ComputeStiffnessMat{
	
  public:
	ComputeStiffnessMat(pTetMesh_const vol_mesh):_vol_mesh(vol_mesh){
	  setMaterial(_vol_mesh->material()._G,_vol_mesh->material()._lambda);
	}
	template<class T> 
	void setMaterial(const std::vector<T> &G,const std::vector<T> &lambda){
	  assert_eq(G.size(), lambda.size());
	  _G.resize(G.size());
	  _lambda.resize(lambda.size());
	  for (int i = 0; i < G.size(); ++i){
		_G[i] = G[i];
		_lambda[i] = lambda[i];
	  }
	}
	bool prepare();
	const SXMatrix &K(const VectorXd &X);

  protected:
	void computeTetForceDerivX(std::vector<TetDF>&df,const VectorXd &X);
	void forceDerivX_tet(TetDF &df, const int& i, const VectorXd& X);
	void dPF(SXMatrix& deriv,const Matrix3d& dF,const Matrix3d& F,const SX& G,const SX& lambda) const;

  protected:
	pTetMesh_const _vol_mesh;
	DefGradTet _def_grad;
	std::vector<double> _volume;
	std::vector<TetDF> _tet_k;
	SXMatrix _Kx; // K(x)
	std::vector<SX> _G; // Shear modulus.
	std::vector<SX> _lambda; // Lamé's first parameter.
  };
  
  typedef boost::shared_ptr<ComputeStiffnessMat> pComputeStiffnessMat;
  
}//end of namespace

#endif /* _COMPUTESTIFFNESSMAT_H_ */
