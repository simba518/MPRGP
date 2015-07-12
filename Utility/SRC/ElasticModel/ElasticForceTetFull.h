#ifndef _ELASTICFORCETETFULL_H_
#define _ELASTICFORCETETFULL_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Sparse>
#include <DefGradTet.h>

namespace UTILITY{

  struct TetDF{ // dF/dX for one tetrahedron.
    Matrix3d df[4][4];
  };
  typedef Eigen::Matrix<double,3,4> mat3x4;
  typedef std::vector<mat3x4,Eigen::aligned_allocator<mat3x4> > VecMat3x4;
  
  /**
   * @class ElasticForceTetFull base class for full elastic forces of
   * tetrahedron mesh.
   * 
   */
  class ElasticForceTetFull{
	
  public:
	ElasticForceTetFull(){
	  
	}
	ElasticForceTetFull(pTetMesh_const vol_mesh):_vol_mesh(vol_mesh){}
	void setVolMesh(pTetMesh_const vol_mesh){
	  _vol_mesh = vol_mesh;
	}
	virtual bool prepare();
	virtual double energy(const VectorXd &X) = 0;
	virtual void force(const VectorXd &X, VectorXd &out);
	virtual const SparseMatrix<double> &K(const VectorXd &X);
	virtual void Kdx(const VectorXd &dx, const VectorXd &X, VectorXd &out);

  protected:
	virtual void computeTetForces(VecMat3x4 &tet_forces, const VectorXd &X)=0;
	virtual void computeTetForceDerivX(std::vector<TetDF>&df,const VectorXd &X)=0;
	virtual void computeTetForceDerivXdX(VecMat3x4 &tet_kdx,const VectorXd &dx,const VectorXd &X)=0;
	void initKEntry(SparseMatrix<double> &K,std::vector<int> &entries)const;
	void deformationGrad(Matrix3d& F,const int& i,const VectorXd& X,const VectorXd& V) const{
	  const Vector4i& e=_vol_mesh->tets()[i];
	  deformationGradEntry(F,_def_grad.invDm()[i],
						   X.block<3,1>(e[0]*3,0),X.block<3,1>(e[1]*3,0),
						   X.block<3,1>(e[2]*3,0),X.block<3,1>(e[3]*3,0));
	}

	template <typename T,typename T2,typename T3>
	void deformationGradEntry(Eigen::Matrix<T,3,3>& ret,const Eigen::Matrix<T2,3,3>& invDm,const T3& a,const T3& b,const T3& c,const T3& d) const{
	  Eigen::Matrix<T,3,3> Ds;
	  Ds.col(0)=a-d;
	  Ds.col(1)=b-d;
	  Ds.col(2)=c-d;
	  ret=Ds*invDm;
	}
	
	template <typename T>
	inline static double contract(const Eigen::Matrix<T,3,3> &A,const Eigen::Matrix<T,3,3> &B){
	  return A(0,0)*B(0,0)+A(0,1)*B(0,1)+A(0,2)*B(0,2)+
		A(1,0)*B(1,0)+A(1,1)*B(1,1)+A(1,2)*B(1,2)+
		A(2,0)*B(2,0)+A(2,1)*B(2,1)+A(2,2)*B(2,2);
	}

  protected:
	pTetMesh_const _vol_mesh;
	DefGradTet _def_grad;
	std::vector<double> _volume;
	VecMat3x4 _tet_forces;
	std::vector<TetDF> _tet_k;
	VecMat3x4 _tet_kdx;
	SparseMatrix<double> _Kx; // K(x)
	std::vector<int> _entries; 
  };
  
  typedef boost::shared_ptr<ElasticForceTetFull> pElasticForceTetFull;
  
}//end of namespace

#endif /*_ELASTICFORCETETFULL_H_*/
