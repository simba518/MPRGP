#ifndef _DEFGRADTET_H_
#define _DEFGRADTET_H_

#include <boost/shared_ptr.hpp>
#include <TetMesh.h>

namespace UTILITY{

  struct derivF{
	Matrix3d _dF[4][3];
  };

  typedef std::vector<Matrix3d,Eigen::aligned_allocator<Matrix3d> > VectorM3x3;
  
  /**
   * @class DefGradTet deformation gradient of tet mesh.
   * 
   */
  class DefGradTet{
	
  public:
	void prepare(pTetMesh_const vol_mesh){
	  this->_vol_mesh = vol_mesh;
	  VectorXd nodes;
	  _vol_mesh->nodes(nodes);
	  prepare(nodes, _vol_mesh->tets());
	}
	void evalF(VectorM3x3 &F, const VectorXd &X)const{
	  F.resize((int)_vol_mesh->tets().size());
	  for (int i=0;i<(int)_vol_mesh->tets().size();i++) {
		evalFe(F[i],X,i);
	  }
	}
	void evalFe(Matrix3d &Fe, const VectorXd &X, const int& i)const{
	  const Vector4i& e=_vol_mesh->tets()[i];
	  deformationGradEntry(Fe,_invDm[i],
						   X.block<3,1>(e[0]*3,0),X.block<3,1>(e[1]*3,0),
						   X.block<3,1>(e[2]*3,0),X.block<3,1>(e[3]*3,0));
	}
	const VectorM3x3 &invDm()const{
	  return _invDm;
	}
	const std::vector<derivF> &dF()const{
	  return _dF;
	}

  protected:
	void prepare(const VectorXd &nodes,const VVec4i &tets){

	  _invDm.resize( tets.size() );
	  _dF.resize(tets.size());
	  for(int i=0;i<(int)tets.size();i++) {

		const Vector4i& e=tets[i];
		Matrix3d Dm;
		Dm.col(0)=nodes.block<3,1>(e[0]*3,0)-nodes.block<3,1>(e[3]*3,0);
		Dm.col(1)=nodes.block<3,1>(e[1]*3,0)-nodes.block<3,1>(e[3]*3,0);
		Dm.col(2)=nodes.block<3,1>(e[2]*3,0)-nodes.block<3,1>(e[3]*3,0);
		_invDm[i]=Dm.inverse();
		
		for(int j=0;j<4;j++)
		  for(int k=0;k<3;k++)
			calculateDF(_dF[i]._dF[j][k],_invDm[i],j,k);
	  }
	}
	void calculateDF(Matrix3d& dF,const Matrix3d& invDm,const int& i,const int& j) const{
	  Vector3d X[4]={Vector3d::Zero(),Vector3d::Zero(),Vector3d::Zero(),Vector3d::Zero()};
	  X[i]=Vector3d::Unit(j);
	  deformationGradEntry(dF,invDm,X[0],X[1],X[2],X[3]);
	}
	template <typename T>
	void deformationGradEntry(Matrix3d& ret,const Matrix3d& invDm,const T& a,
							  const T& b,const T& c,const T& d) const{

	  Matrix3d Ds;
	  Ds.col(0)=a-d;
	  Ds.col(1)=b-d;
	  Ds.col(2)=c-d;
	  ret=Ds*invDm;
	}

  private:
	pTetMesh_const _vol_mesh;
	VectorM3x3 _invDm;
	std::vector<derivF> _dF;
  };
  
  typedef boost::shared_ptr<DefGradTet> pDefGradTet;
  
}//end of namespace

#endif /*_DEFGRADTET_H_*/
