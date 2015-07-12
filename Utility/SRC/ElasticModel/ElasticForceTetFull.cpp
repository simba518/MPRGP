#include <stdio.h>
#include <Timer.h>
#include "ElasticForceTetFull.h"
using namespace UTILITY;

void ElasticForceTetFull::force(const VectorXd &X, VectorXd &out){

  out.resize(X.size());
  out.setZero();

  // compute elastic forces for each nodes of each tet.
  this->computeTetForces(_tet_forces,X);

  // assemble
  for(int i=0;i<(int)_vol_mesh->tets().size();i++){

	const Vector4i &e=_vol_mesh->tets()[i];
	const mat3x4 &f=_tet_forces[i];
	out.block<3,1>(e[0]*3,0)+=f.col(0);
	out.block<3,1>(e[1]*3,0)+=f.col(1);
	out.block<3,1>(e[2]*3,0)+=f.col(2);
	out.block<3,1>(e[3]*3,0)+=f.col(3);
  }
}

const SparseMatrix<double> &ElasticForceTetFull::K(const VectorXd &X){

  this->computeTetForceDerivX(_tet_k,X);
  for (int i = 0; i < _Kx.nonZeros(); ++i){
	_Kx.valuePtr()[i] = 0.0f;
  }
  for(int i=0,etr=0;i<(int)_vol_mesh->tets().size();i++){
	const TetDF &df = _tet_k[i];
	for(int j=0;j<4;j++)
	  for(int k=0;k<4;k++){
		const Matrix3d& dfM=df.df[j][k];
		for(int r=0;r<3;r++)
		  for(int c=0;c<3;c++){
			_Kx.valuePtr()[_entries[etr++]]+=dfM(r,c);
		  }
	  }
  }
  return _Kx;
}

void ElasticForceTetFull::Kdx(const VectorXd&dx,const VectorXd&X,VectorXd&out){

  out.setZero();
  this->computeTetForceDerivXdX(_tet_kdx,dx,X);
  for(int i=0;i<(int)_vol_mesh->tets().size();i++){

	const Vector4i &e=_vol_mesh->tets()[i];
	const mat3x4 &dfdX = _tet_kdx[i];
		
	out.block<3,1>(e[0]*3,0)+=dfdX.col(0);
	out.block<3,1>(e[1]*3,0)+=dfdX.col(1);
	out.block<3,1>(e[2]*3,0)+=dfdX.col(2);
	out.block<3,1>(e[3]*3,0)+=dfdX.col(3);
  }
}

bool ElasticForceTetFull::prepare(){

  // compute inv(Dm), and dF.
  _def_grad.prepare(_vol_mesh);

  // compute volume
  _volume.resize(_vol_mesh->tets().size());
  for(int i=0;i<(int)_vol_mesh->tets().size();i++) {
	_volume[i]=_vol_mesh->volume(i);
  }
  initKEntry(_Kx, _entries);
  return true;
}

void ElasticForceTetFull::initKEntry(SparseMatrix<double> &Kx,std::vector<int> &entries)const{

  //build Kx with all entries are zero.
  const int dim=(int)_vol_mesh->nodes().size()*3;
  Kx.resize(dim,dim);
  std::vector<Eigen::Triplet<double> > trips;
  for(int i=0;i<(int)_vol_mesh->tets().size();i++){
	const Vector4i& e=_vol_mesh->tets()[i];
	for(int j=0;j<4;j++)
	  for(int k=0;k<4;k++){
		for(int row=0;row<3;row++)
		  for(int col=0;col<3;col++)
			trips.push_back(Eigen::Triplet<double>((int)e[j]*3+row,(int)e[k]*3+col,0.0f));
	  }
  }
  Kx.setFromTriplets(trips.begin(),trips.end());
  Kx.makeCompressed();

  //build reference _entries
  entries.clear();
  for(int i=0;i<(int)_vol_mesh->tets().size();i++){
	const Vector4i& e=_vol_mesh->tets()[i];
	for(int j=0;j<4;j++)
	  for(int k=0;k<4;k++){
		const Vector2i base(e[j]*3,e[k]*3);
		for(int r=0;r<3;r++)
		  for(int c=0;c<3;c++){
			const double& v=Kx.coeffRef((int)base.x()+r,(int)base.y()+c);
			entries.push_back((&v-Kx.valuePtr()));
		  }
	  }
  }
}
