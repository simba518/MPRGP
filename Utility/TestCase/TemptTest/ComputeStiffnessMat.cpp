#include <stdio.h>
#include <map>
#include <Timer.h>
#include "ComputeStiffnessMat.h"
using namespace std;
using namespace UTILITY;
using namespace ELASTIC_OPT;

bool ComputeStiffnessMat::prepare(){

  // compute inv(Dm), and dF.
  _def_grad.prepare(_vol_mesh);

  // compute volume
  _volume.resize(_vol_mesh->tets().size());
  for(int i=0;i<(int)_vol_mesh->tets().size();i++) {
	_volume[i]=_vol_mesh->volume(i);
  }
  return true;
}

const SXMatrix &ComputeStiffnessMat::K(const VectorXd &X){

  const int dim=(int)_vol_mesh->nodes().size()*3;
  _Kx.setZero();
  _Kx.resize(dim,dim);
  this->computeTetForceDerivX(_tet_k,X);
  _Kx.reserve(_tet_k.size()*16);

  std::vector<int> row,col;
  std::vector<SX> d;
  row.reserve(_tet_k.size()*16*9);
  col.reserve(_tet_k.size()*16*9);
  d.reserve(_tet_k.size()*16*9);
  map<pair<int,int>, int> record;

  for(int i=0,etr=0;i<(int)_vol_mesh->tets().size();i++){

	const Vector4i& e=_vol_mesh->tets()[i];
	const TetDF &df = _tet_k[i];
	for(int j=0;j<4;j++)
	  for(int k=0;k<4;k++){
		const SXMatrix& dfM=df.df[j][k];
		for(int r=0;r<3;r++)
		  for(int c=0;c<3;c++){
			const int grow = e[j]*3+r;
			const int gcol = e[k]*3+c;
			const map<pair<int,int>,int>::iterator it = record.find(pair<int,int>(grow,gcol));
			if (it != record.end()){
			  const int index = it->second;
			  assert_in(index,0,d.size());
			  d[index] += dfM.elem(r,c);
			}else{
			  row.push_back(grow);
			  col.push_back(gcol);
			  d.push_back(dfM.elem(r,c));
			  record[pair<int,int>(grow,gcol)] = d.size()-1;
			}
		  }
	  }
  }
  _Kx = SXMatrix::sparse(row, col, d, dim, dim);
  return _Kx;
}

void ComputeStiffnessMat::computeTetForceDerivX(std::vector<TetDF>&df,const VectorXd&X){

  df.resize(_vol_mesh->tets().size());
  for (size_t i = 0; i < _vol_mesh->tets().size(); ++i){
	forceDerivX_tet(df[i],i,X);
  }
}

void ComputeStiffnessMat::forceDerivX_tet(TetDF &df, const int& i, const VectorXd& X){

  for (int k = 0; k < 4; ++k){
    for (int j = 0; j < 4; ++j)
	  df.df[k][j].resize(3,3);
  }

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);

  const derivF& dF=_def_grad.dF()[i];
  for(int x=0;x<4;x++)
	for(int d=0;d<3;d++){
	  const Matrix3d& dFxd=dF._dF[x][d];
	  static SXMatrix deriv;
	  dPF(deriv,dFxd,F,_G[i], _lambda[i]);
	  const MatrixXd dt = _def_grad.invDm()[i].transpose();
	  deriv = -_volume[i]*(deriv.mul( CASADI::convert(dt) ));

	  for(int fx=0;fx<3;fx++)
	  	for(int fd=0;fd<3;fd++){
	  	  df.df[fx][x].elem(fd,d) = deriv.elem(fd,fx);
		}
	}

  df.df[3][0]=-df.df[0][0]-df.df[1][0]-df.df[2][0];
  df.df[3][1]=-df.df[0][1]-df.df[1][1]-df.df[2][1];
  df.df[3][2]=-df.df[0][2]-df.df[1][2]-df.df[2][2];
  df.df[3][3]=-df.df[0][3]-df.df[1][3]-df.df[2][3];
}

void ComputeStiffnessMat::dPF(SXMatrix& deriv,const Matrix3d& dF,const Matrix3d& F,const SX& G,const SX& lambda) const{

  const Matrix3d  E=(F.transpose()*F-Matrix3d::Identity())*0.5f;
  const Matrix3d  derivE=(dF.transpose()*F+F.transpose()*dF)*0.5f;

  const MatrixXd M1 = 2.0f*(dF*E+F*derivE);
  const MatrixXd M2 = E.trace()*dF+derivE.trace()*F;
  const SXMatrix SM1 = CASADI::convert(M1);
  const SXMatrix SM2 = CASADI::convert(M2);

  deriv = G*SM1+lambda*SM2;
}
