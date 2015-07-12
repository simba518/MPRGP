#include <ElasticForceTetFullStVK.h>
using namespace UTILITY;

void ElasticForceTetFullStVK::forceHessian(HessianOpHcIJ& op,const VectorXd& X,const VectorXd& V){

  for(int i=0;i<(int)_vol_mesh->tets().size();i++){

	const Vector4i& tet=_vol_mesh->tets()[i];
	const derivF& dF=_def_grad.dF()[i];
	Matrix3d df[4][4];
	Matrix3d F;
	deformationGrad(F,i,X,V);
	for(int x2=0;x2<4;x2++)
	  for(int d2=0;d2<3;d2++){

		const Matrix3d& dFxd2=dF._dF[x2][d2];
		for(int x=0;x<4;x++)
		  for(int d=0;d<3;d++){
			const Matrix3d& dFxd=dF._dF[x][d];
			Matrix3d deriv;
			ddPF(deriv,dFxd,dFxd2,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i],i);
			deriv=-_volume[i]*deriv*_def_grad.invDm()[i].transpose();
			
			for(int fx=0;fx<3;fx++)
			  for(int fd=0;fd<3;fd++)
				df[fx][x](fd,d)=deriv(fd,fx);
		  }
		df[3][0]=-df[0][0]-df[1][0]-df[2][0];
		df[3][1]=-df[0][1]-df[1][1]-df[2][1];
		df[3][2]=-df[0][2]-df[1][2]-df[2][2];
		df[3][3]=-df[0][3]-df[1][3]-df[2][3];
		
		for(int fc=0;fc<4;fc++)
		  for(int nd=0;nd<4;nd++)
			op(tet[x2]*3+d2,tet[fc],tet[nd],df[fc][nd]);
	  }
  }
}

void ElasticForceTetFullStVK::force_tet(mat3x4& f,const int& i,const VectorXd& X){

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);
  Matrix3d  P;
  PF(P,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
  f.block<3,3>(0,0)=-_volume[i]*P*_def_grad.invDm()[i].transpose();
  f.col(3)=-(f.col(0)+f.col(1)+f.col(2));
}

void ElasticForceTetFullStVK::forceDerivX_tet(TetDF &df, const int& i, const VectorXd& X){

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);

  const derivF& dF=_def_grad.dF()[i];
  for(int x=0;x<4;x++)
	for(int d=0;d<3;d++){

	  const Matrix3d& dFxd=dF._dF[x][d];
	  Matrix3d  deriv;
	  dPF(deriv,dFxd,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
	  deriv=-_volume[i]*deriv*_def_grad.invDm()[i].transpose();

	  for(int fx=0;fx<3;fx++)
		for(int fd=0;fd<3;fd++)
		  df.df[fx][x](fd,d)=deriv(fd,fx);
	}

  df.df[3][0]=-df.df[0][0]-df.df[1][0]-df.df[2][0];
  df.df[3][1]=-df.df[0][1]-df.df[1][1]-df.df[2][1];
  df.df[3][2]=-df.df[0][2]-df.df[1][2]-df.df[2][2];
  df.df[3][3]=-df.df[0][3]-df.df[1][3]-df.df[2][3];
}

void ElasticForceTetFullStVK::forceDerivXdX_tet(mat3x4& dfdX,const int& i,const VectorXd& dx,const VectorXd& X){

  Matrix3d  F;
  _def_grad.evalFe(F,X,i);

  Matrix3d  dF;
  _def_grad.evalFe(dF,dx,i);
		
  Matrix3d  deriv;
  dPF(deriv,dF,F,_vol_mesh->material()._G[i],_vol_mesh->material()._lambda[i]);
		
  dfdX.block<3,3>(0,0)=-_volume[i]*deriv*_def_grad.invDm()[i];
  dfdX.col(3)=-(dfdX.col(0)+dfdX.col(1)+dfdX.col(2));
}

void ElasticForceTetFullStVK::W(double& W,const Matrix3d& F,const double& G,const double& lambda,const double& V) const{

  const Matrix3d  E=0.5f*(F.transpose()*F-Matrix3d::Identity());
  const double tr=E.trace();
  W=G*E.squaredNorm()+lambda*0.5f*tr*tr;
  W*=V;
}

void ElasticForceTetFullStVK::PF(Matrix3d& P,const Matrix3d& F,const double& G,const double& lambda) const{

  const Matrix3d  E=(F.transpose()*F-Matrix3d::Identity())*0.5f;
  P=F*E*(2.0f*G)+(E.trace()*lambda)*F;
}

void ElasticForceTetFullStVK::dPF(Matrix3d& deriv,const Matrix3d& dF,const Matrix3d& F,const double& G,const double& lambda) const{

  const Matrix3d  E=(F.transpose()*F-Matrix3d::Identity())*0.5f;
  const Matrix3d  derivE=(dF.transpose()*F+F.transpose()*dF)*0.5f;
  deriv=2.0f*G*(dF*E+F*derivE)+lambda*(E.trace()*dF+derivE.trace()*F);
}

void ElasticForceTetFullStVK::ddPF(Matrix3d& deriv,const Matrix3d& dF,const Matrix3d& dF2,const Matrix3d& F,const double&G,const double&lambda,const int&i) const{
  const Matrix3d derivE=(dF.transpose()*F+F.transpose()*dF)*0.5f;
  const Matrix3d derivE2=(dF2.transpose()*F+F.transpose()*dF2)*0.5f;
  const Matrix3d derivE12=(dF.transpose()*dF2+dF2.transpose()*dF)*0.5f;
  deriv=2.0f*G*(dF*derivE2+dF2*derivE+F*derivE12)+
	lambda*(derivE2.trace()*dF+derivE12.trace()*F+derivE.trace()*dF2);
}

double ElasticForceTetFullStVK::energy(const VectorXd &X){
  
  VectorM3x3 F;
  _def_grad.evalF(F,X);
  const ElasticMaterial<double>& mtl = _vol_mesh->material();

  double energy_phi = 0.0f;
  for (int i = 0; i < F.size(); ++i){

	const double lambda = mtl._lambda[i];
	const double mu = mtl._G[i];
    Matrix3d E = F[i].transpose()*F[i];
	E(0,0) -= 1.0f;
	E(1,1) -= 1.0f;
	E(2,2) -= 1.0f;
	E *= 0.5f;
	const double traceE = E(0,0)+E(1,1)+E(2,2);
	energy_phi += (lambda*traceE*traceE*0.5f+mu*contract(E,E))*_volume[i];
  }
  return energy_phi;
}
