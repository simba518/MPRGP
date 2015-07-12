#ifndef DINO_MODEL_H
#define DINO_MODEL_H

#include "../FEMSolver.h"

PRJ_BEGIN

void buildDino(FEMSolver& sol,scalar res=32.0f,scalar rho=1000.0f)
{
	scalar sz=1.0f/res;
    boost::filesystem::ifstream is("./dino/mesh.abq");
    sol.getMesh().reset(is,sz);
	sol.getMesh().calcMatDist();
	Vec3 pt0(0.0737585283264465f,-0.485245085442333f,0.0391006922458308f);
	Vec3 len0(0.2f,0.1f,0.2f);
	sol.fixPointIn(BBox<scalar>(pt0-len0,pt0+len0));
}
void buildFreeDino(FEMSolver& sol,scalar res=32.0f,scalar rho=1000.0f)
{
	scalar sz=1.0f/res;
    boost::filesystem::ifstream is("./dino/mesh.abq");
    sol.getMesh().reset(is,sz);
	Mat4 T=Mat4::Identity();
	T.block<3,1>(0,3)=Vec3(0.0f,2.5f,0.0f);
	sol.getMesh().applyTrans(T,0,true,true);
	sol.getMesh().calcMatDist();
}
void buildDinoEnergy(FEMSolver& sol,scalar g,MaterialEnergy::TYPE type=MaterialEnergy::LINEAR,BBox<scalar>* bbm=NULL,scalar stiffness=2500000.0f)
{
	sol.addEnergyMass(Vec3(0.0f,g,0.0f),bbm);
	sol.addEnergyMaterial(stiffness,0.49f,type);
}

PRJ_END

#endif
