#ifndef BEAM_MODEL_H
#define BEAM_MODEL_H

#include "../FEMSolver.h"

PRJ_BEGIN

class FuncSolid : public ImplicitFunc<scalar>
{
public:
	FuncSolid(sizeType dim):_dim(dim){}
	scalar operator()(const Vec3& pos) const
	{
		if(_dim == 3)
		return std::max(
			   std::max(std::max(0.1f-pos.x(),pos.x()-3.0f),
						std::max(0.1f-pos.y(),pos.y()-0.9f)),
						std::max(0.1f-pos.z(),pos.z()-0.9f));
		else if(_dim == 2)
		return std::max(std::max(0.1f-pos.x(),pos.x()-3.0f),
						std::max(0.1f-pos.y(),pos.y()-0.9f));
		else return 0.0f;
	}
	sizeType _dim;
};
void buildBeamEnergy(FEMSolver& sol,scalar g,MaterialEnergy::TYPE type=MaterialEnergy::LINEAR,BBox<scalar>* bbm=NULL,scalar stiffness=2500000.0f,bool invertible=true)
{
	sol.clearEnergy();
	sol.addEnergyMaterial(stiffness,0.49f,type,invertible);
	sol.addEnergyMass(Vec3(0.0f,g,0.0f),bbm);
}
void buildBeam(FEMSolver& sol,scalar res=32.0f)
{
	scalar sz=1.0f/res;
	BBox<scalar> bb(Vec3(0.0f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
	sol.getMesh().reset(bb,FuncSolid(sol.dim()),sz,sz*5.0f);
}

PRJ_END

#endif
