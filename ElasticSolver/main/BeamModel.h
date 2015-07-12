#ifndef BEAM_MODEL_H
#define BEAM_MODEL_H

#include "FEMSystem.h"
#include "FEMReducedSystem.h"
#include "FEMMesh.h"
#include "ImplicitFunc.h"

PRJ_BEGIN

enum SYSTEM_TYPE {
    SYSTEM,
    REDUCED_SYSTEM,
    DFO_REDUCED_SYSTEM,
    WARPPED_REDUCED_SYSTEM,
    //basic RS version
    POST_PROCESS_RS_REDUCED_SYSTEM,
    RS_REDUCED_SYSTEM,
    RIGID_RS_REDUCED_SYSTEM,
    //cubature RS version
    POST_PROCESS_CUBATURE_RS_REDUCED_SYSTEM,
    CUBATURE_RS_REDUCED_SYSTEM,
    RIGID_CUBATURE_RS_REDUCED_SYSTEM,
};
class FuncSolid : public ImplicitFunc<scalar>
{
public:
    FuncSolid(sizeType dim,scalar l):_dim(dim),_l(l) {}
    scalar operator()(const Vec3& pos) const {
        if(_dim == 3)
            return std::max(
                       std::max(std::max(0.1f-pos.x(),pos.x()-_l),
                                std::max(0.1f-pos.y(),pos.y()-0.9f)),
                       std::max(0.1f-pos.z(),pos.z()-0.9f));
        else if(_dim == 2)
            return std::max(std::max(0.1f-pos.x(),pos.x()-_l),
                            std::max(0.1f-pos.y(),pos.y()-0.9f));
        else return 0.0f;
    }
    sizeType _dim;
    scalar _l;
};
void buildBeamEnergy(FEMSolver& sol,scalar g,
                     MaterialEnergy::TYPE type=MaterialEnergy::LINEAR,BBox<scalar>* bbm=NULL,
                     scalar stiffness=2500000.0f,bool invertible=true,SYSTEM_TYPE stype=SYSTEM,bool cubature=false)
{
    switch(stype) {
    case SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMSystem(sol.getMesh().getB(0)));
        break;
    case REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::LMA,cubature);
        break;
    case DFO_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::RIGID_LMA,cubature);
        break;
    case WARPPED_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::WARPPED_LMA,cubature);
        break;

    case POST_PROCESS_RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::POST_PROCESS_RS,cubature);
        break;
    case RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::COUPLED_RS,cubature);
        break;
    case RIGID_RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::RIGID_COUPLED_RS,cubature);
        break;

    case POST_PROCESS_CUBATURE_RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::POST_PROCESS_CUBATURE_RS,cubature);
        break;
    case CUBATURE_RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::CUBATURE_COUPLED_RS,cubature);
        break;
    case RIGID_CUBATURE_RS_REDUCED_SYSTEM:
        sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
        sol.getSystem<FEMReducedSystem>(0).setReducedModel(FEMReducedSystem::RIGID_CUBATURE_COUPLED_RS,cubature);
        break;
    }
    FEMSystem& sys=*(sol.getMesh().getB(0)._system);
    sys.clearEnergy();
    sys.addEnergyMaterial(stiffness,0.49f,type,invertible);
    sys.addEnergyMass(Vec3(0.0f,g,0.0f),bbm);
}
void buildBeam(FEMSolver& sol,scalar res=32.0f,scalar l=3.0f)
{
    scalar sz=1.0f/res;
    BBox<scalar> bb(Vec3(0.0f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
    sol.getMesh().reset(bb,FuncSolid(sol.getMesh().dim(),l),sz,sz*5.0f);
}

PRJ_END

#endif
