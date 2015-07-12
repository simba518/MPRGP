#include "mainConfig.h"
#ifdef MAIN_LMA_BASIS

#include "FEMMesh.h"
#include "FEMReducedSystem.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main(int argc,char** argv)
{
    if(argc < 6) {
        WARNING("Usage: exe [path] [NR] [off] [young] [poisson]")
        exit(EXIT_FAILURE);
    }

    //ABQ file
    istringstream issPath(argv[1]);

    //Nr basis
    sizeType nr;
    istringstream issNr(argv[2]);
    issNr >> nr;
    bool recompute=nr < 0;
    nr=abs(nr);

    //writeVTK offset
    scalar off;
    istringstream issOff(argv[3]);
    issOff >> off;

    //young
    scalar young;
    istringstream issYoung(argv[4]);
    issYoung >> young;

    //poisson
    scalar poisson;
    istringstream issPoisson(argv[5]);
    issPoisson >> poisson;

    //solve
    FEMSolver sol(3);
    //buildBeam(sol,32.0f);
    sol.getMesh().reset(issPath.str(),0.0f);
    sol.getMesh().getB(0)._system.reset(new FEMReducedSystem(sol.getMesh().getB(0)));
    sol.getSystem<FEMReducedSystem>(0).addEnergyMaterial(young,poisson,MaterialEnergy::STVK,false);
    sol.getSystem<FEMReducedSystem>(0).buildU(nr,recompute);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK(boost::filesystem::path(issPath.str()).parent_path().string(),off);
    sol.getMesh().getB(0).writeABQ(issPath.str());
    return 0;
}

#endif