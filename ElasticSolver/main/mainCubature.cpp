#include "mainConfig.h"
#ifdef MAIN_CUBATURE

#include "BeamModel.h"
#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMProfiler.h"
#include "FEMCubatureSolver.h"
#include "FEMCollision.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

void mainCubatureBeam(std::string name="CubatureBeam")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::STVK,NULL,2500000.0f,true,REDUCED_SYSTEM,true);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getMesh().getB(0)._cubature.reset(new Cubature);
    //sol.getMesh().getB(0)._cubature->read(boost::filesystem::ifstream("./meshCubature/mesh.CUBATURE",ios::binary));

    //sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.05f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(5,false);
    sol.getSystem<FEMReducedSystem>(0).buildDerivU(-1);
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getSystem<FEMReducedSystem>(0).onDirty();

    sol._geom.reset(new FEMGeom(2));
    sol._geom->addGeomBox(Mat4::Identity(),BBox<scalar>(Vec3(0.0f,-10.0f,0.0f),Vec3(5.0f,-0.2f,0.0f)));
    sol._geom->assemble();
    sol._geom->writeVTK("./"+name+"/geom.vtk");

    boost::filesystem::create_directory("./"+name+"/frm");
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        sol.advance(1E-2f);
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainCubatureDino(std::string name="CubatureDino")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(3);
    sol.getMesh().reset("./dino/mesh.ABQ",0.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::STVK,NULL,2500000.0f,true,REDUCED_SYSTEM,true);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.15f,-0.45f,-0.25f),Vec3(0.35f,-0.35f,0.25f)));
    sol.getSystem<FEMReducedSystem>(0).addEnergyMaterial("./dino/mesh.elastic",MaterialEnergy::STVK,true);
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.05f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(5,false);
    sol.getSystem<FEMReducedSystem>(0).buildDerivU(-1);
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getSystem<FEMReducedSystem>(0).onDirty();

    boost::filesystem::create_directory("./"+name+"/frm");
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        sol.advance(1E-2f);
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
int main()
{
    mainCubatureBeam();
    //mainCubatureDino();
}

#endif
