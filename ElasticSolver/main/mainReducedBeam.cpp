#include "mainConfig.h"
#ifdef MAIN_REDUCED_BEAM

#include "FEMSystem.h"
#include "FEMGeom.h"
#include "FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

void mainReducedBeam1(scalar beta,scalar beta2,scalar gamma,scalar alphaD,scalar betaD,std::string name="ReducedBeam1")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,16.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,true,REDUCED_SYSTEM);
    sol.resetParam(beta,beta2,gamma);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(alphaD,betaD,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(30);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    sol.getSystem<FEMReducedSystem>(0).debugU();

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainReducedBeam2(std::string name="ReducedBeam2")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    //sol.getMesh().reset("./mesh.ABQ",0.0f);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,true,REDUCED_SYSTEM);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.05f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(2,false);
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    sol.setCollK(1E6f);

    sol._geom.reset(new FEMGeom(2));
    sol._geom->addGeomBox(Mat4::Identity(),BBox<scalar>(Vec3(0.0f,-10.0f,0.0f),Vec3(5.0f,-1.0f,0.0f)));
    sol._geom->assemble();
    sol._geom->writeVTK("./"+name+"/geom.vtk");

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        //if(i == 0)
        //	sol.getMesh().getB(0).writeABQ("./mesh.ABQ");
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainReducedBeam3(std::string name="ReducedBeam3")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,16.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,true,REDUCED_SYSTEM);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).buildU(2);
    sol.getSystem<FEMSystem>(0).resetParam(0.01f,0.0f,1000.0f);

    FEMMesh mesh2=sol.getMesh();
    Mat4 trans=Mat4::Identity();
    trans.block<3,1>(0,3)=Vec3(1.0f,4.0f,0.0f);
    mesh2.applyTrans(trans,0,true,true);
    sol.getMesh()+=mesh2;
    sol.getMesh().getB(1)._system.reset(new FEMSystem(sol.getMesh().getB(1)));
    sol.getSystem<FEMSystem>(1).addEnergyMass(Vec3(0.0f,-50.0f,0.0f));
    sol.getSystem<FEMSystem>(1).addEnergyMaterial(2500000.0f,0.49f,MaterialEnergy::COROTATIONAL_EXACT);

    sol.setSelfColl(true);
    sol.setCollK(1E6f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    sol._geom.reset(new FEMGeom(2));
    sol._geom->addGeomBox(Mat4::Identity(),BBox<scalar>(Vec3(-50.0f,-10.0f,0.0f),Vec3(50.0f,-1.0f,0.0f)));
    sol._geom->assemble();
    sol._geom->writeVTK("./"+name+"/geom.vtk");
    //sol.readFrame(boost::filesystem::ifstream("./frm.dat",ios::binary));

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        //if(i == 100)
        //	sol.writeFrame(boost::filesystem::ofstream("./frm.dat",ios::binary));
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainReducedBeam4(std::string name="ReducedBeam4")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    //sol.getMesh().reset("./mesh.ABQ",0.0f);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::STVK,NULL,2500000.0f,false,REDUCED_SYSTEM);
    sol.resetImplicitEuler(0.1f,3);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.01f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(5,false);
    sol.getSystem<FEMReducedSystem>(0).buildDerivU();
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getSystem<FEMReducedSystem>(0).setMaterialType(MaterialEnergy::STVK);

    FEMMesh mesh2=sol.getMesh();
    Mat4 trans=Mat4::Identity();
    trans.block<3,1>(0,3)=Vec3(1.0f,4.0f,0.0f);
    mesh2.applyTrans(trans,0,true,true);
    sol.getMesh()+=mesh2;
    sol.getMesh().getB(1)._system.reset(new FEMSystem(sol.getMesh().getB(1)));
    sol.getSystem<FEMSystem>(1).addEnergyMass(Vec3(0.0f,-50.0f,0.0f));
    sol.getSystem<FEMSystem>(1).addEnergyMaterial(2500000.0f,0.49f,MaterialEnergy::COROTATIONAL_EXACT);

    sol.setSelfColl(true);
    sol.setCollK(1E6f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    sol._geom.reset(new FEMGeom(2));
    sol._geom->addGeomBox(Mat4::Identity(),BBox<scalar>(Vec3(-50.0f,-10.0f,0.0f),Vec3(50.0f,-1.0f,0.0f)));
    sol._geom->assemble();
    sol._geom->writeVTK("./"+name+"/geom.vtk");

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        //sol.getMesh().getB(0).writeABQ("./mesh.ABQ");
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainReducedBeam5(std::string name="ReducedBeam5")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-300.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,true,WARPPED_REDUCED_SYSTEM);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.05f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(5);
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    sol.setCollK(1E6f);

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainReducedBeam6(bool coupled,std::string name="ReducedBeam6")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-300.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,false,coupled?RS_REDUCED_SYSTEM:POST_PROCESS_RS_REDUCED_SYSTEM);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.1f,2.0f,0.0f)));
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.05f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(10);
    sol.getSystem<FEMReducedSystem>(0).setLocalBasis(0.6f,0.45f,true);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    sol.setCollK(1E6f);

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        {
            scalarD EK,EP;
            sol.getSystemEnergy(EK,EP);
            INFOV("EK: %f EP: %f",EK,EP)
        }
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void main()
{
    //mainReducedBeam1(1.0f,0.5f,1.0f,0.0f,0.0f,"ReducedBeamImplicit");
    //mainReducedBeam1(0.25f,0.25f,0.5f,0.0f,0.0f,"ReducedBeamNewmark");
    //mainReducedBeam1(0.25f,0.25f,0.5f,0.1f,0.0f,"ReducedBeamNewmarkAlpha");
    //mainReducedBeam1(0.25f,0.25f,0.5f,0.0f,0.1f,"ReducedBeamNewmarkBeta");
    //mainReducedBeam2();
    //mainReducedBeam3();
    //mainReducedBeam4();
    //mainReducedBeam5();
    //mainReducedBeam6(false);
    //mainReducedBeam6(true);
}

#endif
