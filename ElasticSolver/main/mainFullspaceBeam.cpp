#include "mainConfig.h"
#ifdef MAIN_FULLSPACE_BEAM

#include "FEMSystem.h"
#include "FEMGeom.h"
#include "FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

void mainCorotationalBeam1(scalar beta,scalar beta2,scalar gamma,scalar alphaD,scalar betaD,std::string name="CorotationalBeam1")
{
    boost::filesystem::create_directory("./"+name);

    FEMSolver sol(2);
    buildBeam(sol,16.0f);
    buildBeamEnergy(sol,-300.0f,MaterialEnergy::LINEAR,NULL,2500000.0f);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    sol.getSystem<FEMSystem>(0).resetParam(alphaD,betaD,1000.0f);
    sol.resetParam(beta,beta2,gamma,1E-4f,1);
    //sol.getSystem<FEMSystem>(0).debugU();

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(300,&sol);
    for(sizeType i=0; i<300; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainCorotationalBeam2(MaterialEnergy::TYPE type,std::string name="CorotationalBeam2")
{
    boost::filesystem::create_directory("./"+name);

    FEMSolver sol(2);
    buildBeam(sol,16.0f);
    buildBeamEnergy(sol,-300.0f,type,NULL,2500000.0f);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    FEMMesh mesh2=sol.getMesh();
    Mat4 trans=Mat4::Identity();
    trans.block<3,1>(0,3)=Vec3::Unit(0)*6.0f;
    mesh2.applyTrans(trans,0,true,true);
    sol.getMesh()+=mesh2;
    sol.getSystem<FEMSystem>(1).clearConstraint();
    sol.getSystem<FEMSystem>(1).addConstraintPoint(FEMInterp(5,Vec4::Constant(1.0f/3.0f)),1E6f);
    sol.setSelfColl(false);

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(200,&sol);
    for(sizeType i=0; i<200; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainCorotationalBeam3(std::string name="CorotationalBeam3")
{
    boost::filesystem::create_directory("./"+name);

    FEMSolver sol(2);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::COROTATIONAL,NULL,2500000.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    FEMMesh mesh2=sol.getMesh();
    Mat4 trans=Mat4::Identity();
    trans.block<3,1>(0,3)=Vec3(1.0f,3.0f,0.0f);
    mesh2.applyTrans(trans,0,true,true);
    sol.getMesh()+=mesh2;
    sol.setCollK(1E6f);

    sol._geom.reset(new FEMGeom(2));
    sol._geom->addGeomPlane(Mat4::Identity(),Vec4(0.0f,1.0f,0.0f,4.0f));
    sol._geom->assemble();
    sol._geom->writeVTK("./"+name+"/geom.vtk");

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(200,&sol);
    for(sizeType i=0; i<200; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainCorotationalSTVKBeam4(std::string name="CorotationalBeam4")
{
    boost::filesystem::create_directory("./"+name);

    FEMSolver sol(3);
    buildBeam(sol,16.0f);
    buildBeamEnergy(sol,-300.0f,MaterialEnergy::COROTATIONAL,NULL,2500000.0f,false);
    sol.addConstraintPoint(BBox<scalar>(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f)));
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    FEMMesh mesh2=sol.getMesh();
    Mat4 trans=Mat4::Identity();
    trans.block<3,1>(0,3)=Vec3::Unit(0)*6.0f;
    mesh2.applyTrans(trans,0,true,true);
    sol.getMesh()+=mesh2;
    sol.getSystem<FEMSystem>(1).clearConstraint();
    sol.getSystem<FEMSystem>(1).addConstraintPoint(FEMInterp(5,Vec4::Constant(1.0f/3.0f)),1E6f);
    sol.setSelfColl(false);
    sol.getSystem<FEMSystem>(0).setMaterialType(MaterialEnergy::STVK);
    //sol.getSystem<FEMSystem>(0).energys().debugEnergy();

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        if(i == 100)
            sol.writeFrame(boost::filesystem::ofstream("./frm.dat",ios::binary));
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void main()
{
    //mainCorotationalBeam1(1.0f,0.5f,1.0f,0.0f,0.0f,"CorotationalBeamImplicit");
    //mainCorotationalBeam1(0.25f,0.25f,0.5f,0.0f,0.0f,"CorotationalBeamNewmark");
    //mainCorotationalBeam1(0.25f,0.25f,0.5f,0.1f,0.0f,"CorotationalBeamNewmarkAlpha");
    //mainCorotationalBeam1(0.25f,0.25f,0.5f,0.0f,0.1f,"CorotationalBeamNewmarkBeta");
    //mainCorotationalBeam2(MaterialEnergy::COROTATIONAL,"CorotationalBeam2");
    //mainCorotationalBeam2(MaterialEnergy::COROTATIONAL_EXACT,"CorotationalExactBeam2");
    mainCorotationalBeam3();
    //mainCorotationalSTVKBeam4();
}

#endif
