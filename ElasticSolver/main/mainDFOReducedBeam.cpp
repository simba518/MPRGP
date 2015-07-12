#include "mainConfig.h"
#ifdef MAIN_DFO_REDUCED_BEAM

#include "FEMGeom.h"
#include "FEMProfiler.h"
#include "BeamModel.h"
#include "MakeMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

void mainDFOReducedBeam1(std::string name="DFOReducedBeam1")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    //sol.getMesh().reset("./mesh.ABQ",0.0f);
    buildBeam(sol,16.0f);
    BBox<scalar> bb(Vec3(-0.2f,-10.0f,0.0f),Vec3(0.2f,10.0f,0.0f));
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,&bb,250000.0f,true,DFO_REDUCED_SYSTEM);
    sol.resetImplicitEuler(1E-4f,10);
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.0f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(30,false);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    //sol.getSystem<FEMReducedSystem>(0).debugU();
    //sol.readFrame(boost::filesystem::ifstream("./frm.dat",ios::binary));

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
    profiler.reset(500,&sol);
    for(sizeType i=0; i<500; i++) {
        INFOV("%lu",i)
        profiler.beginFrame();
        sol.advance(1E-2f);
        // if(i == 100)
        //     sol.writeFrame(boost::filesystem::ofstream("./frm.dat",ios::binary));
        profiler.endFrame();
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainDFOReducedBeam2(std::string name="DFOReducedBeam2")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,32.0f);
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,true,DFO_REDUCED_SYSTEM);
    sol.resetImplicitEuler(1E-4f,100);
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.01f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(2,false);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.setSelfColl(true);
    sol.setCollK(1E7f);

    ObjMesh m;
    boost::shared_ptr<FEMGeom> geom(new FEMGeom(2));
    MakeMesh::makeCapsule2D(m,0.5f,2.0f,16);
    Mat4 I=Mat4::Identity();
    for(scalar i=-20.0f; i<30.0f; i+=2.0f) {
        I.block<3,1>(0,3)=Vec3(i,-8.0f,0.0f);
        geom->addGeomMesh(I,m,0.1f);
    }
    geom->assemble();
    geom->writeVTK("./"+name+"/geom.vtk");
    sol._geom=geom;

    FEMMesh mesh2=sol.getMesh();
    Mat4 R=Mat4::Identity();
    R.block<3,1>(0,3)=Vec3::Unit(1)*-3.0f;
    mesh2.applyTrans(R,-1,true,true);
    sol.getMesh()+=mesh2;
    //sol.getSystem<FEMReducedSystem>(1).clearEnergy(typeid(FixedForceEnergy).name());
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");

    //basis
    boost::filesystem::create_directory("./"+name+"/frm");
    for(sizeType i=0; i<3000; i++) {
        INFOV("%lu",i)
        sol.advance(1E-2f);
        if(i%100 == 0 && i > 0)
            sol.getMesh()+=mesh2;
        if(i%155 == 0 && i > 0)
            sol.getMesh()-=0;
        ostringstream ossm;
        ossm << "./"+name+"/frm/m" << i << ".vtk";
        sol.getMesh().writeVTK(ossm.str());
    }
}
void mainDFOReducedBeam3(std::string name="DFOReducedBeam3")
{
    boost::filesystem::create_directory("./"+name);
    FEMSolver sol(2);
    buildBeam(sol,16.0f);
    BBox<scalar> bb(Vec3(-0.2f,-2.0f,0.0f),Vec3(0.2f,2.0f,0.0f));
    buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,&bb,250000.0f,true,RIGID_RS_REDUCED_SYSTEM);
    sol.resetImplicitEuler(1E-4f,10);
    sol.getSystem<FEMReducedSystem>(0).resetParam(0.0f,0.0f,1000.0f);
    sol.getSystem<FEMReducedSystem>(0).buildU(10,false);
    sol.getSystem<FEMReducedSystem>(0).writeBasisVTK("./"+name,1.0f);
    sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
    sol.getMesh().writePSetVTK("./"+name+"/pSet.vtk");
    sol.setSelfColl(false);
    //sol.getSystem<FEMReducedSystem>(0).debugU();

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

int main()
{
    //mainDFOReducedBeam1();
    //mainDFOReducedBeam2();
    mainDFOReducedBeam3();
}

#endif
