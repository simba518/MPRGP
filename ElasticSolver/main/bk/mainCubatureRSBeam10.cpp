#include "mainConfig.h"
#ifdef MAIN_CUBATURE_RS_BEAM_10

#include "../RSCoupledFEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include "MakeMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCubatureRS10";
	boost::filesystem::create_directory("./"+name);

	CubatureRotateReducedFEMSolver sol(2);
	sol.resetParam(0.01f,0.0f,1E7f,1000.0f,10);
	buildBeam(sol,32.0f);
//#define READ
#ifdef READ
	sol.getMesh().read(boost::filesystem::ifstream("./"+name+"/mesh.dat",ios::binary));
	sol.getMesh().clearVel();
	sol.getMesh().resetPos();
	sol.getMesh().calcMatDist();
#endif
	buildBeamEnergy(sol,-9.81f,MaterialEnergy::LINEAR,NULL,1000000.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.setUseRigid(true);
	
	//cubature
#ifndef READ
	sol.buildU(30);
#define COMPUTE_TRAINING_DATA
#ifdef COMPUTE_TRAINING_DATA
    {boost::filesystem::ofstream os("./"+name+"/training.dat",ios::binary);
    sol.computeTrainingData(os,1000,100.0f,true);}
#endif
    boost::filesystem::ifstream is("./"+name+"/training.dat",ios::binary);
    sol.computeCubature(is,0.01f,10);
	Mat4 T=Mat4::Identity();
	T.block<3,1>(0,3)[1]=7.0f;
	sol.getMesh().applyTrans(T,0,true,true);
	sol.getMesh().clearVel();
	sol.getMesh().resetPos();
	sol.getMesh().calcMatDist();
	sol.getMesh().updateMesh();
    {boost::filesystem::ofstream os("./"+name+"/mesh.dat",ios::binary);
    sol.getMesh().write(os);}
	sol.writeCubatureVTK("./"+name+"/cubature/");
#endif
	FEMMesh mesh2=sol.getMesh();
	/*{
		sol.getMesh()+=mesh2;
		sol.getMesh()+=mesh2;
		sol.clearEnergy();
		buildBeamEnergy(sol,-9.81f,MaterialEnergy::LINEAR,NULL,1000000.0f);
	}
	FEMFrameData dat;
	dat.read(boost::filesystem::ifstream("./"+name+"/m200.dat",ios::binary));
	sol.setFrame(dat);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");*/

	ObjMesh m;
	boost::shared_ptr<FEMGeom> geom(new FEMGeom(2));
	MakeMesh::makeCapsule2D(m,0.5f,2.0f,16);
	for(scalar i=-20.0f;i<30.0f;i+=2.0f)
	{
		Mat4 I=Mat4::Identity();
		I.block<3,1>(0,3)=Vec3(i,-5.0f,0.0f);
		geom->addGeomMesh(I,m,0.1f);
	}
	geom->assemble();
	geom->writeVTK("./geom.vtk");
	sol._geom=geom;

	//basis
	boost::filesystem::create_directory("./"+name+"/frm");
	FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
	profiler.reset(1000,&sol);
	for(sizeType i=0;i<1000;i++)
	{
        INFOV("%lu",i)
		profiler.beginFrame();
		sol.advance(1E-2f);
		profiler.endFrame();
		if(i%100 == 0 && i > 0)
		{
			sol.getMesh()+=mesh2;
			sol.clearEnergy();
			buildBeamEnergy(sol,-9.81f,MaterialEnergy::LINEAR,NULL,1000000.0f);
		}
		if(i == 200)
		{
			FEMFrameData dat;
			sol.getFrame(dat);
            boost::filesystem::ofstream os("./"+name+"/m200.dat",ios::binary);
            dat.write(os);
		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
