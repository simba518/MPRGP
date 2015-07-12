#include "mainConfig.h"
#ifdef MAIN_WARPPED_REDUCED_BEAM

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamWarppedReduced";
	boost::filesystem::create_directory("./"+name);

	WarppedReducedFEMSolver sol(2);
	sol.resetParam(1.0f,0.0f,1E8f,1000.0f,1);
//#define READ
#ifndef READ
	buildBeam(sol);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.fixPointIn(BBox<scalar>(Vec3(0.0f,-2.0f,-2.0f),Vec3(0.2f,2.0f,2.0f)));
	buildBeamEnergy(sol,-50.0f);
	sol.buildU(30);
	sol.writeBasisVTK(1.0f);
    boost::filesystem::ofstream os("./mesh.dat",ios::binary);
    sol.getMesh().write(os);
#else
	sol.getMesh().read(boost::filesystem::ifstream("./mesh.dat",ios::binary));
	sol.fixPointIn(BBox<scalar>(Vec3(0.0f,-2.0f,-2.0f),Vec3(0.2f,2.0f,2.0f)));
	buildBeamEnergy(sol,-50.0f);
#endif
	FEMMesh mesh2=sol.getMesh();
	{
		FEMFrameData dat;
        boost::filesystem::ifstream is("./"+name+"/frm50.dat",ios::binary);
        dat.read(is);
		sol.setFrame(dat);
	}

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
		if(i == 100)
		{
			Mat4 T=Mat4::Identity();
			T.block<3,1>(0,3)=Vec3::Unit(0)*5.0f;
			mesh2.applyTrans(T,0,true,true);
			sol.getMesh()+=mesh2;

			sol.clearEnergy();
			sol.fixPointIn(BBox<scalar>(Vec3(0.0f,-2.0f,-2.0f),Vec3(0.2f,2.0f,2.0f)));
			sol.fixPointIn(BBox<scalar>(Vec3(5.0f,-2.0f,-2.0f),Vec3(5.2f,2.0f,2.0f)));
			buildBeamEnergy(sol,-50.0f);
		}
		if(i == 50)
		{
			FEMFrameData dat;
			sol.getFrame(dat);
            boost::filesystem::ofstream os("./"+name+"/frm50.dat",ios::binary);
            dat.write(os);
		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
