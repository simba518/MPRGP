#include "mainConfig.h"
#ifdef MAIN_DFOREDUCED_BEAM

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamDFOReduced";
	boost::filesystem::create_directory("./"+name);

	DFOReducedFEMSolver sol(2);
	sol.resetParam(0.0f,0.0f,1E8f,1000.0f,1);
	buildBeam(sol);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.getMesh().calcMatDist();
	BBox<scalar> bb(Vec3(2.5f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
	buildBeamEnergy(sol,-9.81f,MaterialEnergy::LINEAR,&bb,25000.0f);
	sol.buildU(30);
	FEMMesh mesh2=sol.getMesh();
//	{
//		Mat4 T=Mat4::Identity();
//		T.block<3,1>(0,3)=-Vec3::Unit(1)*5.0f;
//		mesh2.applyTrans(T,0,true,true);
//		sol.getMesh()+=mesh2;

//		FEMFrameData dat;
//        boost::filesystem::ifstream is("./"+name+"/frm230.dat",ios::binary);
//        dat.read(is);
//		sol.setFrame(dat);
//		sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
//	}

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
			T.block<3,1>(0,3)=-Vec3::Unit(1)*5.0f;
			mesh2.applyTrans(T,0,true,true);
			sol.getMesh()+=mesh2;

			sol.clearEnergy();
			BBox<scalar> bb(Vec3(2.5f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
			buildBeamEnergy(sol,-9.81f,MaterialEnergy::LINEAR,&bb,25000.0f);
		}
//		if(i == 230)
//		{
//			FEMFrameData dat;
//			sol.getFrame(dat);
//			dat.write(boost::filesystem::ofstream("./"+name+"/frm230.dat",ios::binary));
//		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
