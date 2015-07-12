#include "mainConfig.h"
#ifdef MAIN_RS_BEAM_LARGE_FORCE

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamRSLargeForce";
	boost::filesystem::create_directory("./"+name);

	RotateReducedFEMSolver sol(2);
	sol.resetParam(1.0f,0.0f,1E8f,1000.0f,1);
	buildBeam(sol);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	{
		FEMMesh mesh2=sol.getMesh();
		Mat4 T=Mat4::Identity();
		T.block<3,1>(0,3)=Vec3::Unit(0)*5.0f;
		mesh2.applyTrans(T,0,true,true);
		sol.getMesh()+=mesh2;
	}
	
	sol.fixPointIn(BBox<scalar>(Vec3(0.0f,-2.0f,-2.0f),Vec3(0.2f,2.0f,2.0f)));
	sol.fixPointIn(BBox<scalar>(Vec3(7.0f,-2.0f,-2.0f),Vec3(7.2f,2.0f,2.0f)));

	buildBeamEnergy(sol,-300.0f);
	sol.buildU(30,sol.getMesh().getB(0));
	sol.buildU(5,sol.getMesh().getB(1));
	sol.writeBasisVTK(1.0f);
    boost::filesystem::ofstream os("./mesh.dat",ios::binary);
    sol.getMesh().write(os);
	
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
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str(),NULL,&(sol.getX()));
	}
    return 0;
}
#endif
