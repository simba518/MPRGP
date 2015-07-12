#include "mainConfig.h"
#ifdef MAIN_COROTATIONAL_BEAM

#include "../RSCoupledFEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCorotational";
	boost::filesystem::create_directory("./"+name);

	VariationalFEMSolver sol(2);
	sol.resetParam(0.01f,0.0f,1E5f,1000.0f,10);
	buildBeam(sol,16.0f);
	buildBeamEnergy(sol,-300.0f,MaterialEnergy::COROTATIONAL,NULL,2500000.0f);
	sol.fixPointIn(BBox<scalar>(Vec3(0.0f,-2.0f,-2.0f),Vec3(0.2f,2.0f,2.0f)));
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	FEMMesh mesh2=sol.getMesh();

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
		/*if(i == 100)
		{
			Mat4 T=Mat4::Identity();
			T.block<3,1>(0,3)=Vec3::Unit(0)*5.0f;
			mesh2.applyTrans(T,0,true,true);
			sol.getMesh()+=mesh2;

			sol.clearEnergy();
			sol.fixPointIn(BBox<scalar>(Vec3(5.0f,-2.0f,-2.0f),Vec3(5.2f,2.0f,2.0f)));
			buildBeamEnergy(sol,-300.0f,MaterialEnergy::COROTATIONAL_EXACT,NULL,2500000.0f);
		}*/
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
