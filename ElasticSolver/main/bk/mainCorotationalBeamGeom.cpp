#include "mainConfig.h"
#ifdef MAIN_GEOMETRY_COROTATIONAL_BEAM

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCorotationalGeom";
	boost::filesystem::create_directory("./"+name);

	FEMSolver sol(3);
	buildBeam(sol);
	buildBeamEnergy(sol,-50.0f,MaterialEnergy::COROTATIONAL,NULL,25000000.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	
	sol._geom.reset(new FEMGeom(3));
	sol._geom->addGeomPlane(Mat4::Identity(),Plane(Vec3(0.0f,-5.0f,0.0f),Vec3(0.0f,1.0f,0.0f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3(0.5f,-4.0f,0.0f),Vec3(1.0f,2.0f,1.0f)));
	sol._geom->addGeomSphere(Vec3(3.0f,-5.0f,1.0f),1.5f);
	sol._geom->assemble();
	sol._geom->writeVTK("./geom.vtk");

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
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
