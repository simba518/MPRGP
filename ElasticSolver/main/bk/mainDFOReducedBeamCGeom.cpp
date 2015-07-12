#include "mainConfig.h"
#ifdef MAIN_COMPLEX_GEOMETRY_DFOREDUCED_BEAM

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include "MakeMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamDFOReducedCGeom";
	boost::filesystem::create_directory("./"+name);

	boost::shared_ptr<FEMGeom> geom(new FEMGeom(3));
	geom->addGeomPlane(Mat4::Identity(),Plane(Vec3(0.0f,-5.0f,0.0f),Vec3(0.0f,1.0f,0.0f)));
	Mat4 I=Mat4::Identity();
	I.block<3,1>(0,3)=Vec3(0.0f,-4.0f,0.0f);
	geom->addGeomMesh(I,"./bunny.obj",0.1f);
	geom->assemble();
	geom->writeVTK("./geom.vtk");

	DFOReducedFEMSolver sol(3);
	sol.resetParam(0.0f,0.0f,1E8f,1000.0f,1000);
	buildBeam(sol);
	buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.buildU(30);
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
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
	return 0;
}
#endif
