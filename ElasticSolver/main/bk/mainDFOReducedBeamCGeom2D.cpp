#include "mainConfig.h"
#ifdef MAIN_COMPLEX_GEOMETRY_DFOREDUCED_BEAM_2D

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include "MakeMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamDFOReducedCGeom2D";
	boost::filesystem::create_directory("./"+name);

	boost::shared_ptr<FEMGeom> geom(new FEMGeom(2));
	ObjMesh m;
	MakeMesh::makeCapsule2D(m,0.5f,2.0f,16);
	Mat4 I=Mat4::Identity();
	for(scalar i=-20.0f;i<30.0f;i+=2.0f)
	{
		I.block<3,1>(0,3)=Vec3(i,-8.0f,0.0f);
		geom->addGeomMesh(I,m,0.1f);
	}
	geom->assemble();
	geom->writeVTK("./geom.vtk");

	DFOReducedFEMSolver sol(2);
	sol.resetParam(0.0f,0.0f,1E8f,1000.0f,1000);
	buildBeam(sol);
	buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.buildU(30);
	sol._geom=geom;
	sol.getMesh().calcMatDist();
	FEMMesh mesh2=sol.getMesh();

	//basis
	boost::filesystem::create_directory("./"+name+"/frm");
	FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
	profiler.reset(1000,&sol);
	for(sizeType i=0;i<3000;i++)
	{
        INFOV("%lu",i)
		profiler.beginFrame();
		sol.advance(1E-2f);
		profiler.endFrame();
		if(i%100 == 0 && i > 0)
		{
			sol.getMesh()+=mesh2;
			sol.clearEnergy();
			buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f);
		}
		if(i%155 == 0 && i > 0)
		{
			sol.getMesh()-=0;
			sol.clearEnergy();
			buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f);
		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
