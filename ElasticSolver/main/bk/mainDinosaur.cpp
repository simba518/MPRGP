#include "mainConfig.h"
#ifdef MAIN_DINOSAUR

#include "../FEMSolver.h"
#include "../FEMProfiler.h"
#include "DinoModel.h"
#include "ParticleCD.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="dinosaurCorotational";
	boost::filesystem::create_directory("./"+name);

	DFOReducedFEMSolver sol(3);
	sol.resetParam(0.0f,0.0f,1E5f,1000.0f,1);
    buildFreeDino(sol);
	buildDinoEnergy(sol,-9.81f,MaterialEnergy::LINEAR);
	sol.buildU(30);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	FEMMesh mesh2=sol.getMesh();

	sol._geom.reset(new FEMGeom(3));
	sol._geom->addGeomPlane(Mat4::Identity(),Plane(Vec3(0.0f,-1.0f,0.0f),Vec3(0.0f,1.0f,0.0f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3(-0.35f,-1.0f, 0.35f),Vec3(0.1f,4.0f,0.1f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3(-0.35f,-1.0f,-0.35f),Vec3(0.1f,4.0f,0.1f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3( 0.35f,-1.0f, 0.35f),Vec3(0.1f,4.0f,0.1f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3( 0.35f,-1.0f,-0.35f),Vec3(0.1f,4.0f,0.1f)));
	sol._geom->addGeomSphere(Vec3(0.0f,-1.0f,0.0f),0.5f);
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
		if(i%100 == 0 && i > 0)
		{
			sol.getMesh()+=mesh2;
			sol.clearEnergy();
			buildDinoEnergy(sol,-9.81f,MaterialEnergy::LINEAR);
		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
