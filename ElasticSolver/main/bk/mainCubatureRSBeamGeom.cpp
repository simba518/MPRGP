#include "mainConfig.h"
#ifdef MAIN_GEOMETRY_CUBATURE_RS_BEAM

#include "../RSCoupledFEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCubatureRS10Geom";
	boost::filesystem::create_directory("./"+name);

	CubatureRotateReducedFEMSolver sol(3);
	buildFreeBeam(sol,-50.0f,MaterialEnergy::LINEAR,NULL,2500000.0f,32.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.buildU(30);
	sol.setUseRigid(true);

	sol._geom.reset(new FEMGeom(3));
	sol._geom->addGeomPlane(Mat4::Identity(),Plane(Vec3(0.0f,-5.0f,0.0f),Vec3(0.0f,1.0f,0.0f)));
	sol._geom->addGeomBox(OBB3D(Mat3::Identity(),Vec3(0.5f,-4.0f,0.0f),Vec3(1.0f,2.0f,1.0f)));
	sol._geom->addGeomSphere(Vec3(3.0f,-5.0f,1.0f),1.5f);
	sol._geom->assemble();
	sol._geom->writeVTK("./geom.vtk");

	//cubature
#define COMPUTE_CUBATURE
#ifdef COMPUTE_CUBATURE
#define REGENERATE_TRAINING_DATA
#ifdef REGENERATE_TRAINING_DATA
	sol.computeTrainingData(boost::filesystem::ofstream("./"+name+"/training.dat",ios::binary),1000,50.0f,true);
#endif
    boost::filesystem::ifstream is("./"+name+"/training.dat",ios::binary);
    sol.computeCubature(is,0.01f,10);
	sol.writeCubatureVTK("./"+name+"/cubature/");
	sol.writeCubature(boost::filesystem::ofstream("./"+name+"/cubature.dat",ios::binary));
#else 
	sol.readCubature(boost::filesystem::ifstream("./"+name+"/cubature.dat",ios::binary));
#endif 

	//basis
	boost::filesystem::create_directory("./"+name+"/frm");
	FEMProfiler profiler("./"+name+"/profile.txt",1E-2f);
	profiler.reset(1000,&sol);
	for(sizeType i=0;i<1000;i++)
	{
		INFOV("%d",i)
		profiler.beginFrame();
		sol.advance(1E-2f);
		profiler.endFrame();
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
