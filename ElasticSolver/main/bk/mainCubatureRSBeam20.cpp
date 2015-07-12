#include "mainConfig.h"
#ifdef MAIN_CUBATURE_RS_BEAM_20

#include "../RSCoupledFEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCubatureRS20";
	boost::filesystem::create_directory("./"+name);

	CubatureRotateReducedFEMSolver sol(2);
	buildBeam(sol,-50.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.buildU(30);

	//cubature
#define COMPUTE_CUBATURE
#ifdef COMPUTE_CUBATURE
#define REGENERATE_TRAINING_DATA
#ifdef REGENERATE_TRAINING_DATA
   { boost::filesystem::ofstream os("./"+name+"/training.dat",ios::binary);
    sol.computeTrainingData(os,1000,5.0f,true);}
#endif
    boost::filesystem::ifstream is("./"+name+"/training.dat",ios::binary);
    sol.computeCubature(is,0.01f,20);
    sol.writeCubatureVTK("./"+name+"/cubature/");
//    sol.writeCubature(boost::filesystem::ofstream("./"+name+"/cubature.dat",ios::binary));
#else 
	sol.readCubature(boost::filesystem::ifstream("./"+name+"/cubature.dat",ios::binary));
#endif 

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
