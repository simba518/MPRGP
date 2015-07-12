#include "mainConfig.h"
#ifdef MAIN_COUPLED_RS_BEAM

#include "../RSCoupledFEMSolver.h"
#include "../FEMProfiler.h"
#include "BeamModel.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	std::string name="beamCoupledRS";
	boost::filesystem::create_directory("./"+name);

	CoupledRotateReducedFEMSolver sol(2);
	sol.resetParam(0.0f,0.0f,1E7f,1000.0f,1000);
	buildBeam(sol,16.0f);
	sol.getMesh().calcMatDist();
	BBox<scalar> bb(Vec3(2.5f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
	buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,&bb,2500000.0f);
	sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	sol.getMesh().getPSet().writeVTK("./"+name+"/pSet.vtk");
	sol.buildU(30);
	sol.setUseRigid(true);
	FEMMesh mesh2=sol.getMesh();
	{
		Mat4 T=Mat4::Identity();
		T.block<3,1>(0,3)=-Vec3::Unit(1)*5.0f;
		mesh2.applyTrans(T,0,true,true);
		sol.getMesh()+=mesh2;

		FEMFrameData dat;
        boost::filesystem::ifstream is("./"+name+"/frm95.dat",ios::binary);
        dat.read(is);
		sol.setFrame(dat);
		sol.getMesh().writeVTK("./"+name+"/mesh.vtk");
	}

	sol._geom.reset(new FEMGeom(2));
	sol._geom->addGeomPlane(Mat4::Identity(),Plane(Vec3(0.0f,-10.0f,0.0f),Vec3(0.0f,1.0f,0.0f)));
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
		if(i == 50)
		{
			Mat4 T=Mat4::Identity();
			T.block<3,1>(0,3)=-Vec3::Unit(1)*5.0f;
			mesh2.applyTrans(T,0,true,true);
			sol.getMesh()+=mesh2;

			sol.clearEnergy();
			BBox<scalar> bb(Vec3(2.5f,0.0f,0.0f),Vec3(4.0f,1.0f,1.0f));
			buildBeamEnergy(sol,-50.0f,MaterialEnergy::LINEAR,&bb,2500000.0f);
		}
		if(i == 95)
		{
			FEMFrameData dat;
			sol.getFrame(dat);
            boost::filesystem::ofstream os("./"+name+"/frm95.dat",ios::binary);
            dat.write(os);
		}
		ostringstream ossm;ossm << "./"+name+"/frm/m" << i << ".vtk";
		sol.getMesh().writeVTK(ossm.str());
	}
    return 0;
}
#endif
