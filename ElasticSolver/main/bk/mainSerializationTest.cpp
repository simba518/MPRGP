#include "mainConfig.h"
#ifdef MAIN_SERIALIZATION_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 3
#define SZ 0.05f
#define RAD 0.5f
	FEMMesh mesh(DIM);
#define WRITE
#ifdef WRITE
	ParticleSetN pSetBall;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-RAD*0.91f;k<RAD*0.9f;k+=SZ)
#endif
	for(scalar i=-RAD*0.91f;i<RAD*0.9f;i+=SZ)
	for(scalar j=-RAD*0.91f;j<RAD*0.9f;j+=SZ)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
			pSetBall.addParticle(p);
	}
#define SZS 5.0f
#define SP 15
	ParticleSetN pSet;
#if DIM==2
#else
	for(scalar k=0;k<SP;k+=SZS)
#endif
	for(scalar i=0;i<SP;i+=SZS)
	for(scalar j=0;j<SP;j+=SZS)
		for(sizeType p=0;p<pSetBall.size();p++)
		{
			ParticleN<scalar> pp=pSetBall[p];
			pp._pos+=Vec3(i,j,k);
			pSet.addParticle(pp);
		}
	pSet.writeVTK("./pset.vtk");
	mesh.reset(DIM,0.2f,pSet);
    boost::filesystem::ofstream os("f:/mesh.dat",ios::binary);
    mesh.write(os);
#endif
    boost::filesystem::ifstream is("f:/mesh.dat",ios::binary);
    mesh.read(is);
	FEMMesh tmp2=mesh;
	for(sizeType i=0;i<tmp2.nrV();i++)
		tmp2.getV(i)._pos+=Vec3::Constant(40.0f);
	tmp2.updateMesh();

	mesh+=tmp2;
	mesh.writeVTK("./mesh.vtk");
	mesh.writeBVH();
}
#endif
