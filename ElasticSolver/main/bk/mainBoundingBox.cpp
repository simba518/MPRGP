#include "mainConfig.h"
#ifdef MAIN_BOUNDING_BOX_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 3
#define SZ 0.25f
#define RAD 0.5f
	ParticleSetN pSetBall;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-RAD;k<RAD;k+=SZ)
#endif
	for(scalar i=-RAD;i<RAD;i+=SZ)
	for(scalar j=-RAD;j<RAD;j+=SZ)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
			pSetBall.addParticle(p);
	}
#define SZS 2.0f
#define SP 15.0f
	ParticleSetN pSet;
#if DIM==2
#else
	for(scalar k=0;k<SP;k+=SZS)
#endif
	for(scalar i=0;i<SP;i+=SZS)
	for(scalar j=0;j<SP;j+=SZS)
	{
		Vec3 rand=Vec3::Random()*0.25f;
		for(sizeType p=0;p<pSetBall.size();p++)
		{
			ParticleN<scalar> pp=pSetBall[p];
			pp._pos+=Vec3(i,j,k)+rand;
			pSet.addParticle(pp);
		}
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM);
	mesh.reset(DIM,0.2f,pSet);
	mesh.writeVTK("./mesh.vtk");
	mesh.writeBVH();
}
#endif
