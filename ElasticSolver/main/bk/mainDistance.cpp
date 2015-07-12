#include "mainConfig.h"
#ifdef MAIN_DISTANCE_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 3
#define SZ 0.05f
	ParticleSetN pSet;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-1;k<1;k+=SZ)
#endif
	for(scalar i=-1;i<1;i+=SZ)
	for(scalar j=-1;j<1;j+=SZ)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < 1)
			pSet.addParticle(p);
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM);
	mesh.reset(DIM,SZ*2,pSet);
	mesh.calcMatDist();
	mesh.writeVTK("./mesh.vtk");
}
#endif
