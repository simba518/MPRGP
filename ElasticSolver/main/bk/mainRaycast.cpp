#include "mainConfig.h"
#ifdef MAIN_RAYCAST_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 2
#define SZ 0.1f
//#define SINGLE_BALL
	ParticleSetN pSet;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-1;k<1;k+=SZ)
#endif
	for(scalar i=-1;i<1;i+=SZ)
	for(scalar j=-1;j<1;j+=SZ)
	{
		ParticleN<scalar> p,pp;
		p._pos=Vec3(i,j,k);
#ifdef SINGLE_BALL
		if(p._pos.norm() < 1)
			pSet.addParticle(p);
#else
		if(p._pos.norm() < 0.75)
		{
			pp=p;
			pp._pos.block<DIM,1>(0,0)-=
				Vec3::Constant(1).block<DIM,1>(0,0);
			pSet.addParticle(pp);

			pp=p;
			pp._pos.block<DIM,1>(0,0)+=
				Vec3::Constant(1).block<DIM,1>(0,0);
			pSet.addParticle(pp);
		}
#endif
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM);
	mesh.reset(DIM,0.2f,pSet);
	mesh.writeVTK("./mesh.vtk");

#define NR_TEST 100
#define BODYID 0
	boost::filesystem::create_directory("./Raycast");
	for(sizeType i=0;i<NR_TEST;i++)
	{
		Vec3 r=Vec3::Random();
#if DIM==2
		r[2]=0.0f;
#endif
		r=r.normalized()*3;
		LineSeg l(r,Vec3::Zero());
		FEMMesh::Interp cd,cd2;
		scalar dist,dist2;

		{
			dist=numeric_limits<scalar>::max();
			mesh.rayCastMeshBF(l,cd,dist,BODYID);
			dist2=numeric_limits<scalar>::max();
			mesh.rayCastMesh(l,cd2,dist2,BODYID);
			ASSERT(dist == dist2 && cd._id == cd2._id && cd._coef == cd2._coef);

			vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
			vss.push_back(r);
			if(dist == numeric_limits<scalar>::max())
				vss.push_back(Vec3::Zero());
			else vss.push_back(mesh.getC(cd._id).getVert(cd));

			ostringstream oss;oss << "./Raycast/fmesh" << i << ".vtk";
			VTKWriter<scalar> os("Raycast",oss.str(),true);
			os.appendPoints(vss.begin(),vss.end());
			os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
						   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
						   VTKWriter<scalar>::LINE);
		}

		{
			dist=numeric_limits<scalar>::max();
			mesh.rayCastPSetBF(l,SZ,cd,dist,BODYID);
			dist2=numeric_limits<scalar>::max();
			mesh.rayCastPSet(l,SZ,cd2,dist2,BODYID);
			ASSERT(dist == dist2 && cd._id == cd2._id && cd._coef == cd2._coef);

			ostringstream oss;oss << "./Raycast/fpset" << i << ".vtk";
			VTKWriter<scalar> os("Raycast",oss.str(),true);
			if(dist < numeric_limits<scalar>::max())
			{
				vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
				vss.push_back(mesh.getPSet()[cd._id]._pos);
				os.appendPoints(vss.begin(),vss.end());
				os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
							   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,0,0),
							   VTKWriter<scalar>::POINT);
			}
		}
	}
}
#endif
