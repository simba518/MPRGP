#include "mainConfig.h"
#ifdef MAIN_SELFCD_TEST

#include "../FEMMesh.h"
#include "../FEMGeom.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 2
#define BID 0
#define SZ 0.1f
#define RAD 0.5f
	ParticleSetN pSet;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-RAD;k<RAD;k+=SZ)
#endif
	for(scalar i=-RAD;i<RAD;i+=SZ)
	for(scalar j=-RAD;j<RAD;j+=SZ)
	{
		Vec3 dv=Vec3::Zero();
		dv.block(0,0,DIM,1).setConstant(RAD);

		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
		{
			pSet.addParticle(p);
			p._pos+=dv*2.5f;
			pSet.addParticle(p);
		}
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM);
	mesh.reset(DIM,0.2f,pSet);
	Mat4 T=Mat4::Identity();
	T.block<3,3>(0,0)=makeRotation<scalar>(Vec3(0.0f,0.0f,rand()));
	mesh.applyTrans(T,0,true,true);
	mesh.writeVTK("./mesh.vtk");

	boost::filesystem::create_directory("./CollSelfCD/");
	sizeType j=0;
	for(scalar i=0;i<RAD*5.0f;i+=0.02f)
	{
		Vec3 dv=Vec3::Zero();
		dv.block(0,0,DIM,1).setConstant(i);
		for(sizeType c=0;c<(sizeType)mesh.nrC();c++)
			if(mesh.getC(c)._bodyId==BID)
			for(sizeType vid=0;vid<DIM+1;vid++)
			{
				boost::shared_ptr<FEMMesh::Vertex> v=mesh.getC(c)._v[vid];
				v->_pos=v->_pos0+dv;
			}
		mesh.updateMesh();
		
		ostringstream oss;oss << "./CollSelfCD/mesh" << j << ".vtk";
		mesh.writeVTK(oss.str());

		ostringstream oss2;oss2 << "./CollSelfCD/pSet" << j << ".vtk";
		mesh.getPSet().writeVTK(oss2.str());
		
		ostringstream oss3;oss3 << "./CollSelfCD/coll" << j << ".vtk";
		FEMMesh::DebugCollider coll(oss3.str(),DIM);
		mesh.collideMesh(mesh,coll,true);
		j++;
	}
}
#endif
