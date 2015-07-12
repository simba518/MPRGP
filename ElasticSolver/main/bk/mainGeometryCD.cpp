#include "mainConfig.h"
#ifdef MAIN_GEOMETRYCD_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	FEMGeom geom(3);
	Mat4 T=Mat4::Identity();
	T.block<3,3>(0,0)=makeRotation<scalar>(Vec3::Random());
	T.block<3,1>(0,3).setConstant(100.0f);
	geom.addGeomBox(T,Vec3::Constant(20.0f));
	geom.addGeomPlane(Mat4::Identity(),Vec4(0.0f,0.0f,-1.0f,0.0f));
	geom.addGeomSphere(Vec3::Constant(50.0f),50.0f);
	geom.assemble();
	geom.writeVTK("./geom.vtk");

#define SZ 0.1f
#define RAD 0.5f
	ParticleSetN pSet;
	for(scalar k=-RAD;k<RAD;k+=SZ)
	for(scalar i=-RAD;i<RAD;i+=SZ)
	for(scalar j=-RAD;j<RAD;j+=SZ)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
		{
			p._pos.array()-=RAD*5.0f;
			pSet.addParticle(p);
		}
	}
	pSet.writeVTK("./pset.vtk");
	FEMMesh mesh(3);
	mesh.reset(3,0.2f,pSet);
	mesh.writeVTK("./mesh.vtk");

	boost::filesystem::create_directory("./CollGeomCD");
	for(scalarD i=0,j=0;i<RAD*500.0f;i+=1.0f)
	{
		Vec3 dv=Vec3::Zero();
		dv.block(0,0,3,1).setConstant(i);
		for(sizeType c=0;c<(sizeType)mesh.nrV();c++){
			FEMMesh::Vertex& v=mesh.getV(c);
			v._pos=v._pos0+dv;
		}
		mesh.updateMesh();
	
		ostringstream oss;oss << "./CollGeomCD/mesh" << j << ".vtk";
		mesh.writeVTK(oss.str());

		ostringstream oss2;oss2 << "./CollGeomCD/pSet" << j << ".vtk";
		mesh.getPSet().writeVTK(oss2.str());
		
		ostringstream oss3;oss3 << "./CollGeomCD/coll" << j << ".vtk";
		FEMMesh::DebugCollider coll(oss3.str(),3);
		geom.collideMesh(mesh,coll);
		j++;
	}
}
#endif
