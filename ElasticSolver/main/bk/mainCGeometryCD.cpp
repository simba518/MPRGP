#include "mainConfig.h"
#ifdef MAIN_COMPLEX_GEOMETRYCD_TEST

#include "../FEMGeom.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	FEMGeom geom(3);
	geom.addGeomMesh(Mat4::Identity(),"./bunny.obj",0.1f);
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
	for(scalarD i=0,j=0;i<RAD*5.0f;i+=0.01f)
	{
		Vec3 dv=Vec3::Zero();
		dv.setConstant(i);
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
