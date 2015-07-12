#include "mainConfig.h"
#ifdef MAIN_GEOMETRY3D_TEST

#include "../FEMGeom.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	FEMGeom geom(3);
	//box
	geom.addGeomBox(OBB3D(makeRotation<scalar>(Vec3::Random()),
						  Vec3::Random(),
						  BBox<scalar>(Vec3::Random()-Vec3::Ones(),
									   Vec3::Random()+Vec3::Ones())));
	//plane
	geom.addGeomPlane(Mat4::Identity(),Vec4(0,1,1,0.0f));
	geom.addGeomPlane(Mat4::Identity(),Plane(Vec3(0,50,50),Vec3(0,-1,-1)));
	//sphere
	for(scalar i=10;i<30;i+=5)
	for(scalar j=10;j<30;j+=5)
	for(scalar k=10;k<30;k+=5)
		geom.addGeomSphere(Vec3(i,j,k),1.5f);
	geom.addGeomMesh("./bunny.obj",0.1f);
	geom.assemble();
	//write
	geom.write(boost::filesystem::ofstream("./geom.dat",ios::binary));
	
	FEMGeom geom2(3);
    boost::filesystem::ifstream is("./geom.dat",ios::binary);
    geom2.read(is);
	geom2.writeVTK("./geom.vtk");
	geom2.writeBVH();
}
#endif
