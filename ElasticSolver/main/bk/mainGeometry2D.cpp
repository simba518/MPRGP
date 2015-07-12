#include "mainConfig.h"
#ifdef MAIN_GEOMETRY2D_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	FEMGeom geom(2);
	//box
	scalar theta=rand()*M_PI/(scalar)RAND_MAX;
	Mat2 r;r << cos(theta),sin(theta),-sin(theta),cos(theta);
	geom.addGeomBox(OBB2D(r,Vec2::Random(),Vec2::Random()+Vec2::Ones()));
	//plane
	Vec3 x0=Vec3::Random();x0[2]=0.0f;
	Vec3 n=Vec3::Random();n[2]=0.0f;
	geom.addGeomPlane(Mat4::Identity(),Plane(x0,n));
	//sphere
	for(scalar i=10;i<30;i+=5)
	for(scalar j=10;j<30;j+=5)
		geom.addGeomSphere(Vec3(i,j,0.0f),1.5f);
	geom.assemble();
	//write
    boost::filesystem::ofstream os("./geom.dat",ios::binary);
    geom.write(os);
	
	FEMGeom geom2(2);
    boost::filesystem::ifstream is("./geom.dat",ios::binary);
    geom2.read(is);
	geom2.writeVTK("./geom.vtk");
	geom2.writeBVH();
}
#endif
