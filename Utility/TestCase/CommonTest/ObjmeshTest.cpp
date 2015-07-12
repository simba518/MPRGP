#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <Objmesh.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(ObjmeshTest)

BOOST_AUTO_TEST_CASE(testObjLoadDino){

  const string fname = string(TEST_DATA_DIR)+"/dino.obj";
  Objmesh mesh;
  TEST_ASSERT( mesh.load(fname) );
  ASSERT_EQ (mesh.getVertsNum(),28098);
  ASSERT_EQ (mesh.getFacesNum(),56192);

  Vector3d v0,vT;
  v0 << 1.79623, -4.35718, 0.464406;
  vT << 3.51123,4.55572, -0.326623;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(0),v0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(mesh.getVertsNum()-1),vT,3,1e-6);

  Vector3d vn0,vnT;
  vn0 << 0.0862053, -0.99618,0.0139087;
  vnT << -0.0341522, 0.996979, 0.0697619;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(0),vn0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(mesh.getNormalNum()-1),vnT,3,1e-6); 

  Vector3i f0,fT;
  f0 << 0,1,2;
  fT << 28094,28090,28091;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getFaces(0),f0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getFaces(mesh.getFacesNum()-1),fT,3,1e-6);

  Vector3d Kd,Ka,Tf,Ks;
  Kd << 0.00,0.60,0.00;
  Ka << 0.00,0.10,0.00;
  Ks << 0.35,0.35,0.35;
  const double Ni = 1.00;
  const double Ns = 200;
  
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().diffuse,Kd,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().ambient,Ka,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().specular,Ks,3,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().shininess,Ns,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().ior,Ni,1e-6);
}

BOOST_AUTO_TEST_CASE(testObjLoadBeam){

  const string fname = string(TEST_DATA_DIR)+"beam.obj";
  Objmesh mesh;
  TEST_ASSERT( mesh.load(fname) );
  ASSERT_EQ (mesh.getVertsNum(),1280);
  ASSERT_EQ (mesh.getFacesNum(),822*2);
  ASSERT_EQ (mesh.getFacesNum()*3,mesh.getNormalIndex().size());

  Vector3d v0,vT;
  v0 << -12.48 ,-0.96 ,-2.24;
  vT << 12.48 ,0.96 ,2.24;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(0),v0,3,1e-2);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(mesh.getVertsNum()-1),vT,3,1e-2);

  Vector3d vn0,vnT;
  vn0 << -1 ,0 ,0;
  vnT << 1, 0, 0;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(0),vn0,3,1e-1);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(mesh.getNormalNum()-1),vnT,3,1e-1); 

  Vector3i f0,f1,fT1,fT;
  f0 << 0,1,9;
  f1 << 0,9,8;
  fT1 << 1271-1,1279-1,1280-1;
  fT << 1271-1,1280-1,1272-1;
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(0),f0,3);
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(1),f1,3);
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(mesh.getFacesNum()-2),fT1,3);
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(mesh.getFacesNum()-1),fT,3);

  Vector3i n0,n1,nT1,nT;
  n0 << 0,0,0;
  n1 << 0,0,0;
  nT1 << 822-1,822-1,822-1;
  nT << 822-1,822-1,822-1;
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(0),n0,3);
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(1),n1,3);
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(mesh.getFacesNum()-2),nT1,3);
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(mesh.getFacesNum()-1),nT,3);

  Vector3d Kd,Ka,Tf,Ks;
  Kd << 0.00,0.60,0.00;
  Ka << 0.00,0.10,0.00;
  Ks << 0.35,0.35,0.35;
  const double Ni = 1.00;
  const double Ns = 200;
  
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().diffuse,Kd,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().ambient,Ka,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().specular,Ks,3,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().shininess,Ns,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().ior,Ni,1e-6);

  // const string outfname = "./TestCase/TestData/temptBeam.obj";
  // TEST_ASSERT( mesh.write(outfname) );
}

BOOST_AUTO_TEST_CASE(testObjWriteDino){

  const string fname = string(TEST_DATA_DIR)+"/dino.obj";
  Objmesh mesh;
  TEST_ASSERT( mesh.load(fname) );
  ASSERT_EQ (mesh.getVertsNum(),28098);
  ASSERT_EQ (mesh.getFacesNum(),56192);
  ASSERT_EQ (mesh.getFacesNum()*3,mesh.getNormalIndex().size());

  const string outfname = string(TEST_DATA_DIR)+"/temptDino.obj";
  TEST_ASSERT( mesh.write(outfname) );
  TEST_ASSERT( mesh.load(outfname) );
  ASSERT_EQ (mesh.getVertsNum(),28098);
  ASSERT_EQ (mesh.getFacesNum(),56192);

  Vector3d v0,vT;
  v0 << 1.79623, -4.35718, 0.464406;
  vT << 3.51123,4.55572, -0.326623;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(0),v0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVerts(mesh.getVertsNum()-1),vT,3,1e-6);

  Vector3d vn0,vnT;
  vn0 << 0.0862053, -0.99618,0.0139087;
  vnT << -0.0341522, 0.996979, 0.0697619;
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(0),vn0,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL (mesh.getVertNormal(mesh.getVertsNum()-1),vnT,3,1e-6); 

  Vector3i f0,fT;
  f0 << 0,1,2;
  fT << 28094,28090,28091;
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(0),f0,3);
  ASSERT_EQ_SMALL_VEC (mesh.getFaces(mesh.getFacesNum()-1),fT,3);

  Vector3i n0,n1,nT1,nT;
  n0 << 1-1, 2-1, 3-1;
  nT << 28095-1, 28091-1, 28092-1;
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(0),n0,3);
  ASSERT_EQ_SMALL_VEC (mesh.getNormalIndex(mesh.getFacesNum()-1),nT,3);

  Vector3d Kd,Ka,Tf,Ks;
  Kd << 0.00,0.60,0.00;
  Ka << 0.00,0.10,0.00;
  Ks << 0.35,0.35,0.35;
  const double Ni = 1.00;
  const double Ns = 200;
  
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().diffuse,Kd,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().ambient,Ka,3,1e-6);
  ASSERT_EQ_SMALL_VEC_TOL(mesh.getMtl().specular,Ks,3,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().shininess,Ns,1e-6);
  ASSERT_EQ_TOL(mesh.getMtl().ior,Ni,1e-6);
}

BOOST_AUTO_TEST_CASE(testObjBeamVtKWrite){
  
  const string fname = string(TEST_DATA_DIR)+"beam.obj";
  Objmesh mesh;
  TEST_ASSERT( mesh.load(fname) );

  const string fnameout = string(TEST_DATA_DIR)+"tempt_beam_binary.vtk";
  mesh.writeVTK(fnameout,true);

  const string fnameout2 = string(TEST_DATA_DIR)+"tempt_beam_ascii.vtk";
  mesh.writeVTK(fnameout2,false);
}

BOOST_AUTO_TEST_CASE(testObjCombine){

  Objmesh mesh_beam;
  TEST_ASSERT( mesh_beam.load(string(TEST_DATA_DIR)+"beam.obj") );

  Objmesh mesh_dino;
  TEST_ASSERT( mesh_dino.load(string(TEST_DATA_DIR)+"dino.obj") );

  mesh_beam.combine(mesh_dino);
  TEST_ASSERT(mesh_beam.write(string(TEST_DATA_DIR)+"tempt_test_combine_dino_beam.obj"));

}

BOOST_AUTO_TEST_CASE(testObjToPlanes){
 
  Objmesh mesh_dino;
  TEST_ASSERT( mesh_dino.load(string(TEST_DATA_DIR)+"dino.obj") );
  VVec4d planes;
  mesh_dino.toPlanes(planes);
  ASSERT_EQ(planes.size(), mesh_dino.getFacesNum());

  ASSERT_EQ_TOL(planes[0].head(3).norm(),1.0f,1e-8);
  ASSERT_EQ_TOL(planes[planes.size()-1].head(3).norm(),1.0f,1e-8);

  for (int f = 0; f < planes.size(); ++f){
	const Vector3d v0 = mesh_dino.getVerts(mesh_dino.getFaces(f)[0]);
	const Vector3d v1 = mesh_dino.getVerts(mesh_dino.getFaces(f)[1]);
	const Vector3d v2 = mesh_dino.getVerts(mesh_dino.getFaces(f)[2]);

	ASSERT_EQ_TOL(planes[f].head(3).dot(v0)+planes[f][3],0.0f,1e-8);
	ASSERT_EQ_TOL(planes[f].head(3).dot(v1)+planes[f][3],0.0f,1e-8);
	ASSERT_EQ_TOL(planes[f].head(3).dot(v2)+planes[f][3],0.0f,1e-8);
  }

  VectorXd verts(3*3);
  verts << 0,0,0,
	1,0,0,
	0,0,1;
  VectorXi faces(3);
  faces << 0,1,2;
  Objmesh one_trangle(verts,faces);
  one_trangle.toPlanes(planes);
  ASSERT_EQ(planes.size(),1);
  Vector4d p;
  p << 0,-1,0,0;
  ASSERT_EQ_SMALL_VEC_TOL(planes[0],p,4,1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
