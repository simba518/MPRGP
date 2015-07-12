#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <TetMesh.h>
#include <TetMeshEmbeding.h>
#include "TetMeshFactoryForTest.h"
using namespace std;
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(TetMeshTest)

BOOST_AUTO_TEST_CASE(initializeTest){

  // create a simple mesh with two tets with five nodes
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  TEST_ASSERT (tet_mesh != NULL);

  // check nodes, tets
  ASSERT_EQ (tet_mesh->nodes().size(), 5);
  ASSERT_EQ (tet_mesh->tets().size(), 2);
  { // check face ids
  	const FaceId &faceId = tet_mesh->faceId();
  	ASSERT_EQ (faceId.size(), 7);
  	FaceId::const_iterator iter;
  	iter = faceId.find(HashedId(1,2,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,2,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,3,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,2,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,0);
  	ASSERT_EQ (iter->second.second,1);

  	iter=faceId.find(HashedId(1,2,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,2,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,4,0));
  	TEST_ASSERT ( iter != faceId.end() );
  	ASSERT_EQ (iter->second.first,1);
  	ASSERT_EQ (iter->second.second,-1);

  	iter=faceId.find(HashedId(0,1,8,0));
  	TEST_ASSERT ( iter == faceId.end() );
  }
  
  { // check neighors of tets sharing the same faces
  	const VVec4i &faceNeighTet = tet_mesh->faceNeighTet();
  	ASSERT_EQ (faceNeighTet.size(),2);
  	ASSERT_EQ ( faceNeighTet[0], Vector4i(-1,-1,-1,1));
  	ASSERT_EQ ( faceNeighTet[1], Vector4i(-1,-1,-1,0));
  }
  
  { // check neighors of nodes share the same edges.
  	const VectorUseti &nodeNeighNode = tet_mesh->nodeNeighNode();

  	ASSERT_EQ( nodeNeighNode.size(), 5 );
  	ASSERT_EQ( nodeNeighNode[0].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[1].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[2].size(), 4 );
  	ASSERT_EQ( nodeNeighNode[3].size(), 3 );
  	ASSERT_EQ( nodeNeighNode[4].size(), 3 );

  	TEST_ASSERT( nodeNeighNode[0].find(1) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(2) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(3) != nodeNeighNode[0].end() );
  	TEST_ASSERT( nodeNeighNode[0].find(4) != nodeNeighNode[0].end() );
  }

  // test empty mesh -----------------------------------------------------------
  VVec3d nodes;
  VVec4i tets;
  tet_mesh->reset(nodes, tets);
  ASSERT_EQ (tet_mesh->nodes().size(), 0);
}

BOOST_AUTO_TEST_CASE(IOTest){

  const std::string fname = std::string(TEST_DATA_DIR)+"/fish.abq";
  TetMesh tetMesh;
  TEST_ASSERT(tetMesh.load(fname));
  ASSERT_EQ(tetMesh.nodes().size(),885);
  ASSERT_EQ(tetMesh.tets().size(),3477);

  Vector3d v;
  v << 0.183433000000000, -0.316178000000000, -0.303435000000000;
  ASSERT_EQ_SMALL_VEC_TOL(tetMesh.nodes()[0],v,3,1e-12);

  v << -0.075092799520637,  0.212514873566465,  -0.485906600082581;
  ASSERT_EQ_SMALL_VEC_TOL(tetMesh.nodes()[884],v,3,1e-12);

  Vector4i t;
  t << 572-1, 270-1, 373-1, 39-1;
  ASSERT_EQ_SMALL_VEC_TOL(tetMesh.tets()[0],t,4,1e-12);
  t << 675-1, 479-1, 473-1, 689-1;
  ASSERT_EQ_SMALL_VEC_TOL(tetMesh.tets()[3476],t,4,1e-12);
}

BOOST_AUTO_TEST_CASE(VTK_IOTest){

  const std::string fname = std::string(TEST_DATA_DIR)+"/fish.abq";
  TetMesh tetMesh;
  TEST_ASSERT(tetMesh.load(fname));
  ASSERT_EQ(tetMesh.nodes().size(),885);
  ASSERT_EQ(tetMesh.tets().size(),3477);
  const std::string vtkfname = std::string(TEST_DATA_DIR)+"/tempt_fish_ascii.vtk";
  TEST_ASSERT(tetMesh.writeVTK(vtkfname,false));

  const std::string vtkfname2 = std::string(TEST_DATA_DIR)+"/tempt_fish_binary.vtk";
  TEST_ASSERT(tetMesh.writeVTK(vtkfname2,true));


  TetMesh tetMeshB;
  const std::string fnamebeam = std::string(TEST_DATA_DIR)+"/beam.abq";
  TEST_ASSERT(tetMeshB.load(fnamebeam));

  const std::string beamfout = std::string(TEST_DATA_DIR)+"/tempt_beamtet_binary.vtk";
  TEST_ASSERT(tetMeshB.writeVTK(beamfout,true));
}

BOOST_AUTO_TEST_CASE(testContainingEle){
 
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  Vector3d pos;

  pos << -1,0,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),-1);

  pos << 0,0,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),0);

  pos << 0,-1,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),1);

  pos << 0.1,-0.1,0.1;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),1);

  pos << 0,-6,0;
  ASSERT_EQ (tet_mesh->getContainingElement(pos),-1);
}

BOOST_AUTO_TEST_CASE(testClosestEle){

  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  Vector3d pos;

  pos << 10,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << -1,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,0,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,-1,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),1);

  pos << 0.1,-0.1,0.1;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),0);

  pos << 0,-6,0;
  ASSERT_EQ (tet_mesh->getClosestElement(pos),1);
}

BOOST_AUTO_TEST_CASE(testInterp){
 
  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();

  VectorXd vertices;
  vector<int> nodes;
  VectorXd weights; 
  tet_mesh->buildInterpWeights(vertices,nodes,weights);

  vertices.resize(9);
  vertices << 0,-1,0,  0.1,0.1,0.1, -1,-10,1;
  tet_mesh->buildInterpWeights(vertices,nodes,weights);

  VectorXd u(5*3);
  u << 1,0,0, 1,0,0, 1,0,0, 1,0,0, 1,0,0;

  VectorXd uTarget;
  tet_mesh->interpolate(nodes,weights,u,uTarget);
  ASSERT_EQ(uTarget.size(),9);
  VectorXd correct_uTarget(9);
  correct_uTarget << 1,0,0, 1,0,0, 1,0,0;
  ASSERT_EQ_SMALL_VEC_TOL(uTarget,correct_uTarget,uTarget.size(),1e-12);

}

BOOST_AUTO_TEST_CASE(testInterpMatrix){

  pTetMesh tet_mesh = TetMeshFactoryForTest::tet2();
  VectorXd vertices;
  vector<int> nodes;
  VectorXd weights; 
  vertices.resize(9);
  vertices << 0,-1,0,  0.1,0.1,0.1, -1,-10,1;
  tet_mesh->buildInterpWeights(vertices,nodes,weights);
  SparseMatrix<double> A;
  tet_mesh->buildInterpMatrix(nodes,weights,tet_mesh->nodes().size(),A);

  VectorXd u(5*3);
  u << 1,0,0, 1,0,0, 1,0,0, 1,0,0, 1,0,0;

  VectorXd uTarget_1;
  tet_mesh->interpolate(nodes,weights,u,uTarget_1);

  ASSERT_EQ(A.rows(),uTarget_1.size());
  ASSERT_EQ(A.cols(),u.size());
  const VectorXd uTarget_2 = A*u;
  ASSERT_EQ_SMALL_VEC_TOL(uTarget_1,uTarget_2,uTarget_2.size(),1e-12);
}

BOOST_AUTO_TEST_CASE(testInterpIODino){

  const string tetfname = std::string(TEST_DATA_DIR)+"dino.abq";
  const string objfname = std::string(TEST_DATA_DIR)+"dino.obj";
  const string weightsfname = std::string(TEST_DATA_DIR)+"dinoW.txt";
  TetMeshEmbeding volobj;

  TEST_ASSERT(volobj.loadTetMesh(tetfname));
  TEST_ASSERT(volobj.loadObjMesh(objfname));
  TEST_ASSERT(volobj.loadWeights(weightsfname));

  ASSERT_EQ(volobj.getInterpNodes().size(),28098*4);
  ASSERT_EQ(volobj.getInterpWeights().size(),28098*4);

  ASSERT_EQ(volobj.getInterpNodes()[0],986);
  ASSERT_EQ_TOL(volobj.getInterpWeights()[0],0.068584,1e-6);

  ASSERT_EQ(volobj.getInterpNodes()[28098*4-1],1343);
  ASSERT_EQ_TOL(volobj.getInterpWeights()[28098*4-1],0.282485,1e-6);
}

BOOST_AUTO_TEST_CASE(testInterpIOBeam){

  const string tetfname = std::string(TEST_DATA_DIR)+"beam.abq";
  const string objfname = std::string(TEST_DATA_DIR)+"beam.obj";
  const string weightsfname = std::string(TEST_DATA_DIR)+"beamW.txt";
  TetMeshEmbeding volobj;

  TEST_ASSERT(volobj.loadTetMesh(tetfname));
  TEST_ASSERT(volobj.loadObjMesh(objfname));
  TEST_ASSERT(volobj.loadWeights(weightsfname));
  
  ASSERT_EQ(volobj.getInterpNodes().size(),1280*4);
  ASSERT_EQ(volobj.getInterpWeights().size(),1280*4);

  ASSERT_EQ(volobj.getInterpNodes()[0],147);
  ASSERT_EQ_TOL(volobj.getInterpWeights()[0],-0.0,1e-8);

  ASSERT_EQ(volobj.getInterpNodes()[1280*4-1],136);
  ASSERT_EQ_TOL(volobj.getInterpWeights()[1280*4-1],0.0,1e-6);

}

BOOST_AUTO_TEST_SUITE_END()
