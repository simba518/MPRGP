#ifndef _TETMESHFACTORYFORTEST_H_
#define _TETMESHFACTORYFORTEST_H_

#include <eigen3/Eigen/Dense>
#include <TetMesh.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;

/**
 * @class TetMeshFactoryForTest generate tet meshes for test.
 * 
 */
class TetMeshFactoryForTest{
	
public:
  static pTetMesh tet2(){ 

	// create a simple mesh with two tets with five nodes
	VVec3d nodes;
	VVec4i tets;
	nodes.push_back(Vector3d(0,0,1));
	nodes.push_back(Vector3d(1,0,0));
	nodes.push_back(Vector3d(0,0,0));
	nodes.push_back(Vector3d(0,1,0));
	nodes.push_back(Vector3d(0,-5,0));

	tets.push_back(Vector4i(0,1,2,3));
	tets.push_back(Vector4i(0,2,1,4));

	pTetMesh tet_mesh= pTetMesh(new TetMesh(nodes, tets));
	tet_mesh->setSingleMaterial(50.0f,0.1f,0.45f);
	return tet_mesh;
  }

  static pTetMesh tet1(){

	// create a simple mesh with two tets with five nodes
	VVec3d nodes;
	VVec4i tets;
	nodes.push_back(Vector3d(0,0,1));
	nodes.push_back(Vector3d(1,0,0));
	nodes.push_back(Vector3d(0,0,0));
	nodes.push_back(Vector3d(0,1,0));

	tets.push_back(Vector4i(0,1,2,3));

	pTetMesh tet_mesh= pTetMesh(new TetMesh(nodes, tets));
	tet_mesh->setSingleMaterial(50.0f,0.5,0.45);
	return tet_mesh;
  }
};

#endif /* _TETMESHFACTORYFORTEST_H_ */
