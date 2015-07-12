#ifndef _TETMESHFACTORYFORTEST_H_
#define _TETMESHFACTORYFORTEST_H_

#include <boost/shared_ptr.hpp>
#include <VolMeshTet.h>
using namespace ELASTIC_BASE;

namespace TEST{
  
  /**
   * @class TetMeshFactoryForTest generate tet meshes for testing.
   * 
   */
  class TetMeshFactoryForTest{
	
  public:
	static pVolMeshTet tet2(){ 

	  // create a simple mesh with two tets with five nodes
	  VectorV3 nodes;
	  VectorV4i tets;
	  nodes.push_back(vec3(0,0,1));
	  nodes.push_back(vec3(1,0,0));
	  nodes.push_back(vec3(0,0,0));
	  nodes.push_back(vec3(0,1,0));
	  nodes.push_back(vec3(0,-5,0));

	  tets.push_back(vec4i(0,1,2,3));
	  tets.push_back(vec4i(0,2,1,4));

	  pVolMeshTet tet_mesh= pVolMeshTet(new VolMeshTet(nodes, tets));
	  const vec2 lame=ElasticMaterial::toLameConstant(vec2(5e5,0.45));
	  tet_mesh->setSingleMaterial(50.0f,lame[1],lame[0]);
	  return tet_mesh;
	}

	static pVolMeshTet tet1(){

	  // create a simple mesh with two tets with five nodes
	  VectorV3 nodes;
	  VectorV4i tets;
	  nodes.push_back(vec3(0,0,1));
	  nodes.push_back(vec3(1,0,0));
	  nodes.push_back(vec3(0,0,0));
	  nodes.push_back(vec3(0,1,0));

	  tets.push_back(vec4i(0,1,2,3));

	  pVolMeshTet tet_mesh= pVolMeshTet(new VolMeshTet(nodes, tets));
	  const vec2 lame=ElasticMaterial::toLameConstant(vec2(5e5,0.45));
	  tet_mesh->setSingleMaterial(50.0f,lame[1],lame[0]);
	  return tet_mesh;
	}
	
  };
  
  typedef boost::shared_ptr<TetMeshFactoryForTest> pTetMeshFactoryForTest;
  
}//end of namespace

#endif /*_TETMESHFACTORYFORTEST_H_*/
