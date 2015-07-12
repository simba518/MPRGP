#ifndef _MESHRENDER_H_
#define _MESHRENDER_H_

#include <Objmesh.h>
#include <TetMesh.h>

namespace UTILITY{

  void draw(const Objmesh& obj,const ObjMtl&mtl);

  inline void draw(const Objmesh& obj){
	draw(obj,obj.getMtl());
  }

  inline void draw(pObjmesh_const obj,const ObjMtl&mtl){
	if (obj){
	  draw(*obj,mtl);
	}
  }

  inline void draw(pObjmesh_const obj){
	if (obj) draw(*obj,obj->getMtl());
  }

  void draw(const TetMesh &tetmesh,const ObjMtl&mtl,const double *u=NULL);

  inline void draw(const TetMesh &tetmesh,const double *u=NULL){
	static ObjMtl mtl; // use default material.
	draw(tetmesh,mtl,u);
  }

  inline void draw(pTetMesh_const tetmesh,const ObjMtl&mtl,const double *u=NULL){
	if (tetmesh) draw(*tetmesh,mtl,u);
  }

  inline void draw(pTetMesh_const tetmesh,const double *u=NULL){
	if(tetmesh) draw(*tetmesh,u);
  }

  

}

#endif /* _MESHRENDER_H_ */
