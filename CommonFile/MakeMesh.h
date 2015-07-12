#ifndef MAKE_MESH_H
#define MAKE_MESH_H

#include "ObjMesh.h"

PRJ_BEGIN

class MakeMesh
{
public:
    static void makeTet3D(ObjMesh& m);
    static void makeBox3D(ObjMesh& m,const Vec3& ext);
	static void makeDiscreteBox3D(ObjMesh& m,const Vec3& ext);
    static void makeBox2D(ObjMesh& m,const Vec3& ext);
    static void makeSphere3D(ObjMesh& m,const scalar& rad,const sizeType& slice);
    static void makeSphere2D(ObjMesh& m,const scalar& rad,const sizeType& slice);
    static void makeCapsule3D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice);
    static void makeCapsule2D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice);
    static void makeCylinder3D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice,const sizeType& sliceY=16,bool cap=true);
    static void makeTorus3D(ObjMesh& m,const scalar& rad1,const scalar& rad2,const sizeType& slice1,const sizeType& slice2);
    static void makeRing3D(ObjMesh& m,const scalar& rad1,const scalar& rad2,const scalar& rad3,const sizeType& slice);
	static void makeGrid(ObjMesh& m,const Vec2i& slice);
protected:
    static sizeType GI(sizeType i,sizeType j,sizeType si,sizeType sj){return (i%si)*sj+(j%sj);}
};

PRJ_END

#endif
