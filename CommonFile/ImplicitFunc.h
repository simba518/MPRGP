#ifndef IMPLICIT_FUNC_H
#define IMPLICIT_FUNC_H

#include "GridOp.h"
#include "ObjMesh.h"
#include "CollisionDetection.h"
#include "SemiLagrangianContouring.h"

PRJ_BEGIN

class ImplicitFuncPlane : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncPlane();
    ImplicitFuncPlane(const Vec3& x0,const Vec3& n);
    ImplicitFuncPlane(const Vec3& a,const Vec3& b,const Vec3& c);
    virtual scalar operator()(const Vec3& pos) const;
    PlaneTpl<scalar> _p;
};
class ImplicitFuncCSG : public ImplicitFunc<scalar>
{
public:
    enum OP_TYPE {
        INTERSECT,
        UNION,
        SUBTRACT,
    };
    ImplicitFuncCSG(OP_TYPE op);
    ImplicitFuncCSG(const BBox<scalar,3>& bb,OP_TYPE op);
    ImplicitFuncCSG(const BBox<scalar,2>& bb,OP_TYPE op);
    virtual scalar operator()(const Vec3& pos) const;
    void setAlpha(const scalar& alpha);
    boost::shared_ptr<ImplicitFunc<scalar> > _a;
    boost::shared_ptr<ImplicitFunc<scalar> > _b;
    scalar _alpha;
    OP_TYPE _op;
};
class ImplicitFuncGridRef : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncGridRef(const Grid<scalar,scalar>& ls):_lsRef(ls) {}
    virtual scalar operator()(const Vec3& pos) const;
    const Grid<scalar,scalar>& _lsRef;
};
class ImplicitFuncGrid : public ImplicitFuncGridRef
{
public:
    ImplicitFuncGrid():ImplicitFuncGridRef(_ls) {}
    Grid<scalar,scalar> _ls;
};
class ImplicitFuncReinit : public ImplicitFuncGrid
{
public:
    ImplicitFuncReinit(const Grid<scalar,scalar> &tpl,const ImplicitFunc<scalar>& inner);
    virtual scalar operator()(const Vec3& pos) const;
    Grid<scalar,scalar> _ls;
};
class ImplicitFuncMesh2D : public ImplicitFuncGrid
{
public:
    ImplicitFuncMesh2D(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb);
protected:
    BBox<scalar> _bb;
    boost::shared_ptr<MeshSLC> _mesh;
    boost::shared_ptr<AABBvh> _bvh;
    const scalar _cellSz;
};
class ImplicitFuncMesh3D : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncMesh3D(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb);
    virtual scalar operator()(const Vec3& pos) const;
protected:
    BBox<scalar> _bb;
    boost::shared_ptr<MeshSLC> _mesh;
    boost::shared_ptr<AABBvh> _bvh;
    const scalar _cellSz;
};
class ImplicitFuncMesh3DAccurate : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncMesh3DAccurate(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb);
    virtual scalar operator()(const Vec3& pos) const;
protected:
    BBox<scalar> _bb;
    boost::shared_ptr<MeshSLC> _mesh;
    boost::shared_ptr<AABBvh> _bvh;
    OctreeDistance::TreeBuildHelper _helper;
    boost::shared_ptr<OctreeDistance> _octree;
    const scalar _cellSz;
};
class ImplicitFuncOffset : public ImplicitFunc<scalar>
{
public:
    scalar operator()(const Vec3& pos) const {
        return _inner->operator()(pos)+_off;
    }
    //data
    scalar _off;
    boost::shared_ptr<ImplicitFunc<scalar> > _inner;
};
class ImplicitFuncNegate : public ImplicitFunc<scalar>
{
public:
    scalar operator()(const Vec3& pos) const {
        return -_inner->operator()(pos);
    }
    //data
    boost::shared_ptr<ImplicitFunc<scalar> > _inner;
};
class ImplicitFuncRosy : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncRosy(const Vec3& origin,const Vec3& X,const scalar& step,const scalar& coef);
    virtual scalar operator()(const Vec3& pos) const;
    virtual scalar dist(const Vec2& p) const;
    virtual scalar y(const scalar& x) const =0;
    //axis
    scalar _step;
    scalar _coef;
    Vec3 _origin;
    Vec3 _X;
};

PRJ_END

#endif
