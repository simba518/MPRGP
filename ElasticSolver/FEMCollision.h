#ifndef FEM_COLLISION_H
#define FEM_COLLISION_H

#include "MathBasic.h"
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

template <typename T>
struct Node;
template <typename T>
class LineSegTpl;
struct FEMVertex;
struct FEMCell;
struct FEMBody;
struct FEMInterp;
struct FEMGeomCell;
class FEMCollider;
class FEMMesh;
class FEMGeom;

template <typename A,typename B>
struct FEMCache {
    FEMCache() {}
    FEMCache(boost::shared_ptr<A> a,boost::shared_ptr<B> b):_A(a),_B(b) {}
    FEMCache(boost::shared_ptr<FEMBody> bA,boost::shared_ptr<A> a,boost::shared_ptr<FEMBody> bB,boost::shared_ptr<B> b):_bA(bA),_A(a),_bB(bB),_B(b) {}
    bool operator!=(const FEMCache& other) const {
        return !operator==(other);
    }
    bool operator==(const FEMCache& other) const {
        return _A==other._A && _B==other._B;
    }
    bool operator<(const FEMCache& other) const {
        return _A<other._A || (_B==other._B && _B<other._B);
    }
    boost::shared_ptr<FEMBody> _bA,_bB;
    boost::shared_ptr<A> _A;
    boost::shared_ptr<B> _B;
};
struct FEMSCache {
    FEMSCache() {}
    bool operator!=(const FEMSCache& other) const {
        return !operator==(other);
    }
    bool operator==(const FEMSCache& other) const {
        for(sizeType i=0; i<4; i++)
            if(_v[i]!=other._v[i])return false;
        return true;
    }
    bool operator<(const FEMSCache& other) const {
        for(sizeType i=0; i<4; i++)
            if(_v[i] < other._v[i])
                return true;
            else if(other._v[i] < _v[i])
                return false;
        return false;
    }
    boost::shared_ptr<FEMBody> _b[4];
    boost::shared_ptr<FEMVertex> _v[4];
    Vec3 _n;
    scalar _depth;
};
class FEMCollision
{
public:
    virtual ~FEMCollision() {}
    void setMesh(FEMMesh& mesh) {
        _mesh=&mesh;
    }
    sizeType dim() const;
    virtual boost::shared_ptr<FEMBody> createBody() const;
    virtual boost::shared_ptr<FEMCollision> copy() const;
    virtual void updateMesh();
    virtual void rayCastPSet(const LineSegTpl<scalar>& l,scalar rad,FEMInterp& cd,scalar& dist) const;
    virtual void rayCastMesh(const LineSegTpl<scalar>& l,FEMInterp& cd,scalar& dist) const;
    virtual void collideMesh(FEMCollider& coll,bool surface);
    virtual void collideGeom(const FEMGeom& other,FEMCollider& coll,bool surface);
protected:
    FEMMesh* _mesh;
    std::vector<FEMCache<FEMCell,FEMVertex> > _cache;
};
class BVHFEMCollision : public FEMCollision
{
public:
    typedef std::vector<Node<boost::shared_ptr<FEMBody> > > BVH;
    virtual boost::shared_ptr<FEMBody> createBody() const;
    virtual boost::shared_ptr<FEMCollision> copy() const;
    virtual void updateMesh();
    virtual void rayCastPSet(const LineSegTpl<scalar>& l,scalar rad,FEMInterp& cd,scalar& dist) const;
    virtual void rayCastMesh(const LineSegTpl<scalar>& l,FEMInterp& cd,scalar& dist) const;
    virtual void collideMesh(FEMCollider& coll,bool surface);
    virtual void collideGeom(const FEMGeom& other,FEMCollider& coll,bool surface);
protected:
    std::vector<FEMCache<FEMVertex,FEMGeomCell> > _cacheG;
    boost::shared_ptr<BVH> _bvh;
};
class SBVHFEMCollision : public BVHFEMCollision
{
public:
    virtual boost::shared_ptr<FEMBody> createBody() const;
    virtual boost::shared_ptr<FEMCollision> copy() const;
    virtual void collideMesh(FEMCollider& coll,bool surface);
protected:
    std::vector<FEMSCache> _cacheS;
};

PRJ_END

#endif