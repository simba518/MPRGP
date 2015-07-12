#ifndef STATIC_GEOM_H
#define STATIC_GEOM_H

#include "MathBasic.h"
#include "IO.h"

PRJ_BEGIN

template <typename T,typename BBOX> struct Node;
template <typename T> class ObjMeshTpl;
template <typename T,typename BBOX> class BVHQuery;
template <typename T,int dim> class OBBTpl;
template <typename T> class PlaneTpl;

struct StaticGeomCell : public Serializable {
    StaticGeomCell(sizeType type);
    StaticGeomCell(const Mat4& T,sizeType dim,sizeType type);
    virtual ~StaticGeomCell() {}
    virtual void getMesh(ObjMeshTpl<scalar>& mesh) const;
    virtual BBox<scalar> getBB() const;
    virtual bool dist(const Vec3& pt,Vec3& n) const;
	virtual void closest(const Vec3& pt,Vec3& n) const;
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    sizeType _index;
protected:
    virtual BBox<scalar> getBBInner() const=0;
    virtual bool distInner(const Vec3& pt,Vec3& n) const=0;
	virtual void closestInner(const Vec3& pt,Vec3& n) const;
    //data
    Mat4 _T,_invT;
    sizeType _dim;
};
class StaticGeom : public HasMagic
{
public:
    //method
    StaticGeom(sizeType dim);
    const vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > >& getBVH() const;
    sizeType nrG() const {
        return (sizeType)_css.size();
    }
    const StaticGeomCell& getG(sizeType i) const {
        return *(_css[i]);
    }
    boost::shared_ptr<StaticGeomCell> getGPtr(sizeType i) const {
        return _css[i];
    }
    StaticGeomCell& getG(sizeType i) {
        return *(_css[i]);
    }
    void assemble();
    //geometry
    void addGeomCell(boost::shared_ptr<StaticGeomCell> c);
    void addGeomBox(const Mat4& trans,const BBox<scalar>& bb);
    void addGeomBox(const Mat4& trans,const Vec3& ext);
    void addGeomBox(const OBBTpl<scalar,2>& obb);
    void addGeomBox(const OBBTpl<scalar,3>& obb);
    void addGeomCylinder(const Mat4& trans,scalar rad,scalar y);
    void addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext=100.0f);
    void addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext=100.0f);
    void addGeomSphere(const Vec3& ctr,scalar rad);
    void addGeomMesh(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth=0.0f,bool insideOut=false);
    void addGeomMesh(const Mat4& trans,const std::string& path,scalar depth=0.0f,bool insideOut=false);
    //IO
    void writeVTK(const std::string& path) const;
    void writeObj(const std::string& path) const;
    void writeBVH() const;
    bool write(std::ostream& os) const;
    bool read(std::istream& is);
private:
    //mesh data
    vector<boost::shared_ptr<StaticGeomCell> > _css;
    boost::shared_ptr<vector<Node<boost::shared_ptr<StaticGeomCell>,BBox<scalar> > > > _bvh;
    sizeType _dim;
};

PRJ_END

#endif