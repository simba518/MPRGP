#ifndef FEM_GEOM_H
#define FEM_GEOM_H

#include "MathBasic.h"
#include "IO.h"

PRJ_BEGIN

template <typename T> struct Node;
template <typename T> class ObjMeshTpl;
template <typename T> class BVHQuery;
template <typename T,int dim> class OBBTpl;
template <typename T> class PlaneTpl;

struct FEMGeomCell : public Serializable {
    FEMGeomCell(sizeType type);
    FEMGeomCell(const Mat4& T,sizeType dim,sizeType type);
    virtual ~FEMGeomCell() {}
    virtual void getMesh(ObjMeshTpl<scalar>& mesh) const;
    virtual BBox<scalar> getBB() const;
    virtual bool dist(const Vec3& pt,Vec3& n) const;
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    sizeType _index;
protected:
    virtual BBox<scalar> getBBInner() const=0;
    virtual bool distInner(const Vec3& pt,Vec3& n) const=0;
    //data
    Mat4 _T,_invT;
    sizeType _dim;
};
class FEMGeom : public HasMagic
{
public:
    //method
    FEMGeom(sizeType dim);
    const vector<Node<boost::shared_ptr<FEMGeomCell> > >& getBVH() const;
    sizeType nrG() const {
        return (sizeType)_css.size();
    }
    const FEMGeomCell& getG(sizeType i) const {
        return *(_css[i]);
    }
    boost::shared_ptr<FEMGeomCell> getGPtr(sizeType i) const {
        return _css[i];
    }
    FEMGeomCell& getG(sizeType i) {
        return *(_css[i]);
    }
    void assemble();
    //geometry
    void addGeomCell(boost::shared_ptr<FEMGeomCell> c);
    void addGeomBox(const Mat4& trans,const BBox<scalar>& bb);
    void addGeomBox(const Mat4& trans,const Vec3& ext);
    void addGeomBox(const OBBTpl<scalar,2>& obb);
    void addGeomBox(const OBBTpl<scalar,3>& obb);
    void addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext=100.0f);
    void addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext=100.0f);
    void addGeomSphere(const Vec3& ctr,scalar rad);
    void addGeomMesh(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth=0.0f,bool insideOut=false);
    void addGeomMesh(const Mat4& trans,const std::string& path,scalar depth=0.0f,bool insideOut=false);
    //IO
    void writeVTK(const std::string& path) const;
    void writeBVH() const;
    bool write(std::ostream& os) const;
    bool read(std::istream& is);
private:
    //mesh data
    vector<boost::shared_ptr<FEMGeomCell> > _css;
    boost::shared_ptr<vector<Node<boost::shared_ptr<FEMGeomCell> > > > _bvh;
    sizeType _dim;
};

PRJ_END

#endif