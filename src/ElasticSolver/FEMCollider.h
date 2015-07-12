#ifndef FEM_COLLIDER_H
#define FEM_COLLIDER_H

#include "FEMEnergy.h"
#include "IO.h"
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

struct FEMVertex;
struct FEMCell;
class FEMMesh;

//Collider
class FEMCollider
{
public:
    virtual ~FEMCollider() {}
    virtual void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV) {
        ASSERT_MSG(false,"Not Implemented!")
    }
    virtual void handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary) {
        ASSERT_MSG(false,"Not Implemented!")
    }
    virtual void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n) {
        ASSERT_MSG(false,"Not Implemented!")
    }
};
class DebugFEMCollider : public FEMCollider
{
public:
    DebugFEMCollider(const std::string& path,sizeType dim);
    virtual ~DebugFEMCollider();
    virtual void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV);
    virtual void handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary);
    virtual void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);
private:
    VTKWriter<scalar> _os;
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > _vss;
    std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> > _issv;
    std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> > _issp;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > _issl;
    std::vector<scalar> _css;
    sizeType _dim;
};
class DefaultFEMCollider : public FEMCollider
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    DefaultFEMCollider(const FEMMesh& mesh,scalar K,scalar CF,scalar CH);
    const Vec& getFE() const;
    void getHE(Eigen::SparseMatrix<scalarD,0,sizeType>& H);
    virtual void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV);
    virtual void handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary);
    virtual void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);
protected:
    Vec _fE;
    TRIPS _HTrips;
    const FEMMesh& _mesh;
    scalar _K,_CF,_CH;
};

PRJ_END

#endif