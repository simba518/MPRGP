#ifndef FEM_MESH_H
#define FEM_MESH_H

#include "MathBasic.h"
#include "ParticleSet.h"
#include "ObjMesh.h"
#include "IO.h"
#include <boost/unordered_map.hpp>
#include <boost/property_tree/ptree.hpp>

PRJ_BEGIN

struct SparseReducedBasis;
struct Cubature;
struct FEMInterp;
class FEMLocalBasis;
class FEMSystem;
class FEMCollision;
template <typename T> struct Node;
template <typename T> class ImplicitFunc;

struct FEMVertex : public Serializable {
    typedef Cold Vec;
    FEMVertex();
    bool read(std::istream& is);
    bool write(std::ostream& os) const;
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new FEMVertex);
    }
    Vec3 _pos0,_pos;
    scalar _mass,_matDist;
    sizeType _index;
    bool _surface;
};
struct FEMCell : public Serializable {
    typedef Cold Vec;
    FEMCell();
    Vec3 operator[](sizeType i) const;
    sizeType operator()(sizeType i) const;
    Vec3 get(sizeType i,const Vec* off) const;
    void addParticle(boost::shared_ptr<FEMInterp> I);
    BBox<scalar> getBB() const;
    bool contain(const Vec3& pos,FEMInterp& I) const;
    void findWeight(const Vec3& pos,Vec4& bary,scalar& sqrDist) const;
    Vec3 getVert(const FEMInterp& p,Vec3* N=NULL) const;
    void buildF(Mat3& F,Mat3& FN,Mat3& FR,Eigen::Matrix<scalarD,9,9>* DRDF=NULL,bool forcePositive=true) const;
    void makePositive();
    void buildD();
    bool read(std::istream& is,IOData* dat);
    bool write(std::ostream& os,IOData* dat) const;
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new FEMCell);
    }
    //data
    boost::shared_ptr<FEMInterp> _pSet;
    boost::shared_ptr<FEMVertex> _v[4];
    Mat3d _d;
    scalar _mass;
    sizeType _index;
};
struct FEMInterp : public Serializable {
    typedef Cold Vec;
    FEMInterp();
    FEMInterp(sizeType id);
    FEMInterp(sizeType id,const Vec4& coef);
    bool read(std::istream& is,IOData* dat);
    bool write(std::ostream& os,IOData* dat) const;
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new FEMInterp);
    }
    boost::shared_ptr<FEMCell> _cell;
    boost::shared_ptr<FEMInterp> _next;
    sizeType _id;
    Vec4 _coef;
};
struct FEMBody : public Serializable {
    friend class FEMMesh;
    typedef Cold Vec;
    FEMBody();
    sizeType nrSV() const;
    sizeType dim() const;
    void applyTrans(const Mat4& R,bool pos=true,bool pos0=false);
    void getFace(vector<std::pair<Vec3i,Vec2i> >* surfaceFaces,
                 vector<std::pair<Vec3i,Vec2i> >* internalFaces) const;
    bool read(std::istream& is,IOData* dat);
    bool write(std::ostream& os,IOData* dat) const;
    bool writeABQ(const std::string& path) const;
    void writeObj(ObjMesh& sMesh) const;
    void writeObj(std::ostream& os) const;
    void writeVTK(VTKWriter<scalar>& os,const Vec* off=NULL,const Vec* off0=NULL,sizeType bid=-1,Vec* customV=NULL,Vec* customC=NULL) const;
    void parityCheck() const;
    virtual FEMBody& operator=(const FEMBody& other);
    virtual void assemble();
    virtual void updateMesh();
    void buildG(vector<Eigen::Triplet<scalarD,sizeType> >& trips) const;
    //getter
    void getPos(Vec& X) const;
    void getDPos(Vec& X) const;
    void getMass(Vec& M,bool surface=false) const;
    void getMass(Eigen::DiagonalMatrix<scalarD,-1>& M,bool surface=false) const;
    void getMass(Eigen::SparseMatrix<scalarD,0,sizeType>& M,bool surface=false) const;
    //setter
    void setPos(const Vec& X);
    void setDPos(const Vec& X);
    void clearPos(char dim=-1);
    //other
    sizeType nrV() const {
        return (sizeType)_vss.size();
    }
    sizeType nrC() const {
        return (sizeType)_css.size();
    }
    FEMVertex& getV(sizeType i) {
        return *(_vss[i]);
    }
    FEMCell& getC(sizeType i) {
        return *(_css[i]);
    }
    ParticleSetN& getPSet() {
        return _pSet;
    }
    boost::shared_ptr<FEMVertex> getVPtr(sizeType i) const {
        return _vss[i];
    }
    boost::shared_ptr<FEMCell> getCPtr(sizeType i) const {
        return _css[i];
    }
    const FEMVertex& getV(sizeType i) const {
        return *(_vss[i]);
    }
    const FEMCell& getC(sizeType i) const {
        return *(_css[i]);
    }
    FEMInterp getI(sizeType i) const {
        return *(_evss[i]);
    }
    const ParticleSetN& getPSet() const {
        return _pSet;
    }
    //extra data
    boost::shared_ptr<SparseReducedBasis> _basis;
    boost::shared_ptr<Cubature> _cubatureKinetic,_cubaturePotential;
    boost::shared_ptr<FEMSystem> _system;
    boost::property_tree::ptree _tree;
    Vec _S;
    Vec3 _X0,_W;
    sizeType _offset;
protected:
    vector<boost::shared_ptr<FEMVertex> > _vss;
    vector<boost::shared_ptr<FEMCell> > _css;
    vector<boost::shared_ptr<FEMInterp> > _evss;
    ParticleSetN _pSet;
};
class FEMMesh : public HasMagic
{
public:
    typedef Cold Vec;
    //method
    FEMMesh(sizeType dim,boost::shared_ptr<FEMCollision> coll);
    FEMMesh(const FEMMesh& other);
    void updateMesh(scalar expand=1E-3f);
    void reset(const std::string& path,scalar rad);
    void reset(BBox<scalar> bb,const ImplicitFunc<scalar>& f,scalar rad,scalar cellSz);
    bool reset(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& nodes,
               const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& tets,const ParticleSetN& pset);
    void reset(scalar cellSz,const ParticleSetN& pset,bool usePSet=true);
    void applyTrans(const Mat4& R,sizeType bid,bool pos=true,bool pos0=false);
    FEMCollision& getColl();
    void buildOffset();
    sizeType nrV() const;
    sizeType nrB() const {
        return (sizeType)_bss.size();
    }
    sizeType dim() const {
        return _dim;
    }
    void setCellSz(scalar cellSz) {
        _cellSz=cellSz;
    }
    //body operation
    FEMMesh& operator=(const FEMMesh& other);
    FEMMesh& operator+=(const FEMMesh& other);
    FEMMesh operator+(const FEMMesh& other) const;
    FEMMesh& operator-=(sizeType bid);
    FEMMesh operator-(sizeType bid) const;
    //get bary
    void getBary(const vector<Vec3,Eigen::aligned_allocator<Vec3> >& ps,std::vector<FEMInterp>& evss,scalar sqrDistThres) const;
    FEMBody& getB(sizeType i) {
        return *(_bss[i]);
    }
    boost::shared_ptr<FEMBody> getBPtr(sizeType i) {
        return _bss[i];
    }
    const FEMBody& getB(sizeType i) const {
        return *(_bss[i]);
    }
    //IO
    void writeVTK(const std::string& path) const;
    void writePSetVTK(const std::string& path) const;
    bool write(std::ostream& os) const;
    bool read(std::istream& is);
private:
    void buildBody();
    void buildFace();
    void buildEvss();
    void buildIndex();
    bool assemble();
    //mesh data
    ParticleSetN _pSet;	//just transient data
    vector<boost::shared_ptr<FEMVertex> > _vss;	//just transient data
    vector<boost::shared_ptr<FEMCell> > _css;	//just transient data
    vector<boost::shared_ptr<FEMBody> > _bss;
    boost::shared_ptr<FEMCollision> _coll;
    scalar _cellSz;
    sizeType _dim;
};

PRJ_END

#endif
