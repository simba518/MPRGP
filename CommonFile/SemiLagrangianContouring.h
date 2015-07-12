#ifndef SEMI_LAGRANGIAN_CONTOURING_H
#define SEMI_LAGRANGIAN_CONTOURING_H

#include "Config.h"
#include "MathBasic.h"
#include "CollisionDetection.h"
#include "AABBvh.h"

#include <vector>
#include <stack>
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

#define USE_DOUBLE_FOR_SEMI_LAGRANGIAN_CONTOURING
#ifdef USE_DOUBLE_FOR_SEMI_LAGRANGIAN_CONTOURING
typedef scalarD scalarSLC;
typedef Mat2d Mat2SLC;
typedef Mat3d Mat3SLC;
typedef Vec2d Vec2SLC;
typedef Vec3d Vec3SLC;
#else
typedef scalarF scalarSLC;
typedef Mat2f Mat2SLC;
typedef Mat3f Mat3SLC;
typedef Vec2f Vec2SLC;
typedef Vec3f Vec3SLC;
#endif

using namespace std;

//--------------------------------------------------------------------------------------------------------minimal mesh lib
struct MeshSLC {
    MeshSLC();
    virtual ~MeshSLC() {}
    MeshSLC(vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts,
            vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors,
            vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors,
            vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,
            vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& neighs,
            bool parityCheck=false);
    MeshSLC(vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts,
            vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors,
            vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors,
            vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,
            bool parityCheck=false);
    void parityCheck() const;
    bool search(const sizeType& i,const sizeType& v1,const sizeType& v2) const;
    void searchNeigh();
    void read(istream& is);
    void write(ostream& os) const;
    void writeVTK(const string& path) const;
    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> > _verts;
    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> > _nors;
    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> > _tnors;
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _inds;
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _neighs;
};

//--------------------------------------------------------------------------------------------------------interface to physics simulator
class VelCalc
{
public:
    VelCalc():_dt(0.0f) {}
    virtual ~VelCalc() {}
    virtual Vec3SLC trace(const Vec3SLC& pos) const {
        return pos-vel(pos-vel(pos)*_dt*0.5f)*_dt;
    }
    virtual Vec3SLC vel(const Vec3SLC& pos) const {
        return Vec3SLC::Zero();
    }
    virtual void init(MeshSLC& mesh) const {}
    virtual void setDt(const scalarSLC& dt) {
        _dt=dt;
    }
protected:
    scalarSLC _dt;
};

//--------------------------------------------------------------------------------------------------------AABB bounding box hierarchy
class AABBvh : public AABBvhMeshTpl<scalarSLC>
{
public:
    AABBvh(const MeshSLC* mesh):AABBvhMeshTpl(mesh->_verts,mesh->_inds),_mesh(mesh){}
    const MeshSLC* mesh() const {
        return _mesh;
    }
protected:
    const MeshSLC* _mesh;
};

//--------------------------------------------------------------------------------------------------------octree
class Octree
{
public:
    struct OctNode : public Serializable
    {
        OctNode():Serializable(-1){}
        virtual bool read(istream &is)
        {
            typedef Eigen::Matrix<sizeType,-1,1> Veci;
            readBinaryData(_id,is);
            readBinaryData(_bb,is);
            readBinaryData(_level,is);
            Veci childs;
            readBinaryData(childs,is);
            Eigen::Map<Veci>(_childs,8)=childs;
            readBinaryData(_triIds,is);
            readBinaryData(_flag,is);
            readBinaryData(_valueOff,is);
            return is.good();
        }
        virtual bool write(ostream &os) const
        {
            typedef Eigen::Matrix<sizeType,-1,1> Veci;
            writeBinaryData(_id,os);
            writeBinaryData(_bb,os);
            writeBinaryData(_level,os);
            Veci childs=Eigen::Map<const Veci>(_childs,8);
            writeBinaryData(childs,os);
            writeBinaryData(_triIds,os);
            writeBinaryData(_flag,os);
            writeBinaryData(_valueOff,os);
            return os.good();
        }
        Vec3i _id;
        BBox<scalarSLC> _bb;
        sizeType _level;
        sizeType _childs[8];			//_child[0] >= 0 for internal nodes, _child[0] < 0 for leaves
        Vec2i _triIds;
        unsigned char _flag;
        sizeType _valueOff;				//used by octreeDistance
    };
    struct OctreeProp {
        OctreeProp(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level);
        FORCE_INLINE uint64_t dilate1By2( uint64_t x ) {
            x &= 0x0000ffffLL;
            x = (x ^ (x << 16)) & 0x0000ff0000ffLL;
            x = (x ^ (x << 8))  & 0x00f00f00f00fLL;
            x = (x ^ (x << 4))  & 0x0c30c30c30c3LL;
            x = (x ^ (x << 2))  & 0x249249249249LL;
            return x;
        }
        FORCE_INLINE uint64_t reduce1By2( uint64_t x ) {
            x &= 0x249249249249LL;
            x = (x ^ (x >>  2)) & 0x0c30c30c30c3LL;
            x = (x ^ (x >>  4)) & 0x00f00f00f00fLL;
            x = (x ^ (x >>  8)) & 0x0000ff0000ffLL;
            x = (x ^ (x >> 16)) & 0x0000ffffLL;
            return x;
        }
        FORCE_INLINE uint64_t encodeMorton48( uint64_t x, uint64_t y, uint64_t z ) {
            return dilate1By2(x) | (dilate1By2(y) << 1) | (dilate1By2(z) << 2);
        }
        FORCE_INLINE void decodeMorton48( uint64_t key, uint64_t& x, uint64_t& y, uint64_t& z ) {
            x = reduce1By2( key );
            y = reduce1By2( key >> 1 );
            z = reduce1By2( key >> 2 );
        }
        sizeType _level;
        BBox<scalarSLC> _bb;
        //per block
        Vec3SLC _cellSz;
        Vec3SLC _invCellSz;
        sizeType _nrCellDimPerBlock;
        Vec3i _cellStridePerBlock;
        sizeType _nrCellPerBlock;
        //cross block
        Vec3SLC _blockSz;
        Vec3i _nrBlockDim;
        Vec3i _blockStride;
        sizeType _nrBlock;
        //whole grid
        Vec3SLC _gridSz;
        Vec3i _nrCellDim;
        //whole vertex grid expanded
        Vec3i _nrCellDimVert;
        Vec3i _gridStrideVert;
        sizeType _nrCellVert;
        sizeType _toOriginVert;
        //neigh index LUT
        vector<sizeType,Eigen::aligned_allocator<sizeType> > _childIndexVert;
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _childOffset;
    };
    struct TopDownBuilder {
        virtual void decideSplit(const sizeType& i) =0;
        virtual void insertSplit(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff) =0;
        virtual void postProcessTopDown() =0;
    };
    struct ParityChecker {
        virtual void checkInternal(OctNode& node) =0;
        virtual void checkLeaf(OctNode& node) =0;
    };
    struct TopDownBuilderTriMesh : public Octree::TopDownBuilder {
        TopDownBuilderTriMesh(Octree& tree):_tree(tree) {}
        virtual void decideSplit(const sizeType& i) {
            _tree.decideSplitTriMesh(i);
        }
        virtual void insertSplit(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff) {
            _tree.insertSplitTriMesh(i,ext,coff);
        }
        virtual void postProcessTopDown() {
            _tree.postProcessTopDownTriMesh();
        }
        Octree& _tree;
    };
    struct ParityCheckerTriMesh : public Octree::ParityChecker {
        ParityCheckerTriMesh(Octree& tree):_tree(tree) {}
        virtual void checkInternal(OctNode& node) {
            _tree.parityCheckEveryInternalTriMesh(node);
        }
        virtual void checkLeaf(OctNode& node) {
            _tree.parityCheckEveryLeafTriMesh(node);
        }
        Octree& _tree;
    };
public:
    Octree(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level);
    virtual ~Octree() {}
    virtual void read(istream& is);
    virtual void write(ostream& os) const;
    virtual void writeVTK(const std::string& str) const;
    virtual bool buildFromAABBvhTopDown(const AABBvh* bvh,bool parityCheck);
    const OctreeProp& prop() const {
        return *_prop;
    }
    const vector<sizeType,Eigen::aligned_allocator<sizeType> >& nodePerLevelBeg() const {
        return _nodePerLevelBeg;
    }
    const vector<sizeType,Eigen::aligned_allocator<sizeType> >& nodePerLevelEnd() const {
        return _nodePerLevelEnd;
    }
    const vector<sizeType,Eigen::aligned_allocator<sizeType> >& roots() const {
        return _roots;
    }
    const vector<OctNode,Eigen::aligned_allocator<OctNode> >& nodes() const {
        return _octNodes;
    }
    const MeshSLC* mesh() const {
        return _mesh;
    }
    const AABBvh* bvh() const {
        return _bvh;
    }
    virtual void reset();
    void resetBB(const BBox<scalarSLC>& bb);
public:	//utility functions
    FORCE_INLINE sizeType getChild(const sizeType& i,const sizeType& cid) const {
        return _octNodes[i]._childs[cid];
    }
    FORCE_INLINE sizeType findLeave(const Vec3SLC& pos) const {
        //early out
        if(!_prop->_bb.contain(pos))
            return -1;
        //which finest cell
        const uint64_t x=min<uint64_t>((uint64_t)std::floor((scalarSLC)((pos.x()-_prop->_bb._minC.x())*_prop->_invCellSz.x())),_prop->_nrCellDim.x()-1);
        const uint64_t y=min<uint64_t>((uint64_t)std::floor((scalarSLC)((pos.y()-_prop->_bb._minC.y())*_prop->_invCellSz.y())),_prop->_nrCellDim.y()-1);
        const uint64_t z=min<uint64_t>((uint64_t)std::floor((scalarSLC)((pos.z()-_prop->_bb._minC.z())*_prop->_invCellSz.z())),_prop->_nrCellDim.z()-1);
        const uint64_t morton=_prop->encodeMorton48(x,y,z);
        const uint64_t off=(x>>_prop->_level)*_prop->_blockStride.x()+
                           (y>>_prop->_level)*_prop->_blockStride.y()+
                           (z>>_prop->_level)*_prop->_blockStride.z();
        sizeType root=_roots[off];
        sizeType bitsOffset=(_prop->_level-1)*3;
        while(_octNodes[root]._childs[0] >= 0) {
            root=_octNodes[root]._childs[(morton>>bitsOffset)&7];
            bitsOffset-=3;
        }
        return root;
    }
    FORCE_INLINE sizeType findLeave(const Vec3i& pos) const {
        if(pos.x() < 0 || pos.y() < 0 || pos.z() < 0 ||
                pos.x() >= _prop->_nrCellDim.x() ||
                pos.y() >= _prop->_nrCellDim.y() ||
                pos.z() >= _prop->_nrCellDim.z()) {
            return -1;
        }
        //which finest cell
        const uint64_t morton=_prop->encodeMorton48(pos.x(),pos.y(),pos.z());
        const uint64_t off=(pos.x()>>_prop->_level)*_prop->_blockStride.x()+
                           (pos.y()>>_prop->_level)*_prop->_blockStride.y()+
                           (pos.z()>>_prop->_level)*_prop->_blockStride.z();
        sizeType root=_roots[off];
        sizeType bitsOffset=(_prop->_level-1)*3;
        while(_octNodes[root]._childs[0] >= 0) {
            root=_octNodes[root]._childs[(morton>>bitsOffset)&7];
            bitsOffset-=3;
        }
        return root;
    }
    FORCE_INLINE sizeType findLeave(const Vec3i& pos,const sizeType& l) const {
        if(pos.x() < 0 || pos.y() < 0 || pos.z() < 0 ||
                pos.x() >= _prop->_nrCellDim.x() ||
                pos.y() >= _prop->_nrCellDim.y() ||
                pos.z() >= _prop->_nrCellDim.z()) {
            return -1;
        }
        //which finest cell
        const uint64_t morton=_prop->encodeMorton48(pos.x(),pos.y(),pos.z());
        const uint64_t off=(pos.x()>>_prop->_level)*_prop->_blockStride.x()+
                           (pos.y()>>_prop->_level)*_prop->_blockStride.y()+
                           (pos.z()>>_prop->_level)*_prop->_blockStride.z();
        sizeType root=_roots[off];
        sizeType bitsOffset=(_prop->_level-1)*3;
        while(_octNodes[root]._childs[0] >= 0 && _octNodes[root]._level > l) {
            root=_octNodes[root]._childs[(morton>>bitsOffset)&7];
            bitsOffset-=3;
        }
        return _octNodes[root]._level > l ? -1 : root;
    }
    FORCE_INLINE Vec3i decodeIndex(const sizeType& i,const Vec3i& s) const {
        return Vec3i(i/s.x(),(i%s.x())/s.y(),i%s.y());
    }
    FORCE_INLINE BBox<scalarSLC> getBB(const Vec3i& idMin,const Vec3SLC& len,const Vec3SLC& minC) const {
        BBox<scalarSLC> bb;
        bb._minC=minC+Vec3SLC(idMin.x()*_prop->_cellSz.x(),idMin.y()*_prop->_cellSz.y(),idMin.z()*_prop->_cellSz.z());
        bb._maxC=bb._minC+len;
        return bb;
    }
    FORCE_INLINE Vec3SLC getPt(const Vec3i& id) const {
        return Vec3SLC(id.x()*_prop->_cellSz.x(),
                       id.y()*_prop->_cellSz.y(),
                       id.z()*_prop->_cellSz.z())+_prop->_bb._minC;
    }
protected:	//build from top down
    void buildRoots();
    void buildTreeTopDown(TopDownBuilder& builder);
    void decideSplitTriMesh(const sizeType& i);
    void insertSplitTriMesh(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff);
    void postProcessTopDownTriMesh();
    void buildTriIds(const sizeType& i);
    void parityCheck(ParityChecker& checker);
    void parityCheckEveryInternalTriMesh(OctNode& node);
    void parityCheckEveryLeafTriMesh(OctNode& node);
protected:
    //tree prperty
    boost::shared_ptr<OctreeProp> _prop;
    //tree data
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _triIds;
    vector<OctNode,Eigen::aligned_allocator<OctNode> > _octNodes;
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _roots;
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _nodePerLevelBeg;
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _nodePerLevelEnd;
    //mesh
    const MeshSLC* _mesh;
    const AABBvh* _bvh;
    const Octree* _another;
};

//--------------------------------------------------------------------------------------------------------oct distance tree
enum NEIGH_INDEX {
    NX=0,
    PX,
    NY,
    PY,
    NZ,
    PZ,
    SELF=-1,
};
class OctreeDistance : public Octree
{
public:
    struct TreeBuildHelper {
        //helper data mem
        vector<sizeType,Eigen::aligned_allocator<sizeType> > _ownerNodes;
        vector<sizeType,Eigen::aligned_allocator<sizeType> > _heap;
        vector<unsigned char,Eigen::aligned_allocator<unsigned char> > _ownerVertIds;
        vector<sizeType,Eigen::aligned_allocator<sizeType> > _neighOffs;
        vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> > _neighDxs;
    };
    struct OctreeDistanceProp {
        OctreeDistanceProp(const OctreeProp& prop,TreeBuildHelper& helper);
        //params
        scalarSLC _boxEps;
        sizeType _expInit;
        scalarSLC _edgeCoef;
        //child index
        sizeType _nrNeighBeforeOwner[8];
        Vec3i _neighOffBeforeOwner[8];
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _neighOffBeforeOwnerByLevel;
        //helper data mem
        vector<sizeType,Eigen::aligned_allocator<sizeType> >& _ownerNodes;
        vector<unsigned char,Eigen::aligned_allocator<unsigned char> >& _ownerVertIds;
        vector<sizeType,Eigen::aligned_allocator<sizeType> >& _neighOffs;
        vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& _neighDxs;
    };
    struct VertIter {
        virtual void iter(const Vec3i& id,const sizeType& valueOff) {}
    };
    struct ClearSign : public VertIter {
        ClearSign(OctreeDistance& tree):_tree(tree) {}
        virtual void iter(const Vec3i& id,const sizeType& valueOff) {
            _tree._dProp->_ownerVertIds[valueOff]=0;
        }
        OctreeDistance& _tree;
    };
    struct TopDownBuilderDistance : public OctreeDistance::TopDownBuilder {
        TopDownBuilderDistance(OctreeDistance& tree):_tree(tree) {}
        virtual void decideSplit(const sizeType& i) {
            _tree.decideSplitDistance(i);
        }
        virtual void insertSplit(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff) {
            _tree.insertSplitDistance(i,ext,coff);
        }
        virtual void postProcessTopDown() {}
        OctreeDistance& _tree;
    };
    struct ParityCheckerDistance : public OctreeDistance::ParityChecker {
        ParityCheckerDistance(OctreeDistance& tree):_tree(tree) {}
        virtual void checkInternal(OctNode& node) {
            _tree.parityCheckEveryInternalDistance(node);
        }
        virtual void checkLeaf(OctNode& node) {}
        OctreeDistance& _tree;
    };
    struct SetupSign : public VertIter {
        SetupSign(OctreeDistance& tree,const OctreeDistance& another,const VelCalc& vc)
            :_tree(tree),_another(another),_velCalc(vc) {}
        virtual void iter(const Vec3i& id,const sizeType& valueOff);
        OctreeDistance& _tree;
        const OctreeDistance& _another;
        const VelCalc& _velCalc;
    };
    struct VertDistanceSet : public VertIter {
        VertDistanceSet(OctreeDistance& tree,const OctreeDistance& another,const VelCalc& vc)
            :_tree(tree),_another(another),_velCalc(vc) {}
        virtual void iter(const Vec3i& id,const sizeType& valueOff) {
            _tree._values[valueOff]=_another.getDist(_velCalc.trace(_tree.getPt(id)));
        }
        OctreeDistance& _tree;
        const OctreeDistance& _another;
        const VelCalc& _velCalc;
    };
public:
    OctreeDistance(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level,TreeBuildHelper& helper);
    virtual void read(istream& is);
    virtual void write(ostream& os) const;
    virtual bool buildFromAABBvhTopDown(const AABBvh* bvh,bool parityCheck,bool parityCheckFM=false);
    virtual bool buildFromDistanceTreeTopDown(const OctreeDistance* tree,const VelCalc& vc,MeshSLC* marchTet,bool parityCheck);
    const vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& values() const {
        return _values;
    }
    virtual void reset();
public:	//utility functions
    FORCE_INLINE bool checkVertOwner(const Vec3i& base,const Vec3i& vert,sizeType& off) const {
        off=findLeave(base);
        if(off < 0)
            return false;
        else {
            const OctNode& node=_octNodes[off];
            const sizeType exp=1<<node._level;
            const Vec3i diff=vert-node._id;

            return (diff.x() == 0 || diff.x() == exp) &&
                   (diff.y() == 0 || diff.y() == exp) &&
                   (diff.z() == 0 || diff.z() == exp);
        }
    }
    FORCE_INLINE bool getVertOffset(const Vec3i& pos,sizeType& vid) const {
        for(sizeType c=0; c<8; c++) {
            const sizeType off=findLeave((Vec3i)(pos+_dProp->_neighOffBeforeOwner[c]));
            if(off >= 0) {
                const OctNode& node=_octNodes[off];
                const sizeType exp=1<<node._level;
                const Vec3i diff=pos-node._id;

                if((diff.x() == 0 || diff.x() == exp) &&
                        (diff.y() == 0 || diff.y() == exp) &&
                        (diff.z() == 0 || diff.z() == exp)) {
                    vid=-node._childs[(diff>>node._level).dot(Vec3i(1,2,4))];
                    return true;
                }
            }
        }
        return false;
    }
    FORCE_INLINE scalarSLC getPhi(const Vec3SLC& pos) const {
        const sizeType off=findLeave(pos);
        if(off < 0)
            return ScalarUtil<scalarSLC>::scalar_max;
        else {
            const OctNode& node=_octNodes[off];
            const Vec3SLC rel=((pos-node._bb._minC).array()/node._bb.getExtent().array()).matrix();
            return interp3D(_values[-node._childs[0]],
                            _values[-node._childs[1]],
                            _values[-node._childs[2]],
                            _values[-node._childs[3]],
                            _values[-node._childs[4]],
                            _values[-node._childs[5]],
                            _values[-node._childs[6]],
                            _values[-node._childs[7]],
                            rel.x(),rel.y(),rel.z());
        }
    }
    FORCE_INLINE scalarSLC getDist(const Vec3SLC& pos) const {
        const sizeType off=findLeave(pos);
        if(off < 0)
            return ScalarUtil<scalarSLC>::scalar_max;
        else {
            const OctNode& node=_octNodes[off];
            if(node._level != 0 || node._triIds.y() <= node._triIds.x()) {
                //no triangle in concentric triple
                const Vec3SLC rel=((pos-node._bb._minC).array()/node._bb.getExtent().array()).matrix();
                return interp3D(_values[-node._childs[0]],
                                _values[-node._childs[1]],
                                _values[-node._childs[2]],
                                _values[-node._childs[3]],
                                _values[-node._childs[4]],
                                _values[-node._childs[5]],
                                _values[-node._childs[6]],
                                _values[-node._childs[7]],
                                rel.x(),rel.y(),rel.z());
            } else {
                const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts=_mesh->_verts;
                const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors=_mesh->_nors;
                const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors=_mesh->_tnors;
                const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds=_mesh->_inds;
                const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& neighs=_mesh->_neighs;

                scalarSLC dist=ScalarUtil<scalarSLC>::scalar_max;
                Vec3SLC cpt;
                for(sizeType i=node._triIds.x(); i<node._triIds.y(); i++) {
                    const sizeType idTri=_triIds[i];
                    const Vec3SLC& a=verts[inds[idTri].x()];
                    const Vec3SLC& b=verts[inds[idTri].y()];
                    const Vec3SLC& c=verts[inds[idTri].z()];
                    const Vec3SLC& tn=tnors[idTri];
                    const Vec3SLC& na=nors[inds[idTri].x()];
                    const Vec3SLC& nb=nors[inds[idTri].y()];
                    const Vec3SLC& nc=nors[inds[idTri].z()];
                    const BBox<scalarSLC> bbTri(compMin(compMin(a,b),c),compMax(compMax(a,b),c));
                    updatePointToTriDist(a,b,c,tn,na,nb,nc,
                                         neighs[idTri].x() >= 0 ? (tn+tnors[neighs[idTri].x()]) : tn,
                                         neighs[idTri].y() >= 0 ? (tn+tnors[neighs[idTri].y()]) : tn,
                                         neighs[idTri].z() >= 0 ? (tn+tnors[neighs[idTri].z()]) : tn,
                                         bbTri,pos,dist,cpt);
                }
                return dist;
            }
        }
    }
    FORCE_INLINE void updatePointToTriDist(const Vec3SLC& a,const Vec3SLC& b,const Vec3SLC& c,const Vec3SLC& n,
                                           const Vec3SLC& na,const Vec3SLC& nb,const Vec3SLC& nc,
                                           const Vec3SLC& ena,const Vec3SLC& enb,const Vec3SLC& enc,
                                           const BBox<scalarSLC>& bbTri,const Vec3SLC& pt,
                                           scalarSLC& currVal,Vec3SLC& cpt) const {
        //mem
        scalarSLC newVal;
        Vec3SLC cp,bary;

        //early out
        if(bbTri.distTo(pt) > std::abs(currVal))
            return;

        //calc dist and bary
        TriangleTpl<scalarSLC>(a,b,c).calcPointDist(pt,newVal,cp,bary);
        newVal=sqrt(std::abs(newVal));

        if(newVal >= std::abs(currVal))
            return;
        else cpt=cp;

        //decide sign by bary
        if(bary.x() > 0.0f && bary.y() > 0.0f && bary.z() > 0.0f) {
            currVal=((pt-a).dot(n) < 0.0f) ? -newVal : newVal;
        } else if(bary.x() == 0.0f) {
            if(bary.y() == 0.0f)
                currVal=((pt-a).dot(nc) < 0.0f) ? -newVal : newVal;
            else if(bary.z() == 0.0f)
                currVal=((pt-a).dot(nb) < 0.0f) ? -newVal : newVal;
            else
                currVal=((pt-a).dot(enb) < 0.0f) ? -newVal : newVal;
        } else if(bary.y() == 0.0f) {
            if(bary.x() == 0.0f)
                currVal=((pt-a).dot(nc) < 0.0f) ? -newVal : newVal;
            else if(bary.z() == 0.0f)
                currVal=((pt-a).dot(na) < 0.0f) ? -newVal : newVal;
            else
                currVal=((pt-a).dot(enc) < 0.0f) ? -newVal : newVal;
        } else {
            if(bary.x() == 0.0f)
                currVal=((pt-a).dot(nb) < 0.0f) ? -newVal : newVal;
            else if(bary.y() == 0.0f)
                currVal=((pt-a).dot(na) < 0.0f) ? -newVal : newVal;
            else
                currVal=((pt-a).dot(ena) < 0.0f) ? -newVal : newVal;
        }
    }
    void checkVisibility(sizeType valueOff,scalarSLC& val) const;
    void checkVisibilityInner(const Vec3SLC& pt,scalarSLC& currVal) const;
    scalarSLC computeDistanceExplicit(sizeType valueOff) const;
    scalarSLC computeDistanceExplicit(const Vec3SLC& cpt) const;
    void computeDistanceExplicitInner(const Vec3SLC& cpt,sizeType root,scalarSLC& currVal) const;
protected:	//for tree building
    void decideSplitDistance(const sizeType& i);
    void insertSplitDistance(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff);
    void buildTriIds();
    void accumulateTriId(const sizeType& i);
    void fillTriId(const sizeType& i);
    void parityCheckEveryInternalDistance(OctNode& node);
    void checkTriIds();
protected:	//find vertex list and neighbor
    void buildNeighInfo(VertIter& iter,bool parityCheck);
    void countOwnedVerts(const sizeType& i);
    void decideOwner(OctNode& node,const sizeType& idV);
    void assignOwnedValues(const sizeType& i,const scalarSLC& initVal);
    void assignNonOwnedValuesClearFlag(const sizeType& i);
    void assignNeigh(const sizeType& i,const Vec3SLC& dx);
    void iterVerts(const sizeType& i,VertIter& iter) const;
    void checkVertexInfo() const;	//parity checker
    void checkNeighInfo() const;
    void checkSameVert(const sizeType& off,const Vec3i& id) const;
    void checkNodeNeigh(const sizeType& i) const;
    void checkVertNeigh(const Vec3i& id,const sizeType& i) const;
    void checkClosestVert(const Vec3i& id,const Vec3i& dir,sizeType& neighOff,sizeType& dist) const;
protected:	//for fast marching
    bool vertexFastMarching(VertIter& vi,bool pCheckNeigh,bool pCheckFM);
    void initKnown();
    bool tagIntersectTri(const Vec3SLC& a,const Vec3SLC& b,const Vec3SLC& c);
    void buildInitialDist(const Vec3SLC& a,const Vec3SLC& b,const Vec3SLC& c,const Vec3SLC& n,
                          const Vec3SLC& na,const Vec3SLC& nb,const Vec3SLC& nc,
                          const Vec3SLC& ena,const Vec3SLC& enb,const Vec3SLC& enc);
protected:	//for polygonize
    void polygonize(MeshSLC& mesh,bool parityCheck);
protected:	//for coarsening
    void coarsen();
    void aggregateTriIds(OctNode& node) const;
    bool decideDelete(const OctNode& node) const;
protected:
    //fast marching property
    boost::shared_ptr<OctreeDistanceProp> _dProp;
    //distance data
    vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> > _values;					//distance value list
    //velocity calculator
    const VelCalc* _velCalc;
};

//--------------------------------------------------------------------------------------------------------duplication free octree polygonizer (marching tetrahedron) algorithm
class Polygonizer
{
public:
    Polygonizer(const VelCalc& vc,OctreeDistance& currentTree,
                vector<Octree::OctNode,Eigen::aligned_allocator<Octree::OctNode> >& nodes,
                vector<sizeType,Eigen::aligned_allocator<sizeType> >& infoCache,
                const OctreeDistance& referenceTree,MeshSLC& mesh);
    virtual ~Polygonizer() {}
    void polygonize(bool parityCheck);
protected:
    void countVertPoly(const sizeType& i,const sizeType& idLeave);
    void genVertPoly(const sizeType& i,const sizeType& idLeave);
    void genNeighInfo(const sizeType& i,const sizeType& idLeave);
    void genTriIds(const sizeType& i,const sizeType& idLeave);
    sizeType findVert(const Octree::OctNode& node,const Vec2i& vid) const;
    sizeType findNeighTri(const Octree::OctNode& node,const Vec3i& neighInfo) const;
    Vec3SLC secantSearch(scalarSLC xN1,scalarSLC xN2,
                         scalarSLC fN1,scalarSLC fN2,
                         const Vec3SLC& x0,const Vec3SLC& x1) const;
    bool checkTriNeigh(const sizeType& v0,const sizeType& v1,const Vec3i& b) const;
    sizeType findFace(const sizeType& a,const sizeType& b) const;
protected:
    const VelCalc& _vc;
    OctreeDistance& _currentTree;
    vector<Octree::OctNode,Eigen::aligned_allocator<Octree::OctNode> >& _nodes;
    vector<sizeType,Eigen::aligned_allocator<sizeType> >& _infoCache;
    const OctreeDistance& _referenceTree;
    MeshSLC& _mesh;
protected:
    //params
    scalarSLC _secantTol;
    sizeType _maxIter;
    sizeType _beg;
    sizeType _end;
    //child index
    Vec3i _neighOff[6];
    Vec2i _idEdge[6];
    Vec4i _tetIndices[6];		//which vertices compose the tet.
    Vec3i _neighIndices[6][4];	//(neighbor,tet,face) triple (six tet, four face each)
    Vec2i _vertIndices[6][6];	//(neighbor,vert)
    sizeType _mcNr[16];			//marching tetrahedron index
    Vec3i _mcIndexEdge[16][2];
    Vec3i _mcIndexFace[16][2];
    Vec4i _mcIndexFaceInv[16];
};

//--------------------------------------------------------------------------------------------------------octree vertex fast marching algorithm
class Redistancing
{
public:
    enum FM_STATE {
        UNKNOWN = 0,
        KNOWN = 1,
        CLOSE = 2,
        POSITIVE = 16,
        NEGATIVE = 32,
    };
    Redistancing(vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& values,
                 vector<unsigned char,Eigen::aligned_allocator<unsigned char> >& states,
                 const vector<sizeType,Eigen::aligned_allocator<sizeType> >& neighOffs,
                 const vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& neighDxs);
    virtual ~Redistancing() {}
    bool fastMarch();
protected:
    void tagClose(const sizeType& i);
    sizeType popHeap();
    void pushHeap(sizeType index);
    void updateHeap(sizeType index);
    bool updateClose(const sizeType& i);
    bool updateKnown(const sizeType& i);
protected:
    vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& _values;
    vector<unsigned char,Eigen::aligned_allocator<unsigned char> >& _states;
    const vector<sizeType,Eigen::aligned_allocator<sizeType> >& _neighOffs;
    const vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& _neighDxs;
    //tmp data
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _heapOffsets;
    vector<sizeType,Eigen::aligned_allocator<sizeType> > _heap;
};

//--------------------------------------------------------------------------------------------------------enright test 3D
class VelCalcEnrightTest : public VelCalc
{
public:
    VelCalcEnrightTest():_reverse(false) {}
    void setReverse(bool rev) {
        _reverse=rev;
    }
    virtual Vec3SLC vel(const Vec3SLC& pos) const;
    virtual void init(MeshSLC& mesh) const;
protected:
    bool _reverse;
};

//--------------------------------------------------------------------------------------------------------surface tracking
class SurfaceTracker
{
public:
    SurfaceTracker(const Vec3SLC& cellSz,const sizeType& level);
    virtual ~SurfaceTracker() {}
    void initFromMesh(const MeshSLC& mesh);
    void setJitter(const BBox<scalarSLC>& jitterRange) {
        _jr=jitterRange;
    }
    void advance(const VelCalc& vel);
    void read(istream& is);
    void write(ostream& os) const;
    void writeVTK(const string& path) const;
    bool valid() const {
        return _valid;
    }
protected:
    //param
    Vec3SLC _cellSz;
    Vec3i _nrCellDim;
    BBox<scalarSLC> _bb;
    BBox<scalarSLC> _jr;
    bool _valid;
    //tree
    boost::shared_ptr<OctreeDistance::TreeBuildHelper> _helper;
    boost::shared_ptr<OctreeDistance> _treeA;
    boost::shared_ptr<OctreeDistance> _treeB;
    boost::shared_ptr<MeshSLC> _meshA;
    boost::shared_ptr<MeshSLC> _meshB;
    boost::shared_ptr<AABBvh> _bvh;
};

PRJ_END

#endif
