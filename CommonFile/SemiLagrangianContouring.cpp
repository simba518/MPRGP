#include "SemiLagrangianContouring.h"
#include "IO.h"
#include <set>
#include <algorithm>
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

//--------------------------------------------------------------------------------------------------------minimal mesh lib
MeshSLC::MeshSLC() {}
MeshSLC::MeshSLC
(vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts,
 vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors,
 vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors,
 vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,
 vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& neighs,bool pCheck)
    :_verts(verts),_nors(nors),_tnors(tnors),_inds(inds),_neighs(neighs)
{
    if(pCheck)
        parityCheck();
}
MeshSLC::MeshSLC
(vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts,
 vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors,
 vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors,
 vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,bool pCheck)
    :_verts(verts),_nors(nors),_tnors(tnors),_inds(inds)
{
    searchNeigh();
    if(pCheck)
        parityCheck();
}
void MeshSLC::parityCheck() const
{
    INFOV("%ludu Vertices, %ludu Triangles!",_verts.size(),_inds.size())
    ASSERT_MSG(_nors.size() == _verts.size(),"Vertex Normal Array Size Must Be Same As Vertex Array Size!")
    ASSERT_MSG(_tnors.size() == _inds.size(),"Triangle Normal Array Size Must Be Same As Triangle Index Size!")
    ASSERT_MSG(_neighs.size() == _inds.size(),"Triangle Neigh Array Size Must Be Same As Triangle Index Size!")

    sizeType nrBoundary=0;
    const sizeType nrT=(sizeType)_inds.size();
    for(sizeType i=0; i<nrT; i++) {
        ASSERT(_inds[i].x() != _inds[i].y() &&
               _inds[i].y() != _inds[i].z() &&
               _inds[i].x() != _inds[i].z())

        if(_neighs[i].x() < 0)
            nrBoundary++;
        else
            ASSERT(search(_neighs[i].x(),_inds[i].x(),_inds[i].y()));

        if(_neighs[i].y() < 0)
            nrBoundary++;
        else
            ASSERT(search(_neighs[i].y(),_inds[i].y(),_inds[i].z()));

        if(_neighs[i].z() < 0)
            nrBoundary++;
        else
            ASSERT(search(_neighs[i].z(),_inds[i].z(),_inds[i].x()));
    }

    INFOV("%ludu Boundary Edges Found!",nrBoundary)
    if(nrBoundary > 0)
        WARNING("Non-Watertight, Cannot Fast March!")
    }
bool MeshSLC::search(const sizeType& i,const sizeType& v1,const sizeType& v2) const
{
    sizeType match=0;
    if(_inds[i].x() == v1 || _inds[i].x() == v2)
        match++;
    if(_inds[i].y() == v1 || _inds[i].y() == v2)
        match++;
    if(_inds[i].z() == v1 || _inds[i].z() == v2)
        match++;
    return match == 2;
}
struct key {
    key() {}
    key(const sizeType& a,const sizeType& b) {
        if(a<b) {
            _a=a;
            _b=b;
        } else {
            _a=b;
            _b=a;
        }
    }
    bool operator()(const key& a,const key& b) const {
        return a._a < b._a || (a._a == b._a && a._b < b._b);
    }
    sizeType _a;
    sizeType _b;
};
void MeshSLC::searchNeigh()
{
    map<key,pair<sizeType,sizeType>,key> map;

    _neighs.assign(_inds.size(),Vec3i::Constant(-1));
    for(sizeType i=0; i<(sizeType)_inds.size(); i++) {
        key k=key(_inds[i].x(),_inds[i].y());
        if(map.find(k) != map.end()) {
            _neighs[i].x()=map[k].first;
            _neighs[map[k].first](map[k].second)=i;
        } else map[k]=pair<sizeType,sizeType>(i,0);

        k=key(_inds[i].y(),_inds[i].z());
        if(map.find(k) != map.end()) {
            _neighs[i].y()=map[k].first;
            _neighs[map[k].first](map[k].second)=i;
        } else map[k]=pair<sizeType,sizeType>(i,1);

        k=key(_inds[i].z(),_inds[i].x());
        if(map.find(k) != map.end()) {
            _neighs[i].z()=map[k].first;
            _neighs[map[k].first](map[k].second)=i;
        } else map[k]=pair<sizeType,sizeType>(i,2);
    }
}
void MeshSLC::read(istream& is)
{
    readVector(_verts,is);
    readVector(_nors,is);
    readVector(_tnors,is);
    readVector(_inds,is);
    readVector(_neighs,is);
}
void MeshSLC::write(ostream& os) const
{
    writeVector(_verts,os);
    writeVector(_nors,os);
    writeVector(_tnors,os);
    writeVector(_inds,os);
    writeVector(_neighs,os);
}
void MeshSLC::writeVTK(const string& path) const
{
    VTKWriter<scalarSLC> writer("mesh",path,true);
    writer.appendPoints(_verts.begin(),_verts.end());
    writer.appendCells(_inds.begin(),_inds.end(),VTKWriter<scalarSLC>::TRIANGLE);
}

//--------------------------------------------------------------------------------------------------------octree
Octree::OctreeProp::OctreeProp(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level)
{
    ASSERT(level > 1 && level < 15)
    ASSERT_MSG(sizeof(sizeType) == 2*sizeof(int),"Octree Maybe Very Big, Please Change to int64_t For SizeType!")

    //param
    _level=level;
    _bb=bb;
    //per block
    _cellSz=cellSz;
    _invCellSz=(Vec3SLC::Ones().array()/_cellSz.array()).matrix();
    _nrCellDimPerBlock=(sizeType)(1<<_level);
    _cellStridePerBlock=Vec3i(_nrCellDimPerBlock*_nrCellDimPerBlock,_nrCellDimPerBlock,1);
    _nrCellPerBlock=_nrCellDimPerBlock*_nrCellDimPerBlock*_nrCellDimPerBlock;
    //cross block
    _blockSz=cellSz*(scalarSLC)_nrCellDimPerBlock;
    _nrBlockDim=ceil((Vec3SLC)(_bb.getExtent().array()/_blockSz.array()).matrix());
    _blockStride=Vec3i(_nrBlockDim.y()*_nrBlockDim.z(),_nrBlockDim.z(),1);
    _nrBlock=_nrBlockDim.prod();
    _bb._maxC=_bb._minC+Vec3SLC(_nrBlockDim.x()*_blockSz.x(),
                                _nrBlockDim.y()*_blockSz.y(),
                                _nrBlockDim.z()*_blockSz.z());
    //whole grid
    _gridSz=(Vec3SLC((scalarSLC)_nrBlockDim.x(),(scalarSLC)_nrBlockDim.y(),(scalarSLC)_nrBlockDim.z()).array()*_blockSz.array()).matrix();
    _nrCellDim=_nrBlockDim*_nrCellDimPerBlock;
    //whole vertex grid expanded
    _nrCellDimVert=_nrCellDim+Vec3i::Constant(4);
    _gridStrideVert=Vec3i(_nrCellDimVert.y()*_nrCellDimVert.z(),_nrCellDimVert.z(),1);
    _nrCellVert=_nrCellDimVert.prod();
    _toOriginVert=Vec3i(1,1,1).dot(_gridStrideVert);
    //neigh index LUT
    _childIndexVert.resize((_level+1)*8);
    _childOffset.resize((_level+1)*8);
    for(sizeType i=0; i<=_level; i++) {
        const sizeType off=1<<i;
        for(sizeType c=0; c<8; c++) {
            _childOffset[i*8+c]=Vec3i((c&1),(c&2)>>1,(c&4)>>2)*off;
            _childIndexVert[i*8+c]=_childOffset[i*8+c].dot(_gridStrideVert);
        }
    }
}
Octree::Octree(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level)
{
    //prop
    _prop.reset(new OctreeProp(bb,cellSz,level));
    //tree data
    _triIds.clear();
    _octNodes.clear();
    _roots.resize(_prop->_nrBlock);
    _nodePerLevelBeg.resize(_prop->_level+1);
    _nodePerLevelEnd.resize(_prop->_level+1);
    //data mem
    _mesh=0;
    _bvh=0;
    _another=0;
}
void Octree::read(istream& is)
{
    //parity check
    {
        BBox<scalarSLC> bb;
        readBinaryData(bb,is);

        //reset prop
        {
            ASSERT((bb.getExtent()-_prop->_bb.getExtent()).norm() < EPS)
            const Vec3SLC cellSz=_prop->_cellSz;
            const sizeType level=_prop->_level;
            _prop.reset(new OctreeProp(bb,cellSz,level));
        }

        Vec3SLC cellSz;
        readBinaryData(cellSz,is);
        ASSERT(cellSz == _prop->_cellSz);

        sizeType level;
        readBinaryData(level,is);
        ASSERT(_prop->_level == level);
    }

    //read data
    {
        readVector(_triIds,is);
        readVector(_octNodes,is);
        readVector(_roots,is);
        readVector(_nodePerLevelBeg,is);
        readVector(_nodePerLevelEnd,is);
    }
}
void Octree::write(ostream& os) const
{
    //parity check
    {
        writeBinaryData(_prop->_bb,os);
        writeBinaryData(_prop->_cellSz,os);
        os.write((char*)&(_prop->_level),sizeof(sizeType));
    }

    //write data
    {
        writeVector(_triIds,os);
        writeVector(_octNodes,os);
        writeVector(_roots,os);
        writeVector(_nodePerLevelBeg,os);
        writeVector(_nodePerLevelEnd,os);
    }
}
void Octree::writeVTK(const std::string& str) const
{
    boost::filesystem::create_directory(str);
    std::vector<sizeType,Eigen::aligned_allocator<sizeType> > curr=_roots;
    for(sizeType i=0;!curr.empty();i++)
    {
        std::ostringstream oss;oss << str << "/OCTREE_Level" << i << ".vtk";
        VTKWriter<scalarSLC> os("OCTREE Level",oss.str(),true);
        std::vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> > hex;

        std::vector<sizeType,Eigen::aligned_allocator<sizeType> > nextLv;
        for(sizeType j=0;j<(sizeType)curr.size();j++)
        {
            hex.push_back(_octNodes[curr[j]]._bb._minC);
            hex.push_back(_octNodes[curr[j]]._bb._maxC);
            if(_octNodes[curr[j]]._childs[0] >= 0) {
                nextLv.push_back(_octNodes[curr[j]]._childs[0]);
                nextLv.push_back(_octNodes[curr[j]]._childs[1]);
                nextLv.push_back(_octNodes[curr[j]]._childs[2]);
                nextLv.push_back(_octNodes[curr[j]]._childs[3]);
                nextLv.push_back(_octNodes[curr[j]]._childs[4]);
                nextLv.push_back(_octNodes[curr[j]]._childs[5]);
                nextLv.push_back(_octNodes[curr[j]]._childs[6]);
                nextLv.push_back(_octNodes[curr[j]]._childs[7]);
            }
        }
        os.appendVoxels(hex.begin(),hex.end(),true);
        curr.swap(nextLv);
    }
}
bool Octree::buildFromAABBvhTopDown(const AABBvh* bvh,bool pCheck)
{
    //setup data
    _bvh=bvh;
    _mesh=_bvh->mesh();

    TopDownBuilderTriMesh builder(*this);
    buildTreeTopDown(builder);

    if(pCheck) {
        ParityCheckerTriMesh checker(*this);
        parityCheck(checker);
    }

    return true;
}
void Octree::reset()
{
    _triIds.clear();
    _octNodes.clear();
    _roots.assign(_roots.size(),-1);
    _nodePerLevelBeg.assign(_nodePerLevelBeg.size(),-1);
    _nodePerLevelEnd.assign(_nodePerLevelEnd.size(),-1);
    buildRoots();
    //data mem
    _mesh=0;
    _bvh=0;
    _another=0;
}
void Octree::resetBB(const BBox<scalarSLC>& bb)
{
    const Vec3SLC cellSz=_prop->_cellSz;
    const sizeType level=_prop->_level;
    _prop.reset(new OctreeProp(bb,cellSz,level));
    reset();
}
void Octree::buildRoots()
{
    _octNodes.resize(_prop->_nrBlock);	//reserve space for the root
    {
        const Vec3SLC bbExt=_prop->_cellSz*(1 << _prop->_level);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<_prop->_nrBlock; i++) {
            OctNode& node=_octNodes[i];
            node._id=decodeIndex(i,_prop->_blockStride)*_prop->_nrCellDimPerBlock;
            node._bb=getBB(node._id,bbExt,_prop->_bb._minC);
            node._level=_prop->_level;
            node._childs[0]=-1;
            node._triIds.x()=0;	//root in AABBvh to start from
            node._flag=false;	//tag don't split
            _roots[i]=i;
        }
        _nodePerLevelBeg[_prop->_level]=0;
        _nodePerLevelEnd[_prop->_level]=(sizeType)_octNodes.size();
    }
}
void Octree::buildTreeTopDown(TopDownBuilder& builder)
{
    //generate root
    buildRoots();

    //build tree by level
    for(sizeType l=_prop->_level-1; l>=0; l--) {
        //build node level start
        _nodePerLevelBeg[l]=(sizeType)_octNodes.size();

        //decide how many nodes to split
        OMP_PARALLEL_FOR_
        for(sizeType i=_nodePerLevelBeg[l+1]; i<_nodePerLevelEnd[l+1]; i++)
            builder.decideSplit(i);

        //scan the node for split
        sizeType nrNodeLevel=(sizeType)_octNodes.size();
        for(sizeType i=_nodePerLevelBeg[l+1]; i<_nodePerLevelEnd[l+1]; i++) {
            if(_octNodes[i]._flag) {
                _octNodes[i]._triIds.y()=nrNodeLevel;
                nrNodeLevel+=8;
            }
        }
        _octNodes.resize(nrNodeLevel);

        //insert nodes
        const Vec3i* coff=&(_prop->_childOffset[l*8]);
        const Vec3SLC ext=_prop->_cellSz*(1<<l);
        OMP_PARALLEL_FOR_
        for(sizeType i=_nodePerLevelBeg[l+1]; i<_nodePerLevelEnd[l+1]; i++)
            builder.insertSplit(i,ext,coff);

        //build level node end index
        _nodePerLevelEnd[l]=nrNodeLevel;
    }

    builder.postProcessTopDown();
}
void Octree::decideSplitTriMesh(const sizeType& i)
{
    OctNode& node=_octNodes[i];
    const Vec3SLC ext=node._bb.getExtent();
    const BBox<scalarSLC> bbExp(compMax(node._bb._minC-ext,_prop->_bb._minC),
                                compMin(node._bb._maxC+ext,_prop->_bb._maxC));

    if(node._triIds.x() >= 0 && _bvh->intersect(bbExp,node._triIds.x()))
        node._flag=true;
    else
        node._childs[0]=-1;
}
void Octree::insertSplitTriMesh(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff)
{
    OctNode& node=_octNodes[i];
    if(!node._flag)
        return;

    for(sizeType c=0; c<8; c++) {
        node._childs[c]=node._triIds.y()+c;
        OctNode& child=_octNodes[node._childs[c]];

        child._id=node._id+coff[c];
        child._bb=getBB(child._id,ext,_prop->_bb._minC);
        child._level=node._level-1;
        child._flag=false;
        child._childs[0]=-1;

        //descend root
        const BBox<scalarSLC> expBB(compMax(child._bb._minC-ext,_prop->_bb._minC),
                                    compMin(child._bb._maxC+ext,_prop->_bb._maxC));
        child._triIds.x()=node._triIds.x();
        _bvh->decideRoot(expBB,child._triIds.x());

        //find nrTri intersect
        if(child._level == 0) {
            child._childs[0]=-1;
            if(child._triIds.x() >= 0)
                child._triIds.y()=_bvh->intersectNr(expBB,child._triIds.x());
            else
                child._triIds.y()=0;
        }
    }
}
void Octree::postProcessTopDownTriMesh()
{
    //count number of triIds needed
    sizeType nrTriIds=0;
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++) {
        _octNodes[i]._childs[1]=nrTriIds;
        nrTriIds+=_octNodes[i]._triIds.y();
    }
    _triIds.resize(nrTriIds);

    //build triIds
    OMP_PARALLEL_FOR_
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++)
        buildTriIds(i);
}
void Octree::buildTriIds(const sizeType& i)
{
    typedef vector<sizeType,Eigen::aligned_allocator<sizeType> >::iterator ITER;

    OctNode& node=_octNodes[i];
    const BBox<scalarSLC> expBB(compMax(node._bb._minC-_prop->_cellSz,_prop->_bb._minC),
                                compMin(node._bb._maxC+_prop->_cellSz,_prop->_bb._maxC));

    if(node._triIds.y() > 0) {
        const sizeType from=node._childs[1];
        ITER iter=_triIds.begin()+node._childs[1];
        ASSERT(_bvh->intersect<ITER>(expBB,node._triIds.x(),iter));
        node._triIds=Vec2i(from,from+node._triIds.y());
    } else {
        node._triIds=Vec2i::Constant(-1);
    }
}
void Octree::parityCheck(ParityChecker& checker)
{
    //check integrity of tree
    const sizeType nrNode=(sizeType)_octNodes.size();
    for(sizeType i=0; i<(sizeType)nrNode; i++)
        _octNodes[i]._flag=0;

    stack<sizeType> s;
    for(sizeType i=0; i<_prop->_nrBlock; i++) {
        s.push(_roots[i]);
        ASSERT(_octNodes[_roots[i]]._level == _prop->_level);
        ASSERT((_octNodes[_roots[i]]._bb.getExtent()-_prop->_cellSz*(1 << _prop->_level)).norm() < EPS)
    }

    while(!s.empty()) {
        const sizeType pos=s.top();
        s.pop();
        OctNode& node=_octNodes[pos];
        ASSERT(node._flag == 0);
        node._flag=1;

        BBox<scalarSLC> bb=getBB(node._id,_prop->_cellSz*(1 << node._level),_prop->_bb._minC);
        ASSERT((node._bb._minC-bb._minC).norm() < EPS);
        ASSERT((node._bb._maxC-bb._maxC).norm() < EPS);

        if(node._childs[0] >= 0) {
            for(sizeType c=0; c<8; c++) {
                const OctNode& child=_octNodes[node._childs[c]];
                ASSERT(node._level-1 == child._level);
                ASSERT(node._bb.enlarge(EPS).contain(child._bb))
                ASSERT(node._id+_prop->_childOffset[child._level*8+c] == child._id);
                s.push(node._childs[c]);
            }

            checker.checkInternal(node);
        } else {
            sizeType found=findLeave((Vec3SLC)((node._bb._minC+node._bb._maxC)*0.5f));
            ASSERT(found == pos);
            found=findLeave(node._id);
            ASSERT(found == pos);

            checker.checkLeaf(node);
        }
    }

    for(sizeType i=0; i<(sizeType)nrNode; i++)
        ASSERT(_octNodes[i]._flag == 1);
}
void Octree::parityCheckEveryInternalTriMesh(OctNode& node)
{
    //ready mesh
    const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts=_mesh->_verts;
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds=_mesh->_inds;
    const sizeType nrT=(sizeType)inds.size();

    BBox<scalarSLC> expBB;
    expBB._minC=node._bb._minC-node._bb.getExtent();
    expBB._maxC=node._bb._maxC+node._bb.getExtent();
    expBB._minC=compMax(expBB._minC,_prop->_bb._minC);
    expBB._maxC=compMin(expBB._maxC,_prop->_bb._maxC);

    if(_bvh) {
        ASSERT(_bvh->intersect(expBB,0))
    } else {
        bool inter=false;
        for(sizeType i=0; i<nrT; i++) {
            if(TriangleTpl<scalarSLC>(verts[inds[i].x()],verts[inds[i].y()],verts[inds[i].z()]).intersect(expBB)) {
                inter=true;
                break;
            }
        }
        ASSERT(inter);
    }
}
void Octree::parityCheckEveryLeafTriMesh(OctNode& node)
{
    //ready mesh
    const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts=_mesh->_verts;
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds=_mesh->_inds;
    const sizeType nrT=(sizeType)inds.size();

    BBox<scalarSLC> expBB;
    expBB._minC=node._bb._minC-node._bb.getExtent();
    expBB._maxC=node._bb._maxC+node._bb.getExtent();
    expBB._minC=compMax(expBB._minC,_prop->_bb._minC);
    expBB._maxC=compMin(expBB._maxC,_prop->_bb._maxC);

    if(node._level != 0) {
        if(_bvh) {
            ASSERT(!_bvh->intersect(expBB,0))
        } else {
            for(sizeType i=0; i<nrT; i++)
                ASSERT(!TriangleTpl<scalarSLC>(verts[inds[i].x()],verts[inds[i].y()],verts[inds[i].z()]).intersect(expBB));
        }
    } else if(_bvh) {
        set<sizeType> total;
        for(sizeType it=node._triIds.x(); it < node._triIds.y(); it++)
            total.insert(_triIds[it]);

        set<sizeType> totalExp;
        std::insert_iterator<set<sizeType> > inserter(totalExp,totalExp.end());
        _bvh->intersect<std::insert_iterator<set<sizeType> > >(expBB,0,inserter);

        ASSERT(totalExp == total);
    }
}

//--------------------------------------------------------------------------------------------------------oct distance tree
OctreeDistance::OctreeDistanceProp::OctreeDistanceProp(const OctreeProp& prop,TreeBuildHelper& helper)
    :_ownerNodes(helper._ownerNodes),
     _ownerVertIds(helper._ownerVertIds),
     _neighOffs(helper._neighOffs),
     _neighDxs(helper._neighDxs)
{
    ASSERT_MSG(sizeof(sizeType) == 2*sizeof(int),"Performing Octree Fast Marching Requires int64_t for SizeType!")

    _boxEps=1E-3f;
    _expInit=3;
    _edgeCoef=3.0f;

    _nrNeighBeforeOwner[0]=0;
    _nrNeighBeforeOwner[1]=4;
    _nrNeighBeforeOwner[2]=2;
    _nrNeighBeforeOwner[3]=6;
    _nrNeighBeforeOwner[4]=1;
    _nrNeighBeforeOwner[5]=5;
    _nrNeighBeforeOwner[6]=3;
    _nrNeighBeforeOwner[7]=7;

    _neighOffBeforeOwner[0]=Vec3i( 0, 0, 0);
    _neighOffBeforeOwner[1]=Vec3i( 0, 0,-1);
    _neighOffBeforeOwner[2]=Vec3i( 0,-1, 0);
    _neighOffBeforeOwner[3]=Vec3i( 0,-1,-1);
    _neighOffBeforeOwner[4]=Vec3i(-1, 0, 0);
    _neighOffBeforeOwner[5]=Vec3i(-1, 0,-1);
    _neighOffBeforeOwner[6]=Vec3i(-1,-1, 0);
    _neighOffBeforeOwner[7]=Vec3i(-1,-1,-1);

    _neighOffBeforeOwnerByLevel.resize((prop._level+1)*64);
    for(sizeType l=0; l<=prop._level; l++) {
        const sizeType off=1<<l;
        for(sizeType i=0; i<8; i++) {
            const Vec3i toNode=Vec3i((i&1),(i&2)>>1,(i&4)>>2)*off;
            for(sizeType j=0; j<8; j++)
                _neighOffBeforeOwnerByLevel[l*64+i*8+j]=
                    _neighOffBeforeOwner[j]+toNode;
        }
    }
}
OctreeDistance::OctreeDistance(const BBox<scalarSLC>& bb,const Vec3SLC& cellSz,const sizeType& level,TreeBuildHelper& helper)
    :Octree(bb,cellSz,level),_dProp(new OctreeDistanceProp(*_prop,helper)) {}
void OctreeDistance::SetupSign::iter(const Vec3i& id,const sizeType& valueOff)
{
    const scalarSLC phi=_another.getPhi(_velCalc.trace(_tree.getPt(id)));
    _tree._dProp->_ownerVertIds[valueOff]=(phi > 0.0f) ? Redistancing::POSITIVE : Redistancing::NEGATIVE;
}
void OctreeDistance::read(istream& is)
{
    Octree::read(is);
    readVector(_values,is);
}
void OctreeDistance::write(ostream& os) const
{
    Octree::write(os);
    writeVector(_values,os);
}
bool OctreeDistance::buildFromAABBvhTopDown(const AABBvh* bvh,bool parityCheck,bool parityCheckFM)
{
    Octree::buildFromAABBvhTopDown(bvh,parityCheck);

    ClearSign signClearer(*this);
    if(!vertexFastMarching(signClearer,true,parityCheckFM)) {
        reset();
        return false;
    }
    return true;
}
bool OctreeDistance::buildFromDistanceTreeTopDown(const OctreeDistance* tree,const VelCalc& vc,MeshSLC* marchTet,bool pCheck)
{
    //setup data
    _another=tree;
    _bvh=0;
    _mesh=marchTet;

    //tracer
    _velCalc=&vc;

    TopDownBuilderDistance builder(*this);
    buildTreeTopDown(builder);

    polygonize(*marchTet,pCheck);
    coarsen();

    //parity check
    if(pCheck) {
        sizeType nrTri=0;
        for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++) {
            nrTri+=_octNodes[i]._triIds.y()-_octNodes[i]._triIds.x();
            if(_octNodes[i]._triIds.y()-_octNodes[i]._triIds.x() > 0)
                ASSERT(_octNodes[i]._flag)
                else
                    ASSERT(!_octNodes[i]._flag)
                }
        ASSERT(nrTri == (sizeType)marchTet->_inds.size())

        ParityCheckerDistance checker(*this);
        OctreeDistance::parityCheck(checker);
    }

    buildTriIds();

    //parity check
    if(pCheck)
        checkTriIds();

    SetupSign signSetter(*this,*tree,vc);
    if(!vertexFastMarching(signSetter,false,pCheck)) {
        reset();
        return false;
    }
    return true;
}
void OctreeDistance::reset()
{
    Octree::reset();

    VertIter defaultIter;
    buildNeighInfo(defaultIter,true);
}
void OctreeDistance::checkVisibility(sizeType valueOff,scalarSLC& val) const
{
    const OctNode& node=_octNodes[_dProp->_ownerNodes[valueOff]];
    for(sizeType i=0;i<8;i++){
        if(node._childs[i] == -valueOff){
            Vec3SLC pt;
            pt[0]=(i&1) ? node._bb._maxC[0] : node._bb._minC[0];
            pt[1]=(i&2) ? node._bb._maxC[1] : node._bb._minC[1];
            pt[2]=(i&4) ? node._bb._maxC[2] : node._bb._minC[2];
            return checkVisibilityInner(pt,val);
        }
    }
}
void OctreeDistance::checkVisibilityInner(const Vec3SLC& pt,scalarSLC& currVal) const
{
    Vec3SLC X=Vec3SLC::Unit(0)*_prop->_bb.getExtent()[0];
    Vec3SLC Y=Vec3SLC::Unit(1)*_prop->_bb.getExtent()[1];
    Vec3SLC Z=Vec3SLC::Unit(2)*_prop->_bb.getExtent()[2];

    //rain check the distance with visibility
    bool visible=false;
    if(!_bvh->intersectLineSeg3D(pt,pt+X,0))
        visible=true;
    else if(!_bvh->intersectLineSeg3D(pt,pt-X,0))
        visible=true;
    else if(!_bvh->intersectLineSeg3D(pt,pt+Y,0))
        visible=true;
    else if(!_bvh->intersectLineSeg3D(pt,pt-Y,0))
        visible=true;
    else if(!_bvh->intersectLineSeg3D(pt,pt+Z,0))
        visible=true;
    else if(!_bvh->intersectLineSeg3D(pt,pt-Z,0))
        visible=true;

    if(!visible && currVal > 0.0f)
        currVal*=-1.0f;
    if(visible && currVal < 0.0f)
        currVal*=-1.0f;
}
scalarSLC OctreeDistance::computeDistanceExplicit(sizeType valueOff) const
{
    const OctNode& node=_octNodes[_dProp->_ownerNodes[valueOff]];
    for(sizeType i=0;i<8;i++){
        if(node._childs[i] == -valueOff){
            Vec3SLC pt;
            pt[0]=(i&1) ? node._bb._maxC[0] : node._bb._minC[0];
            pt[1]=(i&2) ? node._bb._maxC[1] : node._bb._minC[1];
            pt[2]=(i&4) ? node._bb._maxC[2] : node._bb._minC[2];
            scalarSLC ret=computeDistanceExplicit(pt);
            checkVisibilityInner(pt,ret);
            return ret;
        }
    }
    ASSERT(false);
    return 0.0f;
}
scalarSLC OctreeDistance::computeDistanceExplicit(const Vec3SLC& pt) const
{
    scalarSLC ret=ScalarUtil<scalarSLC>::scalar_max;
    sizeType root=0;
    computeDistanceExplicitInner(pt,root,ret);
    return ret;
}
void OctreeDistance::computeDistanceExplicitInner(const Vec3SLC& pt,sizeType root,scalarSLC& currVal) const
{
    Vec3SLC cpt;
    const AABBvh::BVHNode& node=_bvh->getNode(root);
    if(node._bb.distTo(pt) > std::abs(currVal))
        return;
    if(node._tid < 0){
        computeDistanceExplicitInner(pt,_bvh->getLeft(root),currVal);
        computeDistanceExplicitInner(pt,_bvh->getRight(root),currVal);
    }else{
        sizeType i=node._tid;
        Vec3i ind=_mesh->_inds[i];
        Vec3i n=_mesh->_neighs[i];
        updatePointToTriDist
            (_mesh->_verts[ind.x()],//vertex angled weighted normals
             _mesh->_verts[ind.y()],
             _mesh->_verts[ind.z()],
             _mesh->_tnors[i],//triangle normals
             _mesh->_nors[ind.x()],//vertex normals
             _mesh->_nors[ind.y()],
             _mesh->_nors[ind.z()],
             n.x() >= 0 ? (_mesh->_tnors[i]+_mesh->_tnors[n.x()]) : _mesh->_tnors[i],//edge normal direction not normalized (0,1),(1,2),(2,0)
             n.y() >= 0 ? (_mesh->_tnors[i]+_mesh->_tnors[n.y()]) : _mesh->_tnors[i],
             n.z() >= 0 ? (_mesh->_tnors[i]+_mesh->_tnors[n.z()]) : _mesh->_tnors[i],
             node._bb,pt,currVal,cpt);
    }
}
void OctreeDistance::decideSplitDistance(const sizeType& i)
{
    OctNode& node=_octNodes[i];
    const Vec3SLC ctr=(node._bb._minC+node._bb._maxC)*0.5f;
    const scalarSLC edgeLen=_prop->_cellSz.minCoeff()*(scalarSLC)(1<<node._level);
    const scalarSLC phi=static_cast<const OctreeDistance*>(_another)->getPhi(_velCalc->trace(ctr));

    if(std::abs(phi) < _dProp->_edgeCoef*edgeLen)
        node._flag=true;
    else
        node._childs[0]=-1;
}
void OctreeDistance::insertSplitDistance(const sizeType& i,const Vec3SLC& ext,const Vec3i* coff)
{
    OctNode& node=_octNodes[i];
    if(!node._flag)
        return;

    for(sizeType c=0; c<8; c++) {
        node._childs[c]=node._triIds.y()+c;
        OctNode& child=_octNodes[node._childs[c]];

        child._id=node._id+coff[c];
        child._bb=getBB(child._id,ext,_prop->_bb._minC);
        child._level=node._level-1;
        child._flag=false;
        child._childs[0]=-1;
    }
}
void OctreeDistance::buildTriIds()
{
    _triIds.clear();
    _dProp->_ownerNodes.resize(_octNodes.size());

    //count number of triIds needed
    OMP_PARALLEL_FOR_
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++)
        accumulateTriId(i);

    sizeType nrTri=0;
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++) {
        _dProp->_ownerNodes[i]=_octNodes[i]._valueOff;
        _octNodes[i]._valueOff=nrTri;
        nrTri+=_dProp->_ownerNodes[i];
    }
    _triIds.resize(nrTri);

    OMP_PARALLEL_FOR_
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++)
        fillTriId(i);

    OMP_PARALLEL_FOR_
    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++) {
        _octNodes[i]._triIds.x()=_octNodes[i]._valueOff;
        _octNodes[i]._triIds.y()=_octNodes[i]._triIds.x()+_dProp->_ownerNodes[i];
    }
}
void OctreeDistance::accumulateTriId(const sizeType& i)
{
    OctNode& node=_octNodes[i];
    ASSERT(node._level == 0)

    node._valueOff=0;
    for(sizeType x=node._id.x()-1; x<=node._id.x()+1; x++)
        for(sizeType y=node._id.y()-1; y<=node._id.y()+1; y++)
            for(sizeType z=node._id.z()-1; z<=node._id.z()+1; z++) {
                const sizeType off=findLeave(Vec3i(x,y,z));
                if(off >= 0) {
                    const OctNode& neigh=_octNodes[off];
                    if(neigh._level == 0)
                        node._valueOff+=neigh._triIds.y()-neigh._triIds.x();
                }
            }
    if(node._valueOff)
        node._flag=true;
    else
        node._flag=false;
}
void OctreeDistance::fillTriId(const sizeType& i)
{
    OctNode& node=_octNodes[i];
    ASSERT(node._level == 0)

    if(node._flag) {
        sizeType* triIds=&(_triIds[node._valueOff]);
        sizeType index=0;
        for(sizeType x=node._id.x()-1; x<=node._id.x()+1; x++)
            for(sizeType y=node._id.y()-1; y<=node._id.y()+1; y++)
                for(sizeType z=node._id.z()-1; z<=node._id.z()+1; z++) {
                    const sizeType off=findLeave(Vec3i(x,y,z));
                    if(off >= 0) {
                        const OctNode& neigh=_octNodes[off];
                        if(neigh._level == 0)
                            for(sizeType t=neigh._triIds.x(); t<neigh._triIds.y(); t++)
                                triIds[index++]=t;
                    }
                }
    }
}
void OctreeDistance::parityCheckEveryInternalDistance(OctNode& node)
{
    const sizeType ext=1<<node._level;
    const Vec3i minC=node._id-Vec3i::Constant(ext);
    const Vec3i maxC=node._id+Vec3i::Constant(ext)*2;

    unsigned char gen=false;
    for(sizeType x=minC.x(); x<maxC.x(); x++)
        for(sizeType y=minC.y(); y<maxC.y(); y++)
            for(sizeType z=minC.z(); z<maxC.z(); z++) {
                sizeType off=findLeave(Vec3i(x,y,z));
                if(off != -1 && _octNodes[off]._level == 0)
                    if(_octNodes[off]._triIds.y()-_octNodes[off]._triIds.x() > 0)
                        gen=true;
            }
    ASSERT(gen)
}
void OctreeDistance::checkTriIds()
{
    AABBvh bvh(_mesh);

    for(sizeType i=_nodePerLevelBeg[0]; i<_nodePerLevelEnd[0]; i++) {
        const OctNode& node=_octNodes[i];
        const Vec3SLC exp=node._bb.getExtent()-Vec3SLC::Constant(EPS);
        const BBox<scalarSLC> expBB(node._bb._minC-exp,node._bb._maxC+exp);

        set<sizeType> total;
        for(sizeType it=node._triIds.x(); it < node._triIds.y(); it++)
            total.insert(_triIds[it]);

        set<sizeType> totalExp;
        std::insert_iterator<set<sizeType> > inserter(totalExp,totalExp.end());
        bvh.intersect<std::insert_iterator<set<sizeType> > >(expBB,0,inserter);

        ASSERT(totalExp == total);
    }
}
void OctreeDistance::buildNeighInfo(VertIter& iter,bool pCheck)
{
    const sizeType nrNode=(sizeType)_octNodes.size();

    //clear data mem
    _values.clear();
    _dProp->_ownerNodes.clear();
    _dProp->_ownerVertIds.clear();

    //generate duplication free vertex info
    {
        //count how many verts they own
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrNode; i++)
            countOwnedVerts(i);

        //now in every leaf:
        //_childs[1] store the start index from values begin list
        //_childs[2] store the number of owned values
        sizeType nrValues=1;	//we have to offset value by one to ensure leaf._child[0] < 0
        for(sizeType i=0; i<nrNode; i++) {
            OctNode& node=_octNodes[i];
            if(node._childs[0] < 0) {
                node._valueOff=nrValues;
                nrValues+=countBits<sizeType>(node._flag);
            }
        }

        //allocate space
        _values.resize(nrValues);
        _dProp->_ownerNodes.resize(nrValues);
        _dProp->_ownerVertIds.resize(nrValues);
        //allocate neighbor info
        _dProp->_neighOffs.resize(nrValues*6);
        _dProp->_neighDxs.resize(nrValues*6);
        //for safety
        _dProp->_ownerVertIds[0]=-1;

        //assign owned values
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrNode; i++)
            assignOwnedValues(i,ScalarUtil<scalarSLC>::scalar_max);

        //assign non-owned indices
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrNode; i++)
            assignNonOwnedValuesClearFlag(i);
    }
    if(pCheck)
        checkVertexInfo();

    //generate neigh info
    fill(_dProp->_neighOffs.begin(),_dProp->_neighOffs.end(),-1);
    fill(_dProp->_neighDxs.begin(),_dProp->_neighDxs.end(),0.0f);
    {
        //assign info from top down
        for(sizeType l=_prop->_level; l>=0; l--) {
            const Vec3SLC dx=_prop->_cellSz*(1<<l);
            OMP_PARALLEL_FOR_
            for(sizeType i=_nodePerLevelBeg[l]; i<_nodePerLevelEnd[l]; i++)
                assignNeigh(i,dx);
        }
    }
    if(pCheck)
        checkNeighInfo();

    //clear flag
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrNode; i++)
        iterVerts(i,iter);
}
void OctreeDistance::countOwnedVerts(const sizeType& i)
{
    //switch on all owner ship
    OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

    //set self
    node._flag=1<<0;
    node._childs[0]=-(i+1);

    //set all other
    decideOwner(node,1);
    decideOwner(node,2);
    decideOwner(node,3);
    decideOwner(node,4);
    decideOwner(node,5);
    decideOwner(node,6);
    decideOwner(node,7);
}
void OctreeDistance::decideOwner(OctNode& node,const sizeType& idV)
{
    const Vec3i* offs=&(_dProp->_neighOffBeforeOwnerByLevel[node._level*64+idV*8]);
    for(sizeType i=0; i<_dProp->_nrNeighBeforeOwner[idV]; i++) {
        sizeType nid;
        if(checkVertOwner(node._id+offs[i],node._id+offs[0],nid)) {
            node._childs[idV]=-(nid+1);
            return;
        }
    }
    node._flag|=1<<idV;
}
void OctreeDistance::assignOwnedValues(const sizeType& i,const scalarSLC& initVal)
{
    OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

    sizeType offV=node._valueOff;
    for(sizeType v=0; v<8; v++) {
        if(node._flag&(sizeType(1)<<v)) {
            //set node offset
            node._childs[v]=-offV;
            ASSERT(-offV < 0)
            //set value list
            _values[offV]=initVal;
            _dProp->_ownerNodes[offV]=i;
            _dProp->_ownerVertIds[offV++]=(unsigned char)v;
        }
    }
}
void OctreeDistance::assignNonOwnedValuesClearFlag(const sizeType& i)
{
    OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

    for(sizeType v=0; v<8; v++) {
        if(!(node._flag&(sizeType(1)<<v))) {
            const OctNode& neigh=_octNodes[-node._childs[v]-1];
            const Vec3i relId=(node._id-neigh._id)+_prop->_childOffset[node._level*8+v];
            node._childs[v]=neigh._childs[(relId>>neigh._level).dot(Vec3i(1,2,4))];
        }
    }
}
void OctreeDistance::assignNeigh(const sizeType& i,const Vec3SLC& dx)
{
    OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

#define SET_NEIGH(from,to,n,p,d)								\
{_dProp->_neighOffs[-node._childs[to]*6+n]=-node._childs[from];	\
 _dProp->_neighDxs[-node._childs[to]*6+n]=d;}					\
{_dProp->_neighOffs[-node._childs[from]*6+p]=-node._childs[to];	\
 _dProp->_neighDxs[-node._childs[from]*6+p]=d;}

    OMP_CRITICAL_ {
        SET_NEIGH(0,1,NX,PX,dx.x())
        SET_NEIGH(1,3,NY,PY,dx.y())
        SET_NEIGH(2,3,NX,PX,dx.x())
        SET_NEIGH(0,2,NY,PY,dx.y())

        SET_NEIGH(4,5,NX,PX,dx.x())
        SET_NEIGH(5,7,NY,PY,dx.y())
        SET_NEIGH(6,7,NX,PX,dx.x())
        SET_NEIGH(4,6,NY,PY,dx.y())

        SET_NEIGH(0,4,NZ,PZ,dx.z())
        SET_NEIGH(1,5,NZ,PZ,dx.z())
        SET_NEIGH(3,7,NZ,PZ,dx.z())
        SET_NEIGH(2,6,NZ,PZ,dx.z())
    }

#undef SET_NEIGH
}
void OctreeDistance::iterVerts(const sizeType& i,VertIter& iter) const
{
    const OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

    const Vec3i* offs=&(_prop->_childOffset[node._level*8]);
    for(sizeType v=0; v<8; v++) {
        if(node._flag&(sizeType(1)<<v))
            iter.iter(node._id+offs[v],-node._childs[v]);
    }
}
struct LessVec3i {
    bool operator()(const Vec3i& a,const Vec3i& b) const {
        return a.x() < b.x() ||
               (a.x() == b.x() && a.y() < b.y()) ||
               (a.x() == b.x() && a.y() == b.y() && a.z() < b.z());
    }
};
void OctreeDistance::checkVertexInfo() const
{
    const sizeType nrNode=(sizeType)_octNodes.size();

    {
        set<Vec3i,LessVec3i> unique;
        for(sizeType i=0; i<nrNode; i++) {
            const OctNode& node=_octNodes[i];
            if(node._childs[0] < 0) {
                unique.insert(node._id+_prop->_childOffset[node._level*8+0]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+1]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+2]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+3]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+4]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+5]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+6]);
                unique.insert(node._id+_prop->_childOffset[node._level*8+7]);
            }
        }
        ASSERT(_values.size() == unique.size()+1);
    }

    {
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrNode; i++) {
            const OctNode& node=_octNodes[i];
            if(node._childs[0] < 0) {
                checkSameVert(-node._childs[0],node._id+_prop->_childOffset[node._level*8+0]);
                checkSameVert(-node._childs[1],node._id+_prop->_childOffset[node._level*8+1]);
                checkSameVert(-node._childs[2],node._id+_prop->_childOffset[node._level*8+2]);
                checkSameVert(-node._childs[3],node._id+_prop->_childOffset[node._level*8+3]);
                checkSameVert(-node._childs[4],node._id+_prop->_childOffset[node._level*8+4]);
                checkSameVert(-node._childs[5],node._id+_prop->_childOffset[node._level*8+5]);
                checkSameVert(-node._childs[6],node._id+_prop->_childOffset[node._level*8+6]);
                checkSameVert(-node._childs[7],node._id+_prop->_childOffset[node._level*8+7]);
            }
        }
    }
}
void OctreeDistance::checkNeighInfo() const
{
    const sizeType nrNode=(sizeType)_octNodes.size();
    for(sizeType i=0; i<nrNode; i++)
        checkNodeNeigh(i);
}
void OctreeDistance::checkSameVert(const sizeType& off,const Vec3i& id) const
{
    const OctNode& node=_octNodes[_dProp->_ownerNodes[off]];
    ASSERT(id == node._id+_prop->_childOffset[node._level*8+_dProp->_ownerVertIds[off]]);

    sizeType voff;
    ASSERT(getVertOffset(id,voff))
    ASSERT(voff == off);
}
void OctreeDistance::checkNodeNeigh(const sizeType& i) const
{
    const OctNode& node=_octNodes[i];
    if(node._childs[0] >= 0)
        return;

    const Vec3SLC dx=_prop->_cellSz*(1<<node._level);

#define PARITY_CHECK_1(from,to,n,p,dxx)												\
{																					\
	ASSERT(_dProp->_neighDxs[-node._childs[from]*6+p] <= dxx)						\
	if(_dProp->_neighDxs[-node._childs[from]*6+p] == dxx)							\
		ASSERT(_dProp->_neighOffs[-node._childs[from]*6+p] == -node._childs[to])	\
	else																			\
		ASSERT(_dProp->_neighOffs[-node._childs[from]*6+p] != -node._childs[to])	\
																					\
	ASSERT(_dProp->_neighDxs[-node._childs[to]*6+n] <= dxx)							\
	if(_dProp->_neighDxs[-node._childs[to]*6+n] == dxx)								\
		ASSERT(_dProp->_neighOffs[-node._childs[to]*6+n] == -node._childs[from])	\
	else																			\
		ASSERT(_dProp->_neighOffs[-node._childs[to]*6+n] != -node._childs[from])	\
}

    PARITY_CHECK_1(0,1,NX,PX,dx.x())
    PARITY_CHECK_1(1,3,NY,PY,dx.y())
    PARITY_CHECK_1(2,3,NX,PX,dx.x())
    PARITY_CHECK_1(0,2,NY,PY,dx.y())

    PARITY_CHECK_1(4,5,NX,PX,dx.x())
    PARITY_CHECK_1(5,7,NY,PY,dx.y())
    PARITY_CHECK_1(6,7,NX,PX,dx.x())
    PARITY_CHECK_1(4,6,NY,PY,dx.y())

    PARITY_CHECK_1(0,4,NZ,PZ,dx.z())
    PARITY_CHECK_1(1,5,NZ,PZ,dx.z())
    PARITY_CHECK_1(3,7,NZ,PZ,dx.z())
    PARITY_CHECK_1(2,6,NZ,PZ,dx.z())

    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+0],-node._childs[0]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+1],-node._childs[1]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+2],-node._childs[2]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+3],-node._childs[3]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+4],-node._childs[4]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+5],-node._childs[5]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+6],-node._childs[6]);
    checkVertNeigh(node._id+_prop->_childOffset[node._level*8+7],-node._childs[7]);

#undef PARITY_CHECK_1
}
void OctreeDistance::checkVertNeigh(const Vec3i& id,const sizeType& valueOff) const
{
    sizeType neighOff;
    sizeType dist;

#define PARITY_CHECK_2(dir,d,dxx)							\
checkClosestVert(id,dir,neighOff,dist);						\
ASSERT(_dProp->_neighOffs[valueOff*6+d] == neighOff)		\
if(neighOff >= 0)											\
	ASSERT(_dProp->_neighDxs [valueOff*6+d] == dist*dxx)	\
else														\
	ASSERT(_dProp->_neighDxs [valueOff*6+d] == 0.0f)

    PARITY_CHECK_2(Vec3i(-1,0,0),NX,_prop->_cellSz.x())
    PARITY_CHECK_2(Vec3i( 1,0,0),PX,_prop->_cellSz.x())
    PARITY_CHECK_2(Vec3i(0,-1,0),NY,_prop->_cellSz.y())
    PARITY_CHECK_2(Vec3i(0, 1,0),PY,_prop->_cellSz.y())
    PARITY_CHECK_2(Vec3i(0,0,-1),NZ,_prop->_cellSz.z())
    PARITY_CHECK_2(Vec3i(0,0, 1),PZ,_prop->_cellSz.z())

#undef PARITY_CHECK_2
}
void OctreeDistance::checkClosestVert(const Vec3i& id,const Vec3i& dir,sizeType& neighOff,sizeType& dist) const
{
    neighOff=-1;
    for(dist=1; dist<=(sizeType(1)<<_prop->_level); dist++) {
        Vec3i neighId=id+dir*dist;

        //check for a corner
        bool found=false;
        for(sizeType nc=0; nc<8; nc++) {
            sizeType neighNodeId=findLeave((Vec3i)(neighId+_dProp->_neighOffBeforeOwner[nc]));
            if(neighNodeId >= 0) {
                const OctNode& node=_octNodes[neighNodeId];
                const sizeType exp=1<<node._level;
                const Vec3i diff=neighId-node._id;

                if((diff.x() == 0 || diff.x() == exp) &&
                        (diff.y() == 0 || diff.y() == exp) &&
                        (diff.z() == 0 || diff.z() == exp)) {
                    found=true;
                    neighOff=-node._childs[(diff>>node._level).dot(Vec3i(1,2,4))];
                    break;
                }
            }
        }
        if(found)
            return;

        //then at least we should be on an edge
        bool onEdge=false;
        for(sizeType nc=0; nc<8; nc++) {
            sizeType neighNodeId=findLeave((Vec3i)(neighId+_dProp->_neighOffBeforeOwner[nc]));
            if(neighNodeId >= 0) {
                const OctNode& node=_octNodes[neighNodeId];
                const sizeType exp=1<<node._level;
                const Vec3i diff=neighId-node._id;

                sizeType match=0;
                if(diff.x() == 0 || diff.x() == exp)match++;
                if(diff.y() == 0 || diff.y() == exp)match++;
                if(diff.z() == 0 || diff.z() == exp)match++;
                if(match == 2)
                    onEdge=true;
            }
        }
        if(!onEdge)	//then we are interior to face or node, exit
            return;
    }
}
bool OctreeDistance::vertexFastMarching(VertIter& vi,bool pCheckNeigh,bool pCheckFM)
{
    buildNeighInfo(vi,pCheckNeigh);
    initKnown();

    bool ret=Redistancing(_values,
                          _dProp->_ownerVertIds,
                          _dProp->_neighOffs,
                          _dProp->_neighDxs).fastMarch();
    if(pCheckFM) {
        const sizeType nrValue=(sizeType)_values.size();
        for(sizeType i=1; i<nrValue; i++) {
            ASSERT(_dProp->_ownerVertIds[i] == Redistancing::KNOWN)
            ASSERT(_values[i] < ScalarUtil<scalarSLC>::scalar_max)
        }

        //check maximum error
        scalarSLC error=0.0f;
        scalarSLC avgError=0.0f;
        sizeType nrPoints=0;

        const sizeType nrNode=(sizeType)_octNodes.size();
        for(sizeType i=0; i<nrNode; i++) {
            OctNode& node=_octNodes[i];
            if(node._childs[0] < 0) {
                Vec3SLC ext=node._bb.getExtent();
                Vec3SLC grad(0.0f,0.0f,0.0f);
                grad.x()=((_values[-node._childs[1]]-_values[-node._childs[0]])+
                          (_values[-node._childs[3]]-_values[-node._childs[2]])+
                          (_values[-node._childs[5]]-_values[-node._childs[4]])+
                          (_values[-node._childs[7]]-_values[-node._childs[6]]))*0.25f/ext.x();
                grad.y()=((_values[-node._childs[2]]-_values[-node._childs[0]])+
                          (_values[-node._childs[3]]-_values[-node._childs[1]])+
                          (_values[-node._childs[6]]-_values[-node._childs[4]])+
                          (_values[-node._childs[7]]-_values[-node._childs[5]]))*0.25f/ext.y();
                grad.z()=((_values[-node._childs[4]]-_values[-node._childs[0]])+
                          (_values[-node._childs[5]]-_values[-node._childs[1]])+
                          (_values[-node._childs[6]]-_values[-node._childs[2]])+
                          (_values[-node._childs[7]]-_values[-node._childs[3]]))*0.25f/ext.z();
                error=max<scalarSLC>(error,std::abs(grad.norm()-1.0f));
                avgError+=std::abs(grad.norm()-1.0f);
                nrPoints++;
            }
        }
        INFOV("Max Error: %f, Avg Error: %f",error,avgError/nrPoints)
    }
    return ret;
}
void OctreeDistance::initKnown()
{
    const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts=_mesh->_verts;
    const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors=_mesh->_nors;
    const vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors=_mesh->_tnors;
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds=_mesh->_inds;
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& neighs=_mesh->_neighs;
    const sizeType nrT=(sizeType)inds.size();
    const sizeType nrVert=(sizeType)_values.size();

    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrVert; i++) {
        _dProp->_ownerVertIds[i]&=~15;
        _dProp->_ownerVertIds[i]|=Redistancing::UNKNOWN;
    }

    //initialize by triangle intersection test
    bool hasKnown=false;
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrT; i++) {
        hasKnown|=tagIntersectTri(verts[inds[i].x()],verts[inds[i].y()],verts[inds[i].z()]);
    }

    //extreme case: there're no known points
    //we set all points to known
    if(!hasKnown) {
        const sizeType nrValue=_values.size();
        const scalarSLC maxLen=_prop->_bb.getExtent().norm();

        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrValue; i++) {
            _values[i]=maxLen;
            _dProp->_ownerVertIds[i]|=Redistancing::KNOWN;
        }
    }

    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrT; i++)
        buildInitialDist(verts[inds[i].x()],//vertex angled weighted normals
                         verts[inds[i].y()],
                         verts[inds[i].z()],
                         tnors[i],//triangle normals
                         nors[inds[i].x()],//vertex normals
                         nors[inds[i].y()],
                         nors[inds[i].z()],
                         neighs[i].x() >= 0 ? (tnors[i]+tnors[neighs[i].x()]) : tnors[i],//edge normal direction not normalized (0,1),(1,2),(2,0)
                         neighs[i].y() >= 0 ? (tnors[i]+tnors[neighs[i].y()]) : tnors[i],
                         neighs[i].z() >= 0 ? (tnors[i]+tnors[neighs[i].z()]) : tnors[i]);

    //raincheck with visibility
    {
        const sizeType nrValue=_values.size();
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrValue; i++){
            if((_dProp->_ownerVertIds[i]&15) == Redistancing::KNOWN)
                checkVisibility(i,_values[i]);
        }
    }
}
bool OctreeDistance::tagIntersectTri(const Vec3SLC& a,const Vec3SLC& b,const Vec3SLC& c)
{
    //extent one cell to avoid too close artifact
    const Vec3SLC minC=compMin(compMin(a,b),c)-_prop->_bb._minC;
    const Vec3SLC maxC=compMax(compMax(a,b),c)-_prop->_bb._minC;
    const Vec3i l=compMax((Vec3i)(floor((Vec3SLC)(minC.array()/_prop->_cellSz.array()).matrix())),Vec3i::Zero());
    const Vec3i h=compMin((Vec3i)ceil((Vec3SLC)(maxC.array()/_prop->_cellSz.array()).matrix()),_prop->_nrCellDim);

    //tag known points
    bool ret=false;
    for(sizeType x=l.x(); x<=h.x(); x++)
        for(sizeType y=l.y(); y<=h.y(); y++)
            for(sizeType z=l.z(); z<=h.z(); z++) {
                sizeType nodeId=findLeave(Vec3i(x,y,z));
                if(nodeId >= 0) {
                    const OctNode& node=_octNodes[nodeId];
                    BBox<scalarSLC> bbExp=node._bb.enlarge(_dProp->_boxEps);
                    if(TriangleTpl<scalarSLC>(a,b,c).intersect(bbExp)) {
                        ret=true;
                        OMP_CRITICAL_ {
                            _dProp->_ownerVertIds[-node._childs[0]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[0]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[1]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[1]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[2]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[2]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[3]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[3]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[4]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[4]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[5]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[5]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[6]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[6]]|=Redistancing::KNOWN;
                            _dProp->_ownerVertIds[-node._childs[7]]&=~15;
                            _dProp->_ownerVertIds[-node._childs[7]]|=Redistancing::KNOWN;
                        }
                    }
                }
            }
    return ret;
}
void OctreeDistance::buildInitialDist(const Vec3SLC& a,const Vec3SLC& b,const Vec3SLC& c,const Vec3SLC& n,
                                      const Vec3SLC& na,const Vec3SLC& nb,const Vec3SLC& nc,
                                      const Vec3SLC& ena,const Vec3SLC& enb,const Vec3SLC& enc)
{
    const BBox<scalarSLC> bbTri(compMin(compMin(a,b),c),compMax(compMax(a,b),c));
    const Vec3i l=compMax((Vec3i)(floor((Vec3SLC)((bbTri._minC-_prop->_bb._minC).array()/_prop->_cellSz.array()).matrix())-Vec3i::Constant(_dProp->_expInit)),Vec3i::Zero());
    const Vec3i h=compMin((Vec3i)(ceil ((Vec3SLC)((bbTri._maxC-_prop->_bb._minC).array()/_prop->_cellSz.array()).matrix())+Vec3i::Constant(_dProp->_expInit)),_prop->_nrCellDim);

    //update known values
    for(sizeType x=l.x(); x<=h.x(); x++)
        for(sizeType y=l.y(); y<=h.y(); y++)
            for(sizeType z=l.z(); z<=h.z(); z++) {
                sizeType valueOff;
                if(getVertOffset(Vec3i(x,y,z),valueOff)) {
                    OMP_CRITICAL_ {
                        if((_dProp->_ownerVertIds[valueOff]&15) == Redistancing::KNOWN) {
                            Vec3SLC cpt;
                            updatePointToTriDist(a,b,c,n,
                            na,nb,nc,
                            ena,enb,enc,
                            bbTri,getPt(Vec3i(x,y,z)),
                            _values[valueOff],cpt);
                        }
                    }
                }
            }
}
void OctreeDistance::polygonize(MeshSLC& mesh,bool pCheck)
{
    //generate verts
    VertDistanceSet setDist(*this,static_cast<const OctreeDistance&>(*_another),*_velCalc);
    buildNeighInfo(setDist,false);

    Polygonizer(*_velCalc,*this,_octNodes,_dProp->_neighOffs,
                static_cast<const OctreeDistance&>(*_another),mesh).polygonize(pCheck);
}
void OctreeDistance::coarsen()
{
    const sizeType nrNode=(sizeType)_octNodes.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrNode; i++)
        _octNodes[i]._valueOff=0;

    //aggregate triangle ids
    for(sizeType l=1; l<=_prop->_level; l++) {
        OMP_PARALLEL_FOR_
        for(sizeType i=_nodePerLevelBeg[l]; i<_nodePerLevelEnd[l]; i++)
            aggregateTriIds(_octNodes[i]);

        OMP_PARALLEL_FOR_
        for(sizeType i=_nodePerLevelBeg[l]; i<_nodePerLevelEnd[l]; i++) {
            const OctNode& node=_octNodes[i];
            if(node._childs[0] >= 0 && decideDelete(node)) {
                _octNodes[node._childs[0]]._valueOff=1;
                _octNodes[node._childs[1]]._valueOff=1;
                _octNodes[node._childs[2]]._valueOff=1;
                _octNodes[node._childs[3]]._valueOff=1;
                _octNodes[node._childs[4]]._valueOff=1;
                _octNodes[node._childs[5]]._valueOff=1;
                _octNodes[node._childs[6]]._valueOff=1;
                _octNodes[node._childs[7]]._valueOff=1;
            }
        }
    }

    //compact
    _dProp->_neighOffs.resize(max<sizeType>(_dProp->_neighOffs.size(),nrNode));
    fill(_dProp->_neighOffs.begin(),_dProp->_neighOffs.end(),-1);
    fill(_roots.begin(),_roots.end(),-1);
    fill(_nodePerLevelBeg.begin(),_nodePerLevelBeg.end(),-1);
    fill(_nodePerLevelEnd.begin(),_nodePerLevelEnd.end(),-1);
    sizeType j=0;
    sizeType currLevel=-1;
    sizeType rootId=0;
    for(sizeType i=0; i<nrNode; i++) {
        if(_octNodes[i]._valueOff == 0) {
            if(_octNodes[i]._level != currLevel) {
                if(currLevel >= 0)
                    _nodePerLevelEnd[currLevel]=j;
                currLevel=_octNodes[i]._level;
                _nodePerLevelBeg[currLevel]=j;
            }

            if(_octNodes[i]._level == _prop->_level)
                _roots[rootId++]=j;

            _dProp->_neighOffs[i]=j;
            _octNodes[j]=_octNodes[i];
            j++;
        }
    }
    _nodePerLevelEnd[currLevel]=j;
    _octNodes.resize(j);
    ASSERT(rootId == _prop->_nrBlock)

    INFOV("%ludu Nodes Before Coarsen, %ludu After!",nrNode,j)

    //remap
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<j; i++) {
        OctNode& node=_octNodes[i];
        if(node._childs[0] >= 0) {
            if(_dProp->_neighOffs[node._childs[0]] < 0) {
                node._childs[0]=-1;
                continue;
            }

            node._childs[0]=_dProp->_neighOffs[node._childs[0]];
            node._childs[1]=_dProp->_neighOffs[node._childs[1]];
            node._childs[2]=_dProp->_neighOffs[node._childs[2]];
            node._childs[3]=_dProp->_neighOffs[node._childs[3]];
            node._childs[4]=_dProp->_neighOffs[node._childs[4]];
            node._childs[5]=_dProp->_neighOffs[node._childs[5]];
            node._childs[6]=_dProp->_neighOffs[node._childs[6]];
            node._childs[7]=_dProp->_neighOffs[node._childs[7]];
        }
    }
}
void OctreeDistance::aggregateTriIds(OctNode& node) const
{
    if(node._childs[0] < 0)
        node._flag=false;
    else
        node._flag=_octNodes[node._childs[0]]._flag | _octNodes[node._childs[1]]._flag |
                   _octNodes[node._childs[2]]._flag | _octNodes[node._childs[3]]._flag |
                   _octNodes[node._childs[4]]._flag | _octNodes[node._childs[5]]._flag |
                   _octNodes[node._childs[6]]._flag | _octNodes[node._childs[7]]._flag;
}
bool OctreeDistance::decideDelete(const OctNode& node) const
{
    const sizeType noff=1<<node._level;
    for(sizeType x=-1; x<=1; x++)
        for(sizeType y=-1; y<=1; y++)
            for(sizeType z=-1; z<=1; z++) {
                const Vec3i nid=node._id+Vec3i(x,y,z)*noff;
                const sizeType off=findLeave(nid,node._level);
                if(off >= 0 && _octNodes[off]._flag)
                    return false;
            }
    return true;
}

//--------------------------------------------------------------------------------------------------------duplication free octree polygonizer (marching tetrahedron) algorithm
Polygonizer::Polygonizer(const VelCalc& vc,
                         OctreeDistance& currentTree,
                         vector<Octree::OctNode,Eigen::aligned_allocator<Octree::OctNode> >& nodes,
                         vector<sizeType,Eigen::aligned_allocator<sizeType> >& infoCache,
                         const OctreeDistance& referenceTree,
                         MeshSLC& mesh)
    :_vc(vc),
     _currentTree(currentTree),
     _nodes(nodes),
     _infoCache(infoCache),
     _referenceTree(referenceTree),
     _mesh(mesh),
     _beg(_currentTree.nodePerLevelBeg()[0]),
     _end(_currentTree.nodePerLevelEnd()[0])
{
    _secantTol=1E-2f*_currentTree.prop()._cellSz.minCoeff();
    _maxIter=100;

    //neigh off
    _neighOff[0]=Vec3i(-1,0,0);
    _neighOff[1]=Vec3i( 1,0,0);
    _neighOff[2]=Vec3i(0,-1,0);
    _neighOff[3]=Vec3i(0, 1,0);
    _neighOff[4]=Vec3i(0,0,-1);
    _neighOff[5]=Vec3i(0,0, 1);

    //face neighbors
    _tetIndices[0]=Vec4i(0,1,2,4);
    _tetIndices[1]=Vec4i(1,2,4,5);
    _tetIndices[2]=Vec4i(4,6,5,2);
    _tetIndices[3]=Vec4i(2,1,3,5);
    _tetIndices[4]=Vec4i(2,5,3,6);
    _tetIndices[5]=Vec4i(6,7,5,3);

    _neighIndices[0][0]=Vec3i(SELF,1,3);
    _neighIndices[0][1]=Vec3i(NX  ,3,0);
    _neighIndices[0][2]=Vec3i(NY  ,4,1);
    _neighIndices[0][3]=Vec3i(NZ  ,2,3);

    _neighIndices[1][0]=Vec3i(SELF,2,1);
    _neighIndices[1][1]=Vec3i(NY  ,5,2);
    _neighIndices[1][2]=Vec3i(SELF,3,2);
    _neighIndices[1][3]=Vec3i(SELF,0,0);

    _neighIndices[2][0]=Vec3i(SELF,4,2);
    _neighIndices[2][1]=Vec3i(SELF,1,0);
    _neighIndices[2][2]=Vec3i(NX  ,5,0);
    _neighIndices[2][3]=Vec3i(PZ  ,0,3);

    _neighIndices[3][0]=Vec3i(PX  ,0,1);
    _neighIndices[3][1]=Vec3i(SELF,4,3);
    _neighIndices[3][2]=Vec3i(SELF,1,2);
    _neighIndices[3][3]=Vec3i(NZ  ,5,3);

    _neighIndices[4][0]=Vec3i(SELF,5,1);
    _neighIndices[4][1]=Vec3i(PY  ,0,2);
    _neighIndices[4][2]=Vec3i(SELF,2,0);
    _neighIndices[4][3]=Vec3i(SELF,3,1);

    _neighIndices[5][0]=Vec3i(PX  ,2,2);
    _neighIndices[5][1]=Vec3i(SELF,4,0);
    _neighIndices[5][2]=Vec3i(PY  ,1,1);
    _neighIndices[5][3]=Vec3i(PZ  ,3,3);

    //six edge:
    //01,02,03,12,13,23
    //ab,ac,ad,bc,bd,cd
    // 0, 1, 2, 3, 4, 5

    //(0,1,2,4)
    _vertIndices[0][0]=Vec2i(SELF,0);//01
    _vertIndices[0][1]=Vec2i(SELF,1);//02
    _vertIndices[0][2]=Vec2i(SELF,2);//04
    _vertIndices[0][3]=Vec2i(SELF,4);//12
    _vertIndices[0][4]=Vec2i(SELF,3);//14
    _vertIndices[0][5]=Vec2i(SELF,5);//24

    //(1,2,4,5)
    _vertIndices[1][0]=Vec2i(SELF,4);//12
    _vertIndices[1][1]=Vec2i(SELF,3);//14
    _vertIndices[1][2]=Vec2i(1   ,2);//15
    _vertIndices[1][3]=Vec2i(SELF,5);//24
    _vertIndices[1][4]=Vec2i(SELF,6);//25
    _vertIndices[1][5]=Vec2i(4   ,0);//45

    //(4,6,5,2)
    _vertIndices[2][0]=Vec2i(4   ,1);//46
    _vertIndices[2][1]=Vec2i(4   ,0);//45
    _vertIndices[2][2]=Vec2i(SELF,5);//42
    _vertIndices[2][3]=Vec2i(4   ,4);//65
    _vertIndices[2][4]=Vec2i(2   ,2);//62
    _vertIndices[2][5]=Vec2i(SELF,6);//52

    //(2,1,3,5)
    _vertIndices[3][0]=Vec2i(SELF,4);//21
    _vertIndices[3][1]=Vec2i(2   ,0);//23
    _vertIndices[3][2]=Vec2i(SELF,6);//25
    _vertIndices[3][3]=Vec2i(1   ,1);//13
    _vertIndices[3][4]=Vec2i(1   ,2);//15
    _vertIndices[3][5]=Vec2i(1   ,5);//35

    //(2,5,3,6)
    _vertIndices[4][0]=Vec2i(SELF,6);//25
    _vertIndices[4][1]=Vec2i(2   ,0);//23
    _vertIndices[4][2]=Vec2i(2   ,2);//26
    _vertIndices[4][3]=Vec2i(1   ,5);//53
    _vertIndices[4][4]=Vec2i(4   ,4);//56
    _vertIndices[4][5]=Vec2i(2   ,3);//36

    //(6,7,5,3)
    _vertIndices[5][0]=Vec2i(6   ,0);//67
    _vertIndices[5][1]=Vec2i(4   ,4);//65
    _vertIndices[5][2]=Vec2i(2   ,3);//63
    _vertIndices[5][3]=Vec2i(5   ,1);//75
    _vertIndices[5][4]=Vec2i(3   ,2);//73
    _vertIndices[5][5]=Vec2i(1   ,5);//53

    //ab,ac,ad,bc,bd,cd
    // 0, 1, 2, 3, 4, 5
    //mc index
    //0000
    _mcNr[0]=0;
    //0001
    _mcNr[1]=1;
    _mcIndexEdge[1][0]=Vec3i(4,2,5);
    //0010
    _mcNr[2]=1;
    _mcIndexEdge[2][0]=Vec3i(1,3,5);
    //0011
    _mcNr[3]=2;
    _mcIndexEdge[3][0]=Vec3i(2,1,3);
    _mcIndexEdge[3][1]=Vec3i(2,3,4);
    //0100
    _mcNr[4]=1;
    _mcIndexEdge[4][0]=Vec3i(0,4,3);
    //0101
    _mcNr[5]=2;
    _mcIndexEdge[5][0]=Vec3i(2,5,3);
    _mcIndexEdge[5][1]=Vec3i(2,3,0);
    //0110
    _mcNr[6]=2;
    _mcIndexEdge[6][0]=Vec3i(0,4,5);
    _mcIndexEdge[6][1]=Vec3i(0,5,1);
    //0111
    _mcNr[7]=1;
    _mcIndexEdge[7][0]=Vec3i(0,2,1);
    //1000
    _mcNr[8]=1;
    _mcIndexEdge[8][0]=Vec3i(0,1,2);
    //1001
    _mcNr[9]=2;
    _mcIndexEdge[9][0]=Vec3i(0,1,5);
    _mcIndexEdge[9][1]=Vec3i(0,5,4);
    //1010
    _mcNr[10]=2;
    _mcIndexEdge[10][0]=Vec3i(0,3,5);
    _mcIndexEdge[10][1]=Vec3i(0,5,2);
    //1011
    _mcNr[11]=1;
    _mcIndexEdge[11][0]=Vec3i(0,3,4);
    //1100
    _mcNr[12]=2;
    _mcIndexEdge[12][0]=Vec3i(2,4,3);
    _mcIndexEdge[12][1]=Vec3i(2,3,1);
    //1101
    _mcNr[13]=1;
    _mcIndexEdge[13][0]=Vec3i(5,3,1);
    //1110
    _mcNr[14]=1;
    _mcIndexEdge[14][0]=Vec3i(4,5,2);
    //1111
    _mcNr[15]=0;

    //index edge
    _idEdge[0]=Vec2i(0,1);
    _idEdge[1]=Vec2i(0,2);
    _idEdge[2]=Vec2i(0,3);
    _idEdge[3]=Vec2i(1,2);
    _idEdge[4]=Vec2i(1,3);
    _idEdge[5]=Vec2i(2,3);

    //insert face index
    for(sizeType i=0; i<16; i++) {
        _mcIndexFaceInv[i]=Vec4i::Constant(15);

        const sizeType nr=_mcNr[i];
        for(sizeType n=0; n<nr; n++) {
            _mcIndexFace[i][n](0)=findFace(_mcIndexEdge[i][n](0),_mcIndexEdge[i][n](1));
            _mcIndexFace[i][n](1)=findFace(_mcIndexEdge[i][n](1),_mcIndexEdge[i][n](2));
            _mcIndexFace[i][n](2)=findFace(_mcIndexEdge[i][n](2),_mcIndexEdge[i][n](0));

            if(_mcIndexFace[i][n](0) != -1)
                _mcIndexFaceInv[i](_mcIndexFace[i][n](0))=n;
            if(_mcIndexFace[i][n](1) != -1)
                _mcIndexFaceInv[i](_mcIndexFace[i][n](1))=n;
            if(_mcIndexFace[i][n](2) != -1)
                _mcIndexFaceInv[i](_mcIndexFace[i][n](2))=n;
        }
    }
}
void Polygonizer::polygonize(bool pCheck)
{
    //allocate tet key and face indices
    //for every leaf cell, we store 4 indices
    //1 for tet key and other 3 for face indices
    const sizeType nrLeaves=_end-_beg;
    _infoCache.resize(max<sizeType>((sizeType)_infoCache.size(),nrLeaves*4));
    fill(_infoCache.begin(),_infoCache.end(),-1);

    //count how many vertices and polygons
    //that will be generated in each leave
    OMP_PARALLEL_FOR_
    for(sizeType i=_beg; i<_end; i++)
        countVertPoly(i,i-_beg);

    //count and allocate stored in _triIds
    sizeType nrVert=0;
    sizeType nrTri=0;
    for(sizeType i=_beg; i<_end; i++) {
        Octree::OctNode& node=_nodes[i];
        ASSERT(node._level == 0)

        node._triIds.x()=nrVert;
        node._triIds.y()=nrTri;

        nrVert+=countBits<sizeType>(node._flag);
        nrTri+=_infoCache[(i-_beg)*4]&15;
    }

    //allocate space
    _mesh._verts.resize(nrVert);
    _mesh._nors.resize(nrVert);
    _mesh._tnors.resize(nrTri);
    _mesh._inds.resize(nrTri);
    _mesh._neighs.resize(nrTri);

    //generate vertices
    OMP_PARALLEL_FOR_
    for(sizeType i=_beg; i<_end; i++)
        genVertPoly(i,i-_beg);

    //generate neigh info
    OMP_PARALLEL_FOR_
    for(sizeType i=_beg; i<_end; i++)
        genNeighInfo(i,i-_beg);

    //generate triIds
    OMP_PARALLEL_FOR_
    for(sizeType i=_beg; i<_end; i++)
        genTriIds(i,i-_beg);

    //build vertex and triangle normal
    fill(_mesh._nors.begin(),_mesh._nors.end(),Vec3SLC::Constant(0.0f));
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrTri; i++) {
        const Vec3SLC& a=_mesh._verts[_mesh._inds[i].x()];
        const Vec3SLC& b=_mesh._verts[_mesh._inds[i].y()];
        const Vec3SLC& c=_mesh._verts[_mesh._inds[i].z()];
        //triangle normal
        _mesh._tnors[i]=(b-a).normalized().cross((c-a).normalized()).normalized();
        //vertex normal
        _mesh._nors[_mesh._inds[i].x()]+=_mesh._tnors[i]*getAngle3D<scalarSLC>(b-a,c-a);
        _mesh._nors[_mesh._inds[i].y()]+=_mesh._tnors[i]*getAngle3D<scalarSLC>(c-b,a-b);
        _mesh._nors[_mesh._inds[i].z()]+=_mesh._tnors[i]*getAngle3D<scalarSLC>(a-c,b-c);
    }
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrVert; i++)
        _mesh._nors[i].normalize();

    if(pCheck) {
        for(sizeType i=0; i<nrTri; i++) {
            ASSERT(_mesh._neighs[i].x() != i)
            ASSERT(_mesh._neighs[i].y() != i)
            ASSERT(_mesh._neighs[i].z() != i)
            ASSERT(checkTriNeigh(_mesh._inds[i].x(),_mesh._inds[i].y(),_mesh._inds[_mesh._neighs[i].x()]))
            ASSERT(checkTriNeigh(_mesh._inds[i].y(),_mesh._inds[i].z(),_mesh._inds[_mesh._neighs[i].y()]))
            ASSERT(checkTriNeigh(_mesh._inds[i].z(),_mesh._inds[i].x(),_mesh._inds[_mesh._neighs[i].z()]))
        }
        for(sizeType i=0; i<nrVert; i++)
            ASSERT(std::abs(_mesh._nors[i].norm()-1.0f) < EPS)
        }
}
void Polygonizer::countVertPoly(const sizeType& i,const sizeType& idLeave)
{
#define DECIDE_VERT(id,mask,a,b)	\
if(neg[a]^neg[b])						\
{										\
	node._flag|=mask;					\
	compressedVID|=(index<<(id*4));		\
	index++;							\
}										\
else									\
{										\
	compressedVID|=(15<<(id*4));		\
}

#define TEST_VID(mask,bits)							\
if(node._flag & mask)								\
{													\
	ASSERT(((compressedVID>>bits)&15) == index);	\
	index++;										\
}													\
else												\
{													\
	ASSERT(((compressedVID>>bits)&15) == 15);		\
}

    Octree::OctNode& node=_nodes[i];
    ASSERT(node._level == 0)

    bool neg[8];
    neg[0]=_currentTree.values()[-node._childs[0]] < 0.0f;
    neg[1]=_currentTree.values()[-node._childs[1]] < 0.0f;
    neg[2]=_currentTree.values()[-node._childs[2]] < 0.0f;
    neg[3]=_currentTree.values()[-node._childs[3]] < 0.0f;
    neg[4]=_currentTree.values()[-node._childs[4]] < 0.0f;
    neg[5]=_currentTree.values()[-node._childs[5]] < 0.0f;
    neg[6]=_currentTree.values()[-node._childs[6]] < 0.0f;
    neg[7]=_currentTree.values()[-node._childs[7]] < 0.0f;

    node._flag=0;
    node._valueOff=0;

    //count how many triangles to generate
    //valueOff is used as following:
    //0 ~3 bit how many tri to generate
    //4 ~7 bit case for first tet
    //8 ~11bit case for second tet
    //12~15bit case for third tet
    //16~19bit case for forth tet
    //20~23bit case for fifth tet
    //24~27bit case for sixth tet
    sizeType& compressedTK=_infoCache[idLeave*4];
    compressedTK=0;
    sizeType nrTri=0;
    for(sizeType tet=0; tet<6; tet++) {
        unsigned char tetCase=0;
        if(neg[_tetIndices[tet](0)])tetCase|=8;
        if(neg[_tetIndices[tet](1)])tetCase|=4;
        if(neg[_tetIndices[tet](2)])tetCase|=2;
        if(neg[_tetIndices[tet](3)])tetCase|=1;
        compressedTK|=tetCase<<((tet+1)*4);
        nrTri+=_mcNr[tetCase];
    }
    compressedTK|=nrTri;

    //count how many verts to generate
    sizeType& compressedVID=node._valueOff;
    compressedVID=0;
    sizeType index=0;
    DECIDE_VERT(0,1 ,0,1)
    DECIDE_VERT(1,2 ,0,2)
    DECIDE_VERT(2,4 ,0,4)
    DECIDE_VERT(3,8 ,4,1)
    DECIDE_VERT(4,16,1,2)
    DECIDE_VERT(5,32,2,4)
    DECIDE_VERT(6,64,2,5)

    index=0;
    TEST_VID(1 ,0)
    TEST_VID(2 ,4)
    TEST_VID(4 ,8)
    TEST_VID(8 ,12)
    TEST_VID(16,16)
    TEST_VID(32,20)
    TEST_VID(64,24)

#undef DECIDE_VERT
}
void Polygonizer::genVertPoly(const sizeType& i,const sizeType& idLeave)
{
#define GEN_VERT(mask,a,b)																					 \
if(node._flag & mask)																						 \
{																											 \
    const Vec3SLC p0=_currentTree.getPt(node._id+_currentTree.prop()._childOffset[node._level*8+a]);		 \
    const Vec3SLC p1=_currentTree.getPt(node._id+_currentTree.prop()._childOffset[node._level*8+b]);		 \
	verts[index]=secantSearch(0,1,																			 \
							  _currentTree.values()[-node._childs[a]],										 \
							  _currentTree.values()[-node._childs[b]],p0,p1);								 \
	index++;																								 \
}

#define GEN_NEIGH_FACE(id,bits)				\
if(faceInv(id) != 15)						\
	neighFace|=(index+faceInv(id))<<bits;	\
else										\
	neighFace|=15<<bits;

    Octree::OctNode& node=_nodes[i];
    ASSERT(node._level == 0)

    //fill verts
    if(node._flag) {
        Vec3SLC* verts=&(_mesh._verts[node._triIds.x()]);
        sizeType index=0;
        GEN_VERT(1 ,0,1)
        GEN_VERT(2 ,0,2)
        GEN_VERT(4 ,0,4)
        GEN_VERT(8 ,4,1)
        GEN_VERT(16,1,2)
        GEN_VERT(32,2,4)
        GEN_VERT(64,2,5)
    }

    //fill inds
    const sizeType tetInfo=_infoCache[idLeave*4];
    if(tetInfo&15) {
        Vec3i* inds=&(_mesh._inds[node._triIds.y()]);
        unsigned short* neighFaces=(unsigned short*)&(_infoCache[idLeave*4+1]);

        sizeType index=0;
        sizeType tetKey=tetInfo>>4;
        for(sizeType tet=0; tet<6; tet++) {
            const sizeType key=tetKey&15;
            const sizeType nrTCurrTet=_mcNr[key];
            const Vec4i& faceInv=_mcIndexFaceInv[key];

            //generate indices
            unsigned short& neighFace=neighFaces[tet];
            neighFace=0;
            GEN_NEIGH_FACE(0,0)
            GEN_NEIGH_FACE(1,4)
            GEN_NEIGH_FACE(2,8)
            GEN_NEIGH_FACE(3,12)

            //generate triangles
            for(sizeType tri=0; tri<nrTCurrTet; tri++) {
                const Vec3i& edgeId=_mcIndexEdge[key][tri];
                Vec3i& ind=inds[index];
                ind(0)=findVert(node,_vertIndices[tet][edgeId(0)]);
                ASSERT(ind(0) < (sizeType)_mesh._verts.size())
                ind(1)=findVert(node,_vertIndices[tet][edgeId(1)]);
                ASSERT(ind(1) < (sizeType)_mesh._verts.size())
                ind(2)=findVert(node,_vertIndices[tet][edgeId(2)]);
                ASSERT(ind(2) < (sizeType)_mesh._verts.size())
                index++;
            }
            tetKey>>=4;
        }
        ASSERT(index == (tetInfo&15));
    }

#undef GEN_VERT
#undef GEN_NEIGH_FACE
}
void Polygonizer::genNeighInfo(const sizeType& i,const sizeType& idLeave)
{
    Octree::OctNode& node=_nodes[i];
    ASSERT(node._level == 0)

    const sizeType tetInfo=_infoCache[idLeave*4];
    if(tetInfo&15) {
        Vec3i* neighs=&(_mesh._neighs[node._triIds.y()]);
        //const Vec3i* inds=&(_mesh._inds[node._triIds.y()]);

        sizeType index=0;
        sizeType tetKey=tetInfo>>4;
        for(sizeType tet=0; tet<6; tet++) {
            const sizeType key=tetKey&15;
            const sizeType nrTCurrTet=_mcNr[key];

            //generate triangles
            for(sizeType t=0; t<nrTCurrTet; t++) {
                //const Vec3i& edgeId=_mcIndexEdge[key][t];
                const Vec3i& faceId=_mcIndexFace[key][t];

                Vec3i& neigh=neighs[index];
                //const Vec3i& ind=inds[index];
                if(faceId.x() != -1)neigh.x()=findNeighTri(node,_neighIndices[tet][faceId.x()]);
                if(faceId.y() != -1)neigh.y()=findNeighTri(node,_neighIndices[tet][faceId.y()]);
                if(faceId.z() != -1)neigh.z()=findNeighTri(node,_neighIndices[tet][faceId.z()]);
                if(nrTCurrTet == 2) {
                    if(t == 0)
                        neigh.z()=node._triIds.y()+(index+1);
                    else
                        neigh.x()=node._triIds.y()+(index-1);
                }

                index++;
            }
            tetKey>>=4;
        }
        ASSERT(index == (tetInfo&15));
    }
}
void Polygonizer::genTriIds(const sizeType& i,const sizeType& idLeave)
{
    Octree::OctNode& node=_nodes[i];
    ASSERT(node._level == 0)

    const sizeType tetInfo=_infoCache[idLeave*4];
    if(tetInfo&15) {
        node._flag=true;
        node._triIds.x()=node._triIds.y();
        node._triIds.y()+=tetInfo&15;
    } else {
        node._flag=false;
        node._triIds=Vec2i::Constant(-1);
    }
}
sizeType Polygonizer::findVert(const Octree::OctNode& node,const Vec2i& vid) const
{
    sizeType off;
    if(vid.x() == -1) {
        off=(node._valueOff>>(vid.y()*4))&15;
        ASSERT(off != 15)
        return node._triIds.x()+off;
    } else {
        const Vec3i nid=node._id+_currentTree.prop()._childOffset[node._level*8+vid.x()];
        const sizeType neighId=_currentTree.findLeave(nid);
        if(neighId == -1 || _nodes[neighId]._level != 0)
            return -1;

        off=(_nodes[neighId]._valueOff>>(vid.y()*4))&15;
        ASSERT(off != 15)
        return _nodes[neighId]._triIds.x()+off;
    }
}
sizeType Polygonizer::findNeighTri(const Octree::OctNode& node,const Vec3i& neighInfo) const
{
    //which node
    sizeType idNode;
    sizeType idLeave;	//index leave
    if(neighInfo.x() == SELF)
        idNode=(sizeType)(&node-&(_nodes[0]));
    else {
        idNode=_currentTree.findLeave((Vec3i)(node._id+_neighOff[neighInfo.x()]));
        ASSERT(idNode >= 0 && _nodes[idNode]._level == 0)
    }
    idLeave=idNode-_beg;

    //which tet
    unsigned short index=
        ((const unsigned short*)&(_infoCache[idLeave*4+1]))[neighInfo.y()];
    index=(index>>(neighInfo.z()*4))&15;
    ASSERT(index != 15)

    //return global index
    return _nodes[idNode]._triIds.y()+index;
}
Vec3SLC Polygonizer::secantSearch(scalarSLC xN1,scalarSLC xN2,
                                  scalarSLC fN1,scalarSLC fN2,
                                  const Vec3SLC& x0,const Vec3SLC& x1) const
{
    Vec3SLC sol;
    scalarSLC xN;
    scalarSLC fN=ScalarUtil<scalarSLC>::scalar_max;

    sizeType iter=0;
    while(std::abs(fN) > _secantTol && std::abs(fN1-fN2) > EPS && iter < _maxIter) {
        xN=xN1-fN1*(xN1-xN2)/(fN1-fN2);
        sol=x0*(1.0f-xN)+x1*(xN);
        fN=_referenceTree.getDist(_vc.trace(sol));

        if(fN*fN1 < 0.0f) {
            xN2=xN;
            fN2=fN;
        } else {
            xN1=xN;
            fN1=fN;
        }

        iter++;
    }

    if(iter == _maxIter) {
        WARNINGV("Exit By Max Iter, Error: %f, Solution: %f.",fN,xN)
        WARNINGV("Range: (%f %f), Solution Range: (%f %f)",fN1,fN2,xN1,xN2)
        WARNINGV("At: (%f %f %f)",sol.x(),sol.y(),sol.z())
    }

    return sol;
}
bool Polygonizer::checkTriNeigh(const sizeType& v0,const sizeType& v1,const Vec3i& b) const
{
    sizeType match=0;
    if(b.x() == v0 || b.x() == v1)match++;
    if(b.y() == v0 || b.y() == v1)match++;
    if(b.z() == v0 || b.z() == v1)match++;
    return match >= 2;
}
sizeType Polygonizer::findFace(const sizeType& a,const sizeType& b) const
{
    set<sizeType> verts;
    verts.insert(_idEdge[a].x());
    verts.insert(_idEdge[a].y());
    verts.insert(_idEdge[b].x());
    verts.insert(_idEdge[b].y());

    if(verts.find(0) == verts.end())return 0;
    if(verts.find(1) == verts.end())return 1;
    if(verts.find(2) == verts.end())return 2;
    if(verts.find(3) == verts.end())return 3;
    return -1;
}

//--------------------------------------------------------------------------------------------------------octree vertex fast marching algorithm
Redistancing::Redistancing(vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& values,
                           vector<unsigned char,Eigen::aligned_allocator<unsigned char> >& states,
                           const vector<sizeType,Eigen::aligned_allocator<sizeType> >& neighOffs,
                           const vector<scalarSLC,Eigen::aligned_allocator<scalarSLC> >& neighDxs)
    :_values(values),
     _states(states),
     _neighOffs(neighOffs),
     _neighDxs(neighDxs)
{
    ASSERT(_values.size() == _states.size())
    ASSERT(_values.size()*6 == _neighOffs.size())
    ASSERT(_values.size()*6 == _neighDxs.size())
    if((_states[1]&POSITIVE) || (_states[1]&NEGATIVE))
        INFO("Contain Sign Info!")
}
bool Redistancing::fastMarch()
{
    //initialize heap
    const sizeType nrValue=(sizeType)_values.size();
    for(sizeType i=1; i<nrValue; i++)
        tagClose(i);

    //add all close to heap
    _heap.clear();
    _heapOffsets.resize(_values.size());
    fill(_heapOffsets.begin(),_heapOffsets.end(),-1);
    for(sizeType i=1; i<nrValue; i++) {
        if((_states[i]&15) == CLOSE) {
            if(!updateClose(i))
                return false;
        }
    }

    //march until heap empty
    while(!_heap.empty()) {
        sizeType i=popHeap();
        {
            _states[i]&=~15;
            _states[i]|=KNOWN;
        }
        if(!updateKnown(i))
            return false;
    }

    return true;
}
void Redistancing::tagClose(const sizeType& i)
{
    if((_states[i]&15) == KNOWN) {
        const sizeType* neighOffs=&(_neighOffs[i*6]);
        OMP_CRITICAL_ {
#define CHECK_SET_CLOSE(id)											\
if(neighOffs[id] >= 0 && (_states[neighOffs[id]]&15) == UNKNOWN)	\
{																	\
	_states[neighOffs[id]]&=~15;									\
	_states[neighOffs[id]]|=CLOSE;									\
}

            CHECK_SET_CLOSE(NX)
            CHECK_SET_CLOSE(PX)
            CHECK_SET_CLOSE(NY)
            CHECK_SET_CLOSE(PY)
            CHECK_SET_CLOSE(NZ)
            CHECK_SET_CLOSE(PZ)

#undef CHECK_SET_CLOSE
        }
    }
}
sizeType Redistancing::popHeap()
{
    if(_heap.empty())
        return -1;

    const sizeType index=_heap[0];
    _heap[0]=_heap.back();
    _heapOffsets[_heap[0]]=0;

    const sizeType back=(sizeType)(_heap.size()-1);
    for(sizeType i=0,j; i<back-1; i=j) {
        sizeType lc=(i<<1)+1;
        sizeType rc=lc+1;
        scalarSLC current=std::abs(_values[_heap[i]]);
        scalarSLC lv,rv;
        if(lc < back) {
            lv=std::abs(_values[_heap[lc]]);
            if(rc<back) {
                rv=std::abs(_values[_heap[rc]]);
                if(lv>rv) {
                    lc=rc;
                    lv=rv;
                }
            }
            if(current>lv) {
                _heap[i]=_heap[lc];
                _heapOffsets[_heap[lc]]=i;
                _heap[lc]=_heap.back();
                _heapOffsets[_heap.back()]=lc;
                j=lc;
            } else break;
        } else break;
    }

    _heap.pop_back();
    return index;
}
void Redistancing::pushHeap(sizeType valueOff)
{
    const sizeType back=(sizeType)_heap.size();
    _heap.push_back(valueOff);
    _heapOffsets[valueOff]=back;

    for(sizeType i=back,j=(i-1)>>1; i>0; i=j,j=(i-1)>>1) {
        if(std::abs(_values[_heap[i]])<std::abs(_values[_heap[j]])) {
            _heap[i]=_heap[j];
            _heapOffsets[_heap[j]]=i;
            _heap[j]=valueOff;
            _heapOffsets[valueOff]=j;
        } else break;
    }
}
void Redistancing::updateHeap(sizeType valueOff)
{
    for(sizeType i=_heapOffsets[valueOff],j=(i-1)>>1; i>0; i=j,j=(i-1)>>1) {
        if(std::abs(_values[_heap[i]])<std::abs(_values[_heap[j]])) {
            _heap[i]=_heap[j];
            _heapOffsets[_heap[j]]=i;
            _heap[j]=valueOff;
            _heapOffsets[valueOff]=j;
        } else break;
    }
}
bool Redistancing::updateClose(const sizeType& valueOff)
{
#define CHECK_DIR(d)													\
exist=false;															\
if(neighOff[N##d] >= 0 && (_states[neighOff[N##d]]&15) == KNOWN)		\
{																		\
	A[off]=1.0f/neighDx[N##d];											\
	A[off]=A[off]*A[off];												\
	B[off]=_values[neighOff[N##d]];										\
	exist=true;															\
}																		\
if(neighOff[P##d] >= 0 && (_states[neighOff[P##d]]&15) == KNOWN)		\
{																		\
	tmp=_values[neighOff[P##d]];										\
	if(exist)															\
	{																	\
		if(neighDx[P##d] < neighDx[N##d] || std::abs(tmp) < std::abs(B[off]))		\
		{																\
			A[off]=1.0f/neighDx[P##d];									\
			A[off]=A[off]*A[off];										\
			B[off]=tmp;													\
		}																\
	}																	\
	else																\
	{																	\
		A[off]=1.0f/neighDx[P##d];										\
		A[off]=A[off]*A[off];											\
		B[off]=tmp;														\
	}																	\
	exist=true;															\
}																		\
if(exist)																\
	off++;

    const sizeType* neighOff=&(_neighOffs[valueOff*6]);
    const scalarSLC* neighDx=&(_neighDxs[valueOff*6]);

    //locals
    scalarSLC A[3];
    scalarSLC B[3];
    sizeType off=0;
    scalarSLC a,b,c,tmp,det;
    bool exist;

    //check every direction
    CHECK_DIR(X)
    CHECK_DIR(Y)
    CHECK_DIR(Z)

    //check same sign
    bool positive=B[0] > 0.0f;
    for(sizeType i=1; i<off; i++) {
        if(B[0]*B[i]<0.0f) {
            if((_states[valueOff]&POSITIVE) || (_states[valueOff]&NEGATIVE)) {
                WARNING("Resort To Original Sign!")
                positive=(_states[valueOff]&POSITIVE) != 0;
                break;
            } else {
                WARNINGV("(%lud) In-Out Test Fail For Fast Marching, Mesh Not Watertight!",valueOff)
                return false;
            }
        }
    }

    //sort the minimal direction
    if(off == 3) {
        if(std::abs(B[1]) < std::abs(B[0])) {
            swap(B[0],B[1]);
            swap(A[0],A[1]);
        }

        if(std::abs(B[2]) < std::abs(B[1])) {
            swap(B[1],B[2]);
            swap(A[1],A[2]);
        }

        if(std::abs(B[1]) < std::abs(B[0])) {
            swap(B[0],B[1]);
            swap(A[0],A[1]);
        }
    } else if(off == 2) {
        if(std::abs(B[1]) < std::abs(B[0])) {
            swap(B[0],B[1]);
            swap(A[0],A[1]);
        }
    }

    //try solving
    while(off>0) {
        //build a,b,c
        a=b=0.0f;
        c=-1.0f;
        for(sizeType i=0; i<off; i++) {
            a+=A[i];
            b-=2.0f*B[i]*A[i];
            c+=B[i]*B[i]*A[i];
        }

        det=b*b-4.0f*a*c;
        if(det < 0.0f) {
            off--;
        } else {
            //solution to the quadratic
            a=max<scalarSLC>(2.0f*a,EPS);
            _values[valueOff]=positive ? ((-b+sqrt(det))/a) : ((-b-sqrt(det))/a);
            _values[valueOff]=std::abs(_values[valueOff]);
            if(!positive)
                _values[valueOff]*=-1.0f;

            //update heap
            if(_heapOffsets[valueOff] < 0)
                pushHeap(valueOff);
            else
                updateHeap(valueOff);

            break;
        }
    }

    if(off == 0) {
        WARNING("Numerically Impossible!")
        return false;
    }
    return true;

#undef CHECK_DIR
}
bool Redistancing::updateKnown(const sizeType& i)
{
#define UPDATE_CLOSE(id)										\
if(neighOffs[id] >= 0 && (_states[neighOffs[id]]&15) != KNOWN)	\
{																\
	_states[neighOffs[id]]&=~15;								\
	_states[neighOffs[id]]|=CLOSE;								\
	if(!updateClose(neighOffs[id]))								\
		return false;											\
}

    const sizeType* neighOffs=&(_neighOffs[i*6]);
    UPDATE_CLOSE(NX)
    UPDATE_CLOSE(PX)
    UPDATE_CLOSE(NY)
    UPDATE_CLOSE(PY)
    UPDATE_CLOSE(NZ)
    UPDATE_CLOSE(PZ)
    return true;

#undef UPDATE_CLOSE
}

//--------------------------------------------------------------------------------------------------------enright test 3D
Vec3SLC VelCalcEnrightTest::vel(const Vec3SLC& pos) const
{
    const scalarSLC SPX=sin(pos.x()*M_PI),S2PX=sin(pos.x()*M_PI*2.0f);
    const scalarSLC SPY=sin(pos.y()*M_PI),S2PY=sin(pos.y()*M_PI*2.0f);
    const scalarSLC SPZ=sin(pos.z()*M_PI),S2PZ=sin(pos.z()*M_PI*2.0f);
    Vec3SLC ret=Vec3SLC(2.0f*SPX*SPX*S2PY*S2PZ,
                        -S2PX*SPY*SPY*S2PZ,
                        -S2PX*S2PY*SPZ*SPZ);
    if(_reverse)
        return ret*-1.0f;
    else
        return ret;
}
void VelCalcEnrightTest::init(MeshSLC& mesh) const
{
    static const Vec3SLC origin(0.35,0.35,0.35);
    static const scalarSLC radius=0.15f;
    static const sizeType SLICE=64;

    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& verts=mesh._verts;
    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& nors=mesh._nors;
    vector<Vec3SLC,Eigen::aligned_allocator<Vec3SLC> >& tnors=mesh._tnors;
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds=mesh._inds;
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& neighs=mesh._neighs;

    verts.clear();
    inds.clear();
    neighs.clear();

    scalarSLC angV=M_PI*0.5f;
    scalarSLC angH;
    sizeType idA,idB;

    //top
    {
        idA=verts.size();
        Vec3SLC pos=origin+Vec3SLC(0.0f,0.0f,1.0f)*radius;
        verts.push_back(pos);

        angV-=M_PI/SLICE;
        angH=0;
        idB=verts.size();
        for(sizeType i=0,j=idB; i<SLICE; i++,j++) {
            verts.push_back(origin+Vec3SLC(cos(angV)*cos(angH),cos(angV)*sin(angH),sin(angV))*radius);
            if(i == SLICE-1)
                inds.push_back(Vec3i(idA,j,idB));
            else
                inds.push_back(Vec3i(idA,j,j+1));
            angH+=M_PI*2.0f/SLICE;
        }
    }

    //vertical strip
    for(sizeType v=1; v<SLICE-1; v++) {
        angV-=M_PI/SLICE;
        angH=0;
        idA=idB;
        idB=verts.size();
        for(sizeType i=0; i<SLICE; i++) {
            verts.push_back(origin+Vec3SLC(cos(angV)*cos(angH),cos(angV)*sin(angH),sin(angV))*radius);
            angH+=M_PI*2.0f/SLICE;
        }

        for(sizeType i=0,j=idA,k=idB; i<SLICE; i++,j++,k++) {
            if(i == SLICE-1) {
                inds.push_back(Vec3i(k,idB,idA));
                inds.push_back(Vec3i(k,idA,j));
            } else {
                inds.push_back(Vec3i(k,k+1,j+1));
                inds.push_back(Vec3i(k,j+1,j));
            }
        }
    }

    //bottom
    {
        angV-=M_PI/SLICE;
        angH=0;
        idA=idB;
        idB=verts.size();
        verts.push_back(origin-Vec3SLC(0.0f,0.0f,1.0f)*radius);

        for(sizeType i=0,j=idA; i<SLICE; i++,j++) {
            if(i == SLICE-1)
                inds.push_back(Vec3i(j,idB,idA));
            else
                inds.push_back(Vec3i(j,idB,j+1));
        }
    }

    nors.resize(verts.size());
    for(sizeType i=0; i<(sizeType)verts.size(); i++)
        nors[i]=(verts[i]-origin).normalized();

    tnors.resize(inds.size());
    for(sizeType i=0; i<(sizeType)tnors.size(); i++)
        tnors[i]=(verts[inds[i].y()]-verts[inds[i].x()]).cross(verts[inds[i].z()]-verts[inds[i].x()]).normalized();

    mesh.searchNeigh();
    mesh.parityCheck();
}

//--------------------------------------------------------------------------------------------------------surface tracking
SurfaceTracker::SurfaceTracker(const Vec3SLC& cellSz,const sizeType& level)
{
    _cellSz=cellSz;
    _nrCellDim=Vec3i::Constant(1<<level);
    _bb._minC=Vec3SLC::Zero();
    _bb._maxC=Vec3SLC(cellSz.x()*_nrCellDim.x(),
                      cellSz.y()*_nrCellDim.y(),
                      cellSz.z()*_nrCellDim.z());
    _jr._minC=_cellSz*-1.0f;
    _jr._maxC=Vec3SLC::Zero();
    _valid=false;

    _helper.reset(new OctreeDistance::TreeBuildHelper);
    _treeA.reset(new OctreeDistance(_bb,_cellSz,level,*_helper));
    _treeB.reset(new OctreeDistance(_bb,_cellSz,level,*_helper));

    _meshA.reset(new MeshSLC);
    _meshB.reset(new MeshSLC);
}
void SurfaceTracker::initFromMesh(const MeshSLC& mesh)
{
    *_meshA=mesh;
    _bvh.reset(new AABBvh(_meshA.get()));
    _valid=_treeA->buildFromAABBvhTopDown(_bvh.get(),false);
}
void SurfaceTracker::advance(const VelCalc& vel)
{
    if(!_valid) {
        WARNING("Invalid Tree State!")
        return;
    }

    //perform jitter
    BBox<scalarSLC> bb=_bb;
    Vec3SLC rv(rand()/(scalarSLC)RAND_MAX,
               rand()/(scalarSLC)RAND_MAX,
               rand()/(scalarSLC)RAND_MAX);
    Vec3SLC randOff(rv.x()*_jr._minC.x()+(1.0f-rv.x())*_jr._maxC.x(),
                    rv.y()*_jr._minC.y()+(1.0f-rv.y())*_jr._maxC.y(),
                    rv.z()*_jr._minC.z()+(1.0f-rv.z())*_jr._maxC.z());
    bb._minC+=randOff;
    bb._maxC+=randOff;
    _treeB->resetBB(bb);

    _valid=_treeB->buildFromDistanceTreeTopDown(_treeA.get(),vel,_meshB.get(),false);
    _bvh.reset((AABBvh*)0);
    if(!_valid)
        WARNING("Entering Invalid State!")

        _treeA.swap(_treeB);
    _meshA.swap(_meshB);
}
void SurfaceTracker::read(istream& is)
{
    _meshA->read(is);
    _treeA->read(is);
}
void SurfaceTracker::write(ostream& os) const
{
    _meshA->write(os);
    _treeA->write(os);
}
void SurfaceTracker::writeVTK(const string& path) const
{
    _meshA->writeVTK(path);
}
