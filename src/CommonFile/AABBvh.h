#ifndef AABBVH_H
#define AABBVH_H

#include "MathBasic.h"
#include "CollisionDetection.h"
#include "IO.h"
#include <boost/filesystem.hpp>
PRJ_BEGIN

template <typename T>
class AABBvhTpl
{
public:
    typedef typename ScalarUtil<T>::ScalarVec3 PT3;
    typedef typename ScalarUtil<T>::ScalarVec2 PT2;
    struct BVHNode {
        BBox<T> _bb;	        //the bounding box of current subtree
        sizeType _tid;			//tag, -1 for internal, index for leave
        sizeType _escape;		//subtree size (number of nodes)
    };
    AABBvhTpl():_eps(1E-3f){}
    AABBvhTpl(const vector<PT3,Eigen::aligned_allocator<PT3> >& verts,const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,bool edge=false,sizeType d=3):_eps(1E-3f){reset(verts,inds,edge,d);}
    AABBvhTpl(const vector<BBox<T> >& bbs,sizeType d=3){reset(bbs,d);}
    virtual bool write(std::ostream& os) const{
        sizeType nr=(sizeType)_internals.size();
        writeBinaryData(nr,os);
        for(sizeType i=0;i<nr;i++)
        {
            writeBinaryData(_internals[i]._bb,os);
            writeBinaryData(_internals[i]._tid,os);
            writeBinaryData(_internals[i]._escape,os);
        }
        return os.good();
    }
    virtual bool read(std::istream& is){
        sizeType nr;
        readBinaryData(nr,is);
        _internals.resize(nr);
        for(sizeType i=0;i<nr;i++)
        {
            readBinaryData(_internals[i]._bb,is);
            readBinaryData(_internals[i]._tid,is);
            readBinaryData(_internals[i]._escape,is);
        }
        return is.good();
    }
    void reset(const vector<PT3,Eigen::aligned_allocator<PT3> >& verts,const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,bool edge,sizeType d=3)
    {
        const sizeType nrTri=(sizeType)inds.size();
        sizeType currIndex=0;

        //set all leaves
        _leaves.resize(nrTri);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrTri; i++) {
            _leaves[i]._bb.setPoints(verts[inds[i][0]],verts[inds[i][1]],verts[inds[i][edge?1:2]]);
            _leaves[i]._bb.enlarged(_eps,d);
            _leaves[i]._tid=i;
            _leaves[i]._escape=1;
        }

        //build tree
        _internals.resize(nrTri*2);
        buildTree(0,nrTri,currIndex);
        ASSERT(!nrTri || (currIndex == nrTri*2-1))

        //force clear
        vector<BVHNode,Eigen::aligned_allocator<BVHNode> > tmp;
        _leaves.swap(tmp);
    }
    void reset(const vector<BBox<T> >& bbs,sizeType d=3)
    {
        const sizeType nrBB=(sizeType)bbs.size();
        sizeType currIndex=0;

        //set all leaves
        _leaves.resize(nrBB);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrBB; i++) {
            _leaves[i]._bb=bbs[i];
            _leaves[i]._bb.enlarged(_eps,d);
            _leaves[i]._tid=i;
            _leaves[i]._escape=1;
        }

        //build tree
        _internals.resize(nrBB*2);
        buildTree(0,nrBB,currIndex);
        ASSERT(!nrBB || (currIndex == nrBB*2-1))

        //force clear
        vector<BVHNode,Eigen::aligned_allocator<BVHNode> > tmp;
        _leaves.swap(tmp);
    }
    virtual ~AABBvhTpl() {}
    //intersect as AABB tree
    void decideRoot(const BBox<T>& bb,sizeType& rootIndex) const
    {
        if(_internals.empty()) {
            rootIndex=-1;
            return;
        }

        sizeType leftId,rightId;
        while(_internals[rootIndex]._tid < 0) {
            leftId=getLeft(rootIndex);
            rightId=getRight(rootIndex);

            if(!bb.intersect(_internals[leftId]._bb))
                rootIndex=rightId;
            else if(!bb.intersect(_internals[rightId]._bb))
                rootIndex=leftId;
            else
                break;
        }
    }
    const vector<BVHNode,Eigen::aligned_allocator<BVHNode> >& internals() const {
        return _internals;
    }
public:	//utility functions
    FORCE_INLINE sizeType getLeft(const sizeType& id) const {
        ASSERT(_internals[id]._tid == -1);
        return id+1;
    }
    FORCE_INLINE sizeType getRight(const sizeType& id) const {
        ASSERT(_internals[id]._tid == -1);
        return id+_internals[id+1]._escape+1;
    }
    FORCE_INLINE sizeType getLeavesOnSubtree(const sizeType& id) const {
        return ((_internals[id]._escape+1) >> 1);
    }
    FORCE_INLINE const BVHNode& getNode(const sizeType& id) const{return _internals[id];}
    template <typename T2> FORCE_INLINE void addAllChild(sizeType id,T2& iter) const {
        if(_internals[id]._tid == -1) {
            addAllChild(getLeft(id),iter);
            addAllChild(getRight(id),iter);
        } else {
            *(iter++)=_internals[id]._tid;
        }
    }
    void writeVTK(const std::string& path) const
    {
        boost::filesystem::create_directory(path);
        std::vector<sizeType> curr(1,0);
        for(sizeType i=0;!curr.empty();i++)
        {
            std::ostringstream oss;oss << path << "/BVH_Level" << i << ".vtk";
            VTKWriter<T> os("BVH Level",oss.str(),true);
            std::vector<PT3,Eigen::aligned_allocator<PT3> > hex;

            std::vector<sizeType> nextLv;
            for(sizeType j=0;j<(sizeType)curr.size();j++)
            {
                hex.push_back(getNode(curr[j])._bb._minC);
                hex.push_back(getNode(curr[j])._bb._maxC);
                if(getNode(curr[j])._escape > 1)
                {
                    nextLv.push_back(getLeft(curr[j]));
                    nextLv.push_back(getRight(curr[j]));
                }
            }
            os.appendVoxels(hex.begin(),hex.end(),true);
            curr.swap(nextLv);
        }
    }
protected:	//tree building
    void buildTree(sizeType from,sizeType to,sizeType& currIndex)
    {
        const sizeType nrTri=to-from;
        if(!nrTri)
            return;

        const sizeType currId=currIndex;
        BVHNode& node=_internals[currIndex];

        if(nrTri == 1) {
            node=_leaves[from];
            currIndex++;
            return;
        } else {
            node._bb.reset();
            node._tid=-1;
        }

        T mean;
        sizeType splitAxis=calcSplitAxisAndBBox(from,to,mean);
        sizeType splitIndex=calcSplitIndex(from,to,splitAxis,mean);

        for(sizeType i=from; i<to; i++)
            node._bb.setUnion(_leaves[i]._bb);

        currIndex++;
        buildTree(from,splitIndex,currIndex);
        buildTree(splitIndex,to,currIndex);

        node._escape=currIndex-currId;
    }
    sizeType calcSplitAxisAndBBox(const sizeType& from,const sizeType& to,T& mean)
    {
        PT3 var(0.0f,0.0f,0.0f);
        PT3 ex(0.0f,0.0f,0.0f);
        PT3 ex2(0.0f,0.0f,0.0f);

        for(sizeType i=from; i<to; i++) {
            const PT3 ctr=(_leaves[i]._bb._minC+_leaves[i]._bb._maxC)*0.5f;
            ex+=ctr;
            ex2+=(ctr.array()*ctr.array()).matrix();
        }
        ex2/=(T)(to-from);
        ex /=(T)(to-from);
        var=ex2-(ex.array()*ex.array()).matrix();

        sizeType maxVar=0;
        if(var(1) > var(maxVar))
            maxVar=1;
        if(var(2) > var(maxVar))
            maxVar=2;

        mean=ex(maxVar);
        return maxVar;
    }
    sizeType calcSplitIndex(const sizeType& from,const sizeType& to,const sizeType& splitAxis,T& mean)
    {
        const sizeType nrTri=to-from;
        sizeType splitIndex=from;
        for(sizeType i=from; i<to; i++) {
            const T ctr=(_leaves[i]._bb._minC(splitAxis)+_leaves[i]._bb._maxC(splitAxis))*0.5f;
            if(ctr > mean) {
                swap(_leaves[i],_leaves[splitIndex]);
                splitIndex++;
            }
        }

        //a simple port of bullet's unbalance adjustment
        sizeType rangeBalancedIndices = nrTri/3;
        bool unbalanced=( (splitIndex <= (from+rangeBalancedIndices)) || (splitIndex >= (to-1-rangeBalancedIndices)) );
        if(unbalanced)
            splitIndex=from+(nrTri>>1);
        bool unbal=(splitIndex==from || splitIndex == to);
        ASSERT(!unbal);

        return splitIndex;
    }
protected:
    //data mem
    vector<BVHNode,Eigen::aligned_allocator<BVHNode> > _leaves;
    vector<BVHNode,Eigen::aligned_allocator<BVHNode> > _internals;
    scalar _eps;
};

template <typename T>
class AABBvhMeshTpl : public AABBvhTpl<T>
{
public:
  typedef typename ScalarUtil<T>::ScalarVec3 PT3; // add by zsw
    AABBvhMeshTpl(const vector<PT3,Eigen::aligned_allocator<PT3> >& verts,const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& inds,bool edge=false)
        :AABBvhTpl<T>(verts,inds,edge),_verts(verts),_inds(inds){}
    bool intersect(const BBox<T>& bb,const sizeType& rootIndex) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return false;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(bb.intersect(currNode._bb)) {
            if(bb.contain(currNode._bb))
                return true;

            if(currNode._tid < 0) {
              return intersect(bb,AABBvhTpl<T>::getLeft(rootIndex)) ||
                  intersect(bb,AABBvhTpl<T>::getRight(rootIndex));
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(bb))
                    return true;
            }
        }

        return false;
    }
    bool intersectLineSeg2D(const PT3& a,const PT3& b,const sizeType& rootIndex) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return false;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.intersect(a,b,2)) {
            if(currNode._tid < 0) {
              return intersectLineSeg2D(a,b,AABBvhTpl<T>::getLeft(rootIndex)) ||
                  intersectLineSeg2D(a,b,AABBvhTpl<T>::getRight(rootIndex));
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                T s,t;
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(LineSegTpl<T>(a,b),s,t))
                    return true;
            }
        }

        return false;
    }
    bool intersectLineSeg3D(const PT3& a,const PT3& b,const sizeType& rootIndex) const
    {
        T t=ScalarUtil<T>::scalar_max;
        return distanceToRay3D(a,b,rootIndex,t,false);
    }
    bool distanceToRay3D(const PT3& a,const PT3& b,const sizeType& rootIndex,T& t,bool inf) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return false;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.intersect(a,b,3)) {
            if(currNode._tid < 0) {
                return distanceToRay3D(a,b,AABBvhTpl<T>::getLeft(rootIndex),t,inf) ||
                       distanceToRay3D(a,b,AABBvhTpl<T>::getRight(rootIndex),t,inf);
            } else {
                const Vec3i& ind=_inds[currNode._tid];
				T tTmp;
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(LineSegTpl<T>(a,b),tTmp,inf) && tTmp < t){
					t=tTmp;
					return true;
				}
            }
        }

        return false;
    }
    void distanceToMesh3D(const PT3& a,const sizeType& rootIndex,sizeType& id,PT3& closest,PT3& bary,scalar& minDist) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.distTo(a,3) < std::abs(minDist)) {
            if(currNode._tid < 0) {
                distanceToMesh3D(a,AABBvhTpl<T>::getLeft(rootIndex),id,closest,bary,minDist);
                distanceToMesh3D(a,AABBvhTpl<T>::getRight(rootIndex),id,closest,bary,minDist);
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                TriangleTpl<T> tri(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]);
                scalar sqrDist;
                PT3 cp,b;
                tri.calcPointDist(a,sqrDist,cp,b);
                sqrDist=sqrt(sqrDist);
                if(sqrDist < std::abs(minDist))
                {
                    id=currNode._tid;
                    closest=cp;
                    bary=b;
                    minDist=sqrDist*(scalar)sgn((a-cp).dot(tri.normal()));
                }
            }
        }
    }
    void distanceToMesh2D(const PT3& a,const sizeType& rootIndex,sizeType& id,PT3& closest,PT3& bary,scalar& minDist) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.distTo(a,2) < std::abs(minDist)) {
            if(currNode._tid < 0) {
                distanceToMesh2D(a,AABBvhTpl<T>::getLeft(rootIndex),id,closest,bary,minDist);
                distanceToMesh2D(a,AABBvhTpl<T>::getRight(rootIndex),id,closest,bary,minDist);
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                LineSegTpl<T> lineSeg(_verts[ind.x()],_verts[ind.y()]);
                scalar sqrDist;
                PT3 cp,b;
                lineSeg.calcPointDist(a,sqrDist,cp,b);
                sqrDist=sqrt(sqrDist);
                if(sqrDist < std::abs(minDist))
                {
                    id=currNode._tid;
                    closest=cp;
                    bary=b;
                    minDist=sqrDist*(scalar)sgn((a-cp).dot(lineSeg.normal()));
                }
            }
        }
    }
    sizeType intersectNr(const BBox<T>& bb,const sizeType& rootIndex) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return 0;

        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(bb.intersect(currNode._bb)) {
            if(bb.contain(currNode._bb))
                return (currNode._escape+1) >> 1;

            if(currNode._tid < 0) {
                return intersectNr(bb,AABBvhTpl<T>::getLeft(rootIndex))+
                       intersectNr(bb,AABBvhTpl<T>::getRight(rootIndex));
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(bb))
                    return 1;
            }
        }
        return 0;
    }
    template <typename T2> bool intersect(const BBox<T>& bb,const sizeType& rootIndex,T2& iter) const {
        if(AABBvhTpl<T>::_internals.empty())
            return false;

        //coarse test
        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(bb.intersect(currNode._bb)) {
            //if(bb.contain(currNode._bb))
            //{
            //	addAllChild(rootIndex,iter);
            //	return true;
            //}

            if(currNode._tid < 0) {
                bool left=intersect(bb,AABBvhTpl<T>::getLeft(rootIndex),iter);
                bool right=intersect(bb,AABBvhTpl<T>::getRight(rootIndex),iter);
                return left||right;
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(bb)) {
                    *(iter++)=currNode._tid;
                    return true;
                }
            }
        }
        return false;
    }
    template <typename T2> sizeType intersectLineSeg2D(const PT3& a,const PT3& b,const sizeType& rootIndex,T2& iter) const {
        if(AABBvhTpl<T>::_internals.empty())
            return 0;

        //coarse test
        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.intersect(a,b,2)) {
            if(currNode._tid < 0) {
                sizeType nrL=intersectLineSeg2D(a,b,AABBvhTpl<T>::getLeft(rootIndex),iter);
                sizeType nrR=intersectLineSeg2D(a,b,AABBvhTpl<T>::getRight(rootIndex),iter);
                return nrL+nrR;
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                T s,t;
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).intersect(LineSegTpl<T>(a,b),s,t)) {
                    *(iter++)=PT2(s,t);
                    return 1;
                }
            }
        }
        return 0;
    }
    sizeType intersectPointTri2D(const PT3& pt,const sizeType& rootIndex) const
    {
        if(AABBvhTpl<T>::_internals.empty())
            return -1;

        //coarse test
        const typename AABBvhTpl<T>::BVHNode& currNode=AABBvhTpl<T>::_internals[rootIndex];
        if(currNode._bb.containDim<2>(pt)) {
            if(currNode._tid < 0) {
                sizeType left=intersectPointTri2D(pt,AABBvhTpl<T>::getLeft(rootIndex));
                if(left >= 0)
                    return left;

                sizeType right=intersectPointTri2D(pt,AABBvhTpl<T>::getRight(rootIndex));
                if(right >= 0)
                    return right;
            } else {
                const Vec3i& ind=_inds[currNode._tid];
                if(TriangleTpl<T>(_verts[ind.x()],_verts[ind.y()],_verts[ind.z()]).isPlanePointInside(pt))
                    return currNode._tid;
            }
        }
        return -1;
    }
protected:
    const vector<PT3,Eigen::aligned_allocator<PT3> >& _verts;
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& _inds;
};

PRJ_END

#endif
