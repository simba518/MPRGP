#ifndef PARTICLE_CD_H
#define PARTICLE_CD_H

#include "ParticleSet.h"
#include "GridBasic.h"

PRJ_BEGIN

//collision detection
template<typename P_TYPE>
class CollisionFunc
{
public:
    virtual void operator()(const P_TYPE& p,const sizeType& id) =0;
};

template<typename VEC3>
struct ExtractPosDirect {
    static FORCE_INLINE const VEC3& extract(const VEC3& p) {
        return p;
    }
};

template<typename P_TYPE>
struct ExtractPosParticle {
    typedef typename ScalarUtil<typename P_TYPE::scalarType>::ScalarVec3 Vec3Type;
    static FORCE_INLINE const Vec3Type& extract(const P_TYPE& p) {
        return p._pos;
    }
};

template <typename T,typename PS_TYPE,typename EXTRACT_POS>
class CollisionInterfaceTpl
{
public:
    typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
    virtual ~CollisionInterfaceTpl() {}
    void reset(const BBox<T>& bb,const Vec3i& nrCell) {
        _gridStartEnd.reset(nrCell,bb,pair<sizeType,sizeType>(0,0),true,16);
        _maxIndex=_gridStartEnd.getNrCell()-Vec3i::Ones();
        if(_gridStartEnd.getDim() == 2)
            _gridLen=min<T>(_gridStartEnd.getCellSize().x(),_gridStartEnd.getCellSize().y());
        else
            _gridLen=_gridStartEnd.getCellSize().minCoeff();
        //INFOV("Grid Len: %f",_gridLen);
    }
    void reset(const BBox<T>& bb,const Vec3Type& cellSz) {
        BBox<T> bbAdj=bb;
        Vec3i nrCell=ceil((Vec3Type)(bb.getExtent().array()/cellSz.array()).matrix());
        bbAdj._maxC=bbAdj._minC+Vec3Type(nrCell.x()*cellSz.x(),nrCell.y()*cellSz.y(),nrCell.z()*cellSz.z());
        reset(bbAdj,nrCell);
    }
    void resize(const sizeType& nrParticle) {
        _pid.resize(nrParticle);
        _pidBK.resize(nrParticle);
        _index.resize(nrParticle);
        _indexBK.resize(nrParticle);
    }
    virtual void prepare(const PS_TYPE& pSet) {
        const sizeType n=pSet.size();
        if(n == 0)
            return;
        if(_pid.size() != n)
            resize(pSet.size());

        calcHash(n,pSet);
        radixSort(n,&(_pid[0]),&(_pidBK[0]),&(_index[0]),&(_indexBK[0]));
        calcStartEnd(n);
    }
    void fill3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
        sizeType nOff=(sizeType)std::ceil(sqrt(radSqr)/_gridLen);
        Vec3i base=floor(_gridStartEnd.getIndexFracSafe(pos));
        base=compMin(compMax(base,Vec3i::Ones()),_maxIndex);

        for(sizeType x=base.x()-nOff; x <= base.x()+nOff; x++)
            for(sizeType y=base.y()-nOff; y <= base.y()+nOff; y++)
                for(sizeType z=base.z()-nOff; z <= base.z()+nOff; z++) {
                    Vec3i id(x,y,z);
                    if(_gridStartEnd.isSafeIndex(id))
                    {
                        const pair<sizeType,sizeType>& range=_gridStartEnd.get(id);
                        for(sizeType k=range.first; k < range.second; k++)
                            if((EXTRACT_POS::extract(pSet[_index[k]])-pos).squaredNorm() < radSqr)
                                f(pSet[_index[k]],_index[k]);
                    }
                }
    }
    void fill2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
        sizeType nOff=(sizeType)std::ceil(sqrt(radSqr)/_gridLen);
        Vec3i base=floor(_gridStartEnd.getIndexFracSafe(pos));
        base=compMin(compMax(base,Vec3i::Ones()),_maxIndex);

        for(sizeType x=base.x()-nOff; x <= base.x()+nOff; x++)
            for(sizeType y=base.y()-nOff; y <= base.y()+nOff; y++) {
                Vec3i id(x,y,0);
                if(_gridStartEnd.isSafeIndex(id))
                {
                    const pair<sizeType,sizeType>& range=_gridStartEnd.get(id);
                    for(sizeType k=range.first; k < range.second; k++)
                        if((EXTRACT_POS::extract(pSet[_index[k]])-pos).squaredNorm() < radSqr)
                            f(pSet[_index[k]],_index[k]);
                }
            }
    }
    void fill3DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
        const sizeType n=pSet.size();
        for(sizeType k=0; k < n; k++)
            if((EXTRACT_POS::extract(pSet[k])-pos).squaredNorm() < radSqr)
                f(pSet.get(k),k);
    }
    void fill2DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
        const sizeType n=pSet.size();
        for(sizeType k=0; k < n; k++)
            if((EXTRACT_POS::extract(pSet[k])-pos).squaredNorm() < radSqr)
                f(pSet.get(k),k);
    }
    bool hasNeigh3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
        sizeType nOff=(sizeType)std::ceil(sqrt(radSqr)/_gridLen);
        Vec3i base=floor(_gridStartEnd.getIndexFracSafe(pos));
        base=compMin(compMax(base,Vec3i::Ones()),_maxIndex);

        for(sizeType x=base.x()-nOff; x <= base.x()+nOff; x++)
            for(sizeType y=base.y()-nOff; y <= base.y()+nOff; y++)
                for(sizeType z=base.z()-nOff; z <= base.z()+nOff; z++) {
                    const pair<sizeType,sizeType>& range=_gridStartEnd.getSafe(Vec3i(x,y,z));
                    for(sizeType k=range.first; k < range.second; k++)
                        if((EXTRACT_POS::extract(pSet[_index[k]])-pos).squaredNorm() < radSqr)
                            return true;
                }
        return false;
    }
    bool hasNeigh2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
        sizeType nOff=(sizeType)std::ceil(sqrt(radSqr)/_gridLen);
        Vec3i base=floor(_gridStartEnd.getIndexFracSafe(pos));
        base=compMin(compMax(base,Vec3i::Ones()),_maxIndex);

        for(sizeType x=base.x()-nOff; x <= base.x()+nOff; x++)
            for(sizeType y=base.y()-nOff; y <= base.y()+nOff; y++) {
                const pair<sizeType,sizeType>& range=_gridStartEnd.getSafe(Vec3i(x,y,0));
                for(sizeType k=range.first; k < range.second; k++)
                    if((EXTRACT_POS::extract(pSet[_index[k]])-pos).squaredNorm() < radSqr)
                        return true;
            }
        return false;
    }
    sizeType getDim() const {
        return _gridStartEnd.getDim();
    }
    BBox<T> getBB() const {
        return _gridStartEnd.getBB();
    }
    Vec3i getNrCell() const {
        return _gridStartEnd.getNrCell();
    }
    const Grid<pair<sizeType,sizeType>,T>& getStartEndGrid() const {
        return _gridStartEnd;
    }
    const std::vector<sizeType,Eigen::aligned_allocator<sizeType> >& getIndex() const {
        return _index;
    }
public:
    static void radixSort(const sizeType& n,sizeType* pid,sizeType* pidBK,sizeType* index,sizeType* indexBK) {
        sizeType bucket[8][256];
        memset(bucket,0,sizeof(sizeType)*256*8);

        //count number of values in a bucket
        for(sizeType i=0; i<n; i++) {
            bucket[0][(pid[i]>>0LL )&255LL]++;
            bucket[1][(pid[i]>>8LL )&255LL]++;
            bucket[2][(pid[i]>>16LL)&255LL]++;
            bucket[3][(pid[i]>>24LL)&255LL]++;
            bucket[4][(pid[i]>>32LL)&255LL]++;
            bucket[5][(pid[i]>>40LL)&255LL]++;
            bucket[6][(pid[i]>>48LL)&255LL]++;
            bucket[7][(pid[i]>>56LL)&255LL]++;
            index[i]=i;
        }

        //accumulate buckets
        for(sizeType i=1; i<256; i++) {
            bucket[0][i]+=bucket[0][i-1];
            bucket[1][i]+=bucket[1][i-1];
            bucket[2][i]+=bucket[2][i-1];
            bucket[3][i]+=bucket[3][i-1];
            bucket[4][i]+=bucket[4][i-1];
            bucket[5][i]+=bucket[5][i-1];
            bucket[6][i]+=bucket[6][i-1];
            bucket[7][i]+=bucket[7][i-1];
        }
        for(sizeType i=255; i>=1; i--) {
            bucket[0][i]=bucket[0][i-1];
            bucket[1][i]=bucket[1][i-1];
            bucket[2][i]=bucket[2][i-1];
            bucket[3][i]=bucket[3][i-1];
            bucket[4][i]=bucket[4][i-1];
            bucket[5][i]=bucket[5][i-1];
            bucket[6][i]=bucket[6][i-1];
            bucket[7][i]=bucket[7][i-1];
        }
        bucket[0][0]=0;
        bucket[1][0]=0;
        bucket[2][0]=0;
        bucket[3][0]=0;
        bucket[4][0]=0;
        bucket[5][0]=0;
        bucket[6][0]=0;
        bucket[7][0]=0;

        //redistribute
        for(sizeType p=0; p<8; p++) {
            sizeType k=p*8;
            for(sizeType i=0; i<n; i++) {
                sizeType& pos=bucket[p][(pid[i]>>(sizeType)k)&255LL];
                indexBK[pos]=index[i];
                pidBK[pos]=pid[i];
                pos++;
            }

            std::swap(pid,pidBK);
            std::swap(index,indexBK);
        }
    }
    void calcHash(const sizeType& n,const PS_TYPE& pSet) {
        const sizeType maxVal=_gridStartEnd.getSzLinear();
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<n; i++) {
            const Vec3Type& pos=EXTRACT_POS::extract(pSet[i]);
            if(_gridStartEnd.getBB().contain(pos)) {
                Vec3i base=floor(_gridStartEnd.getIndexFracSafe(pos));
                _pid[i]=_gridStartEnd.getIndex(base);
            } else {
                _pid[i]=maxVal;
            }
            ASSERT(_pid[i] >= 0 && _pid[i] <= maxVal)
        }
    }
    void calcStartEnd(const sizeType& n) {
        const sizeType maxVal=_gridStartEnd.getSzLinear();
        _gridStartEnd.init(pair<sizeType,sizeType>(0,0));
        OMP_PARALLEL_FOR_
        for(sizeType i=1; i<n; i++) {
            if(_pid[i] != _pid[i-1]) {
                if(_pid[i] != maxVal)
                    _gridStartEnd.get(_pid[i]).first=i;
                _gridStartEnd.get(_pid[i-1]).second=i;
            }
        }

        if(_pid[0] != maxVal)
            _gridStartEnd.get(_pid[0]).first=0;
        if(_pid[n-1] != maxVal)
            _gridStartEnd.get(_pid[n-1]).second=n;
    }
    std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _pid;
    std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _pidBK;
    std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _index;
    std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _indexBK;
    Grid<pair<sizeType,sizeType>,T> _gridStartEnd;
    Vec3i _maxIndex;
    T _gridLen;
};

template <typename T,typename PS_TYPE>
class CollisionInterface : public CollisionInterfaceTpl<T,PS_TYPE,ExtractPosParticle<typename PS_TYPE::ParticleType> > {};

typedef CollisionInterfaceTpl<scalarF,std::vector<Vec3f,Eigen::aligned_allocator<Vec3f> >,ExtractPosDirect<Vec3f> > CollisionInterfaceDirectF;
typedef CollisionInterfaceTpl<scalarD,std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> >,ExtractPosDirect<Vec3d> > CollisionInterfaceDirectD;
typedef CollisionInterfaceTpl<scalar,std::vector<Vec3,Eigen::aligned_allocator<Vec3> >,ExtractPosDirect<Vec3> > CollisionInterfaceDirect;

PRJ_END

#endif
