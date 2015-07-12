#ifndef GRID_OP_H
#define GRID_OP_H

#include "Config.h"
#include "MathBasic.h"
#include "GridBasic.h"
#include "ImplicitFuncInterface.h"

#include "IO.h"
#include "Heap.h"
#include "Zero.h"
#include <algorithm>
#include <set>

PRJ_BEGIN

template <typename T,typename TI,typename TV=vector<T,Eigen::aligned_allocator<T> > >
class GridOp
{
public:
    enum FM_STATE {
        UNKNOWN = 0,
        KNOWN = 1,
        CLOSE = 2,
        POSITIVE = 16,
        NEGATIVE = 32,
    };
    typedef Grid<T,TI,TV> GridType;
    typedef MACGrid<T,TI,TV> MACGridType;
    typedef Grid<unsigned char,TI> TagGridType;
    typedef Grid<sizeType,TI> MarkGridType;
    typedef typename GridType::ValueType ValueType;
    typedef typename GridType::IndexType IndexType;
    //pde
    static void solveEikonalPDE(GridType& from,const T& dist);
    static void smoothEikonalPDE3DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr);
    static void smoothEikonalPDE2DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr);
    static void solveEikonalPDE3DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt);
    static void solveEikonalPDE2DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt);
    static void solveEikonalPDE3D(GridType& from,GridType& to,GridType& smooth,const T& dist);
    static void solveEikonalPDE2D(GridType& from,GridType& to,GridType& smooth,const T& dist);
    static FORCE_INLINE T upwindingGrad(T partialPlus,T partialMinus,const T& phi0);
    //fast march
    template<typename VT2> static bool solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& normals,Grid<IndexType,TI>& normalExtra);
    template<typename VT2> static bool solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes);
    template<typename T2> static bool solveEikonalFM(GridType& from,const T& low,const T& high,Grid<T2,TI>* extra=NULL,T thres=-1.0f) {
        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > known;
        const Vec3i nrPoint=from.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const Vec3i id(x,y,z);
                    if(from.get(id) > low && from.get(id) < high)
                        known.push_back(id);
                }
        return solveEikonalFM(from,known,extra,thres);
    }
    template<typename T2> static bool solveEikonalFM(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=-1.0f);
    template<typename T2> static bool solveEikonalFM3D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=-1.0f);
    template<typename T2> static bool solveEikonalFM2D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=-1.0f);
    template<typename T2> static void extrapolate3D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos);
    template<typename T2> static void extrapolate2D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos);
    static void floodFill3D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags);
    static void floodFill2D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags);
    static bool updateClose2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close);
    static bool updateClose3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close);
    static bool updateKnown2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id);
    static bool updateKnown3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id);
    //reinitialize
    static void reinitialize(GridType& from);
    //smooth
    static void smooth(GridType& from,const TagGridType* tagField=0,unsigned char tag=0);
    static void smooth3D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag);
    static void smooth2D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag);
    //io
    static void write2DScalarGridVTK(const string& path,const GridType& grd,bool moveVertex=true);
    static void write2DScalarBarChartVTK(const string& path,const GridType& grd,bool moveVertex=true,const TI& scale=1.0f);
    static void write2DMarkGridVTK(const string& path,const MarkGridType& grd,bool moveVertex=true);
    static void write2DScalarGridGradVTK(const string& path,const GridType& grd,const sizeType& sampleInterval,const TI& len);
    static void write2DScalarGridGradMagVTK(const string& path,const GridType& grd);
    static void write3DScalarGridVTK(const string& path,const GridType& grd);
    static void write2DVectorGridVTK(const string& path,const Grid<IndexType,TI>& vel,const TI& len=0.0f);
    static void write2DMACGridVTK(const string& path,const MACGridType& mac);
    static void write2DVelocityGridVTK(const string& path,const MACGridType& vel,const TI& len=0.0f);
    static void write3DVelocityGridVTK(const string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter);
    static void advectRK2(const IndexType& from,const TI& time,const MACGridType& vel,
                          std::vector<IndexType,Eigen::aligned_allocator<IndexType> >& vertices,
                          std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> >& indices,
                          std::vector<TI>& colorData);
    //miscellaneous
    template<typename T2>
    static void copyVelFromFunc(MACGridType& vel,const VelFunc<T2>& func);
    template<typename T2>
    static void copyVelFromFunc(ScalarField& vel,const VelFunc<T2>& func,const sizeType& a);
    template<typename T2>
    static void copyFromImplictFunc(GridType& to,const ImplicitFunc<T2>& func);
    template<typename T2,typename TI2>
    static void copyFromOtherGrid(GridType& to,const Grid<T2,TI2>& from);
    template<typename T2,typename TI2>
    static void copyFromOtherGridOfSameGeometry(GridType& to,const Grid<T2,TI2>& from);
    template<typename T2,typename TI2>
    static void copyFromOtherGridOfSameGeometry(MACGridType& to,const MACGrid<T2,TI2>& from);
    template <typename COMPARE,typename ALLOC>
    static void getAllDifferentValues(const GridType& g,std::set<T,COMPARE,ALLOC>& vals);
    template <typename COMPARE,typename ALLOC>
    static void getAllDifferentValues(const MACGridType& g,std::set<T,COMPARE,ALLOC>& vals);
    //interpolation
    static void fromFaceToCenter(const MACGridType& from,Grid<ValueType,T>& to);
    static void fromFaceToCenter3D(const MACGridType& from,Grid<ValueType,T>& to);
    static void fromFaceToCenter2D(const MACGridType& from,Grid<ValueType,T>& to);
    static void fromCenterToFace(const Grid<ValueType,T>& from,MACGridType& to);
    static void fromCenterToFace3D(const Grid<ValueType,T>& from,MACGridType& to);
    static void fromCenterToFace2D(const Grid<ValueType,T>& from,MACGridType& to);
};

//pde
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE(GridType& from,const T& dist)
{
    GridType tmp1=from;
    GridType tmp2=from;
    if(from.getDim() == 3)
        solveEikonalPDE3D(from,tmp1,tmp2,dist);
    else if(from.getDim() == 2)
        solveEikonalPDE2D(from,tmp1,tmp2,dist);
    else ASSERT(false);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smoothEikonalPDE3DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr)
{
    for(sizeType y=0; y<nrPoint.y(); y++)
        for(sizeType z=0; z<nrPoint.z(); z++) {
            T ctr=from.get(Vec3i(x,y,z));
            smooth.get(Vec3i(x,y,z))=ctr/sqrt(ctr*ctr+maxCellSzSqr);
        }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smoothEikonalPDE2DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr)
{
    for(sizeType y=0; y<nrPoint.y(); y++) {
        T ctr=from.get(Vec3i(x,y,0));
        smooth.get(Vec3i(x,y,0))=ctr/sqrt(ctr*ctr+maxCellSzSqr);
    }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE3DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt)
{
    for(sizeType y=0; y<nrPoint.y(); y++)
        for(sizeType z=0; z<nrPoint.z(); z++) {
            T ctr=from.get(Vec3i(x,y,z));
            T smoothVal=smooth.get(Vec3i(x,y,z));
            ValueType upNorm=(ValueType(upwindingGrad(from.getSafe(Vec3i(x+1,y,z))-ctr,ctr-from.getSafe(Vec3i(x-1,y,z)),smoothVal),
                                        upwindingGrad(from.getSafe(Vec3i(x,y+1,z))-ctr,ctr-from.getSafe(Vec3i(x,y-1,z)),smoothVal),
                                        upwindingGrad(from.getSafe(Vec3i(x,y,z+1))-ctr,ctr-from.getSafe(Vec3i(x,y,z-1)),smoothVal)).array()*invCellSz.array()).matrix();
            to.get(Vec3i(x,y,z))=ctr-dt*smoothVal*(upNorm.norm()-1.0f);
        }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE2DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt)
{
    for(sizeType y=0; y<nrPoint.y(); y++) {
        T ctr=from.get(Vec3i(x,y,0));
        T smoothVal=smooth.get(Vec3i(x,y,0));
        ValueType upNorm=(ValueType(upwindingGrad(from.getSafe(Vec3i(x+1,y,0))-ctr,ctr-from.getSafe(Vec3i(x-1,y,0)),smoothVal),
                                    upwindingGrad(from.getSafe(Vec3i(x,y+1,0))-ctr,ctr-from.getSafe(Vec3i(x,y-1,0)),smoothVal),0.0f).array()*invCellSz.array()).matrix();
        upNorm.z()=0.0f;
        to.get(Vec3i(x,y,0))=ctr-dt*smoothVal*(upNorm.norm()-1.0f);
    }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE3D(GridType& from,GridType& to,GridType& smooth,const T& dist)
{
    const T dt=from.getCellSize().minCoeff()*0.3f;
    const IndexType invCellSz=from.getInvCellSize();
    const T maxCellSzSqr=pow(from.getCellSize().maxCoeff(),(TI)2.0f);

    const sizeType nrIter=(sizeType)std::ceil(dist/dt);
    const Vec3i nrPoint=from.getNrPoint();

    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
        smoothEikonalPDE3DLayer(x,from,smooth,nrPoint,maxCellSzSqr);

    for(sizeType iter=0; iter<nrIter; iter++) {
        INFOV("%lu",iter)
        OMP_PARALLEL_FOR_
        for(sizeType x=0; x<nrPoint.x(); x++)
            solveEikonalPDE3DLayer(x,from,smooth,to,nrPoint,invCellSz,dt);
        from.swap(to);
    }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE2D(GridType& from,GridType& to,GridType& smooth,const T& dist)
{
    const T dt=min<TI>(from.getCellSize().x(),from.getCellSize().y())*0.3f;
    const IndexType invCellSz=from.getInvCellSize();
    const T maxCellSzSqr=pow(max<TI>(from.getCellSize().x(),from.getCellSize().y()),(TI)2.0f);

    const sizeType nrIter=(sizeType)std::ceil(dist/dt);
    const Vec3i nrPoint=from.getNrPoint();

    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
        smoothEikonalPDE2DLayer(x,from,smooth,nrPoint,maxCellSzSqr);

    for(sizeType iter=0; iter<nrIter; iter++) {
        OMP_PARALLEL_FOR_
        for(sizeType x=0; x<nrPoint.x(); x++)
            solveEikonalPDE2DLayer(x,from,smooth,to,nrPoint,invCellSz,dt);
        from.swap(to);
    }
}

template <typename T,typename TI,typename TV>
T GridOp<T,TI,TV>::upwindingGrad(T partialPlus,T partialMinus,const T& phi0)
{
    T determinePlus=partialPlus*phi0;
    T determineMinus=partialMinus*phi0;

    if(determineMinus <= 0 && determinePlus <= 0) {} //no modify
    else if(determineMinus >= 0 && determinePlus >= 0)
        partialPlus=partialMinus;
    else if(determineMinus <= 0 && determinePlus >= 0)
        partialPlus=0.0;
    else if(determineMinus >= 0 && determinePlus <= 0) {
        if(std::abs(determinePlus) >= std::abs(determineMinus)) {}//no modify
        else partialPlus=partialMinus;
    }
    return partialPlus;
}

//fast march
template <typename T,typename TI,typename TV>
template<typename VT2>
bool GridOp<T,TI,TV>::solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& normals,Grid<IndexType,TI>& normalExtra)
{
#define CHECK_NEIGH_EXTRA																		\
{																								\
IndexType dir=from.getPt(nPos)-IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z());				\
TI dist=std::abs(dir.dot(IndexType(normals[i].x(),normals[i].y(),normals[i].z())));					\
if(from.get(nPos) == -1.0f){																	\
	knowns.push_back(nPos);																		\
	from.get(nPos)=dist;																		\
    normalExtra.get(nPos)=IndexType(normals[i].x(),normals[i].y(),normals[i].z());				\
}else{																							\
	from.get(nPos)=std::min<TI>(dist,from.get(nPos));											\
    normalExtra.get(nPos)=IndexType(normals[i].x(),normals[i].y(),normals[i].z());				\
}																								\
}

    from.init(-1.0f);

    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > knowns;
    for(sizeType i=0; i<(sizeType)nodes.size(); i++) {
        Vec3i base=floor(from.getIndexFracSafe(IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())));
        Vec3i nPos;
        if(from.getDim() == 2) {
            nPos=base;
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,0,0);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,1,0);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(0,1,0);
            CHECK_NEIGH_EXTRA
        } else if(from.getDim() == 3) {
            nPos=base;
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,0,0);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,1,0);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(0,1,0);
            CHECK_NEIGH_EXTRA

            nPos=base+Vec3i(0,0,1);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,0,1);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(1,1,1);
            CHECK_NEIGH_EXTRA
            nPos=base+Vec3i(0,1,1);
            CHECK_NEIGH_EXTRA
        }
    }

    return solveEikonalFM<IndexType>(from,knowns,NULL);
#undef CHECK_NEIGH_EXTRA
}

template <typename T,typename TI,typename TV>
template<typename VT2>
bool GridOp<T,TI,TV>::solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes)
{
#define CHECK_NEIGH																				\
{																								\
TI dist=(from.getPt(nPos)-IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())).norm();			\
if(from.get(nPos) == -1.0f){																	\
	knowns.push_back(nPos);																		\
	from.get(nPos)=dist;																		\
}else{																							\
	from.get(nPos)=std::min<T>(dist,from.get(nPos));											\
}																								\
}

    from.init(-1.0f);

    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > knowns;
    for(sizeType i=0; i<(sizeType)nodes.size(); i++) {
        Vec3i base=floor(from.getIndexFracSafe(IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())));
        Vec3i nPos;
        if(from.getDim() == 2) {
            nPos=base;
            CHECK_NEIGH
            nPos=base+Vec3i(1,0,0);
            CHECK_NEIGH
            nPos=base+Vec3i(1,1,0);
            CHECK_NEIGH
            nPos=base+Vec3i(0,1,0);
            CHECK_NEIGH
        } else if(from.getDim() == 3) {
            nPos=base;
            CHECK_NEIGH
            nPos=base+Vec3i(1,0,0);
            CHECK_NEIGH
            nPos=base+Vec3i(1,1,0);
            CHECK_NEIGH
            nPos=base+Vec3i(0,1,0);
            CHECK_NEIGH

            nPos=base+Vec3i(0,0,1);
            CHECK_NEIGH
            nPos=base+Vec3i(1,0,1);
            CHECK_NEIGH
            nPos=base+Vec3i(1,1,1);
            CHECK_NEIGH
            nPos=base+Vec3i(0,1,1);
            CHECK_NEIGH
        }
    }

    return solveEikonalFM<IndexType>(from,knowns,NULL);
#undef CHECK_NEIGH
}

template <typename T,typename TI,typename TV>
template <typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres)
{
    if(from.getDim() == 3)
        return solveEikonalFM3D<T2>(from,known,extra,thres);
    else if(from.getDim() == 2)
        return solveEikonalFM2D<T2>(from,known,extra,thres);
    else ASSERT(false);
    return false;
}

template <typename T,typename TI,typename TV>
template<typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM3D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres)
{
    //tags and marks
    TagGridType tags;
    tags.makeSameGeometry(from);
    tags.init(UNKNOWN);

    MarkGridType heapOffsets;
    heapOffsets.makeSameGeometry(from);
    heapOffsets.init(-1);

    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > closes;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > heap;

    //initialize heap
    const sizeType nrKnown=(sizeType)known.size();
    for(sizeType i=0; i<nrKnown; i++)
        tags.get(known[i])=KNOWN;
    for(sizeType i=0; i<nrKnown; i++) {
        const Vec3i& id=known[i];
        if(id.x() > 0) {
            Vec3i nPos=Vec3i(id.x()-1,id.y(),id.z());
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.x() < tags.getNrPoint().x()-1) {
            Vec3i nPos=Vec3i(id.x()+1,id.y(),id.z());
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.y() > 0) {
            Vec3i nPos=Vec3i(id.x(),id.y()-1,id.z());
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.y() < tags.getNrPoint().y()-1) {
            Vec3i nPos=Vec3i(id.x(),id.y()+1,id.z());
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.z() > 0) {
            Vec3i nPos=Vec3i(id.x(),id.y(),id.z()-1);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.z() < tags.getNrPoint().z()-1) {
            Vec3i nPos=Vec3i(id.x(),id.y(),id.z()+1);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
    }

    //add all close to heap
    const sizeType nrClose=closes.size();
    for(sizeType i=0; i<nrClose; i++)
        if(!updateClose3D(from,tags,heapOffsets,heap,closes[i]))
            return false;

    //march until heap empty
    while(!heap.empty()) {
        Vec3i pos=popHeapAbs(from,heapOffsets,heap,Vec3i::Constant(-1));
        tags.get(pos)=KNOWN;

        if(thres > 0.0f && std::abs(from.get(pos)) > thres)
            break;

        if(extra)
            extrapolate3D(from,tags,*extra,pos);

        if(!updateKnown3D(from,tags,heapOffsets,heap,pos))
            return false;
    }

    if(!heap.empty())
        floodFill3D(from,heap,tags);
    return true;
}

template <typename T,typename TI,typename TV>
template<typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM2D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres)
{
    //tags and marks
    TagGridType tags;
    tags.makeSameGeometry(from);
    tags.init(UNKNOWN);

    MarkGridType heapOffsets;
    heapOffsets.makeSameGeometry(from);
    heapOffsets.init(-1);

    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > closes;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > heap;

    //initialize heap
    const sizeType nrKnown=(sizeType)known.size();
    for(sizeType i=0; i<nrKnown; i++)
        tags.get(known[i])=KNOWN;
    for(sizeType i=0; i<nrKnown; i++) {
        const Vec3i& id=known[i];
        if(id.x() > 0) {
            Vec3i nPos=Vec3i(id.x()-1,id.y(),0);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.x() < tags.getNrPoint().x()-1) {
            Vec3i nPos=Vec3i(id.x()+1,id.y(),0);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.y() > 0) {
            Vec3i nPos=Vec3i(id.x(),id.y()-1,0);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
        if(id.y() < tags.getNrPoint().y()-1) {
            Vec3i nPos=Vec3i(id.x(),id.y()+1,0);
            unsigned char& tag=tags.get(nPos);
            if(tag == UNKNOWN) {
                tag=CLOSE;
                closes.push_back(nPos);
            }
        }
    }

    //add all close to heap
    const sizeType nrClose=closes.size();
    for(sizeType i=0; i<nrClose; i++)
        if(!updateClose2D(from,tags,heapOffsets,heap,closes[i]))
            return false;

    //march until heap empty
    while(!heap.empty()) {
        Vec3i pos=popHeapAbs(from,heapOffsets,heap,Vec3i::Constant(-1));
        tags.get(pos)=KNOWN;

        if(thres > 0.0f && std::abs(from.get(pos)) > thres)
            break;

        if(extra)
            extrapolate2D(from,tags,*extra,pos);

        if(!updateKnown2D(from,tags,heapOffsets,heap,pos))
            return false;
    }

    if(!heap.empty())
        floodFill2D(from,heap,tags);
    return true;
}

template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::extrapolate3D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos)
{
#define CHECK_DIR_EXTRA(DIR,N,P)									\
has=false;															\
if(pos.DIR() > 0 && tags.get(N) == KNOWN)							\
{																	\
	extraDenomTmp=std::abs(from.get(pos)-from.get(N))*					\
				  invCellSz.DIR();									\
	extraNumTmp=extraRef.get(N)*extraDenomTmp;						\
	has=true;														\
}																	\
if(pos.DIR() < nrPoint.DIR()-1 && tags.get(P) == KNOWN)				\
{																	\
	if(!has || std::abs(from.get(P)) < std::abs(from.get(N)))					\
	{																\
		extraDenomTmp=std::abs(from.get(pos)-from.get(P))*				\
					  invCellSz.DIR();								\
		extraNumTmp=extraRef.get(P)*extraDenomTmp;					\
	}																\
	has=true;														\
}																	\
if(has)																\
{																	\
	extraNum+=extraNumTmp;											\
	extraDenom+=extraDenomTmp;										\
}

    Vec3i nrPoint=from.getNrPoint();
    IndexType invCellSz=from.getInvCellSize();
    Grid<T2,TI>& extraRef=extra;

    T2 extraNum=Zero<T2>::value();
    T extraDenom=0.0f;
    bool has;

    Vec3i NX=pos-Vec3i::Unit(0);
    Vec3i PX=pos+Vec3i::Unit(0);
    Vec3i NY=pos-Vec3i::Unit(1);
    Vec3i PY=pos+Vec3i::Unit(1);
    Vec3i NZ=pos-Vec3i::Unit(2);
    Vec3i PZ=pos+Vec3i::Unit(2);

    T2 extraNumTmp=Zero<T2>::value();
    T extraDenomTmp=0.0f;
    CHECK_DIR_EXTRA(x,NX,PX)
    CHECK_DIR_EXTRA(y,NY,PY)
    CHECK_DIR_EXTRA(z,NZ,PZ)
    extraRef.get(pos)=extraNum/std::max<T>(extraDenom,ScalarUtil<T>::scalar_eps);
}

template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::extrapolate2D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos)
{
    Vec3i nrPoint=from.getNrPoint();
    IndexType invCellSz=from.getInvCellSize();
    Grid<T2,TI>& extraRef=extra;

    T2 extraNum=Zero<T2>::value();
    T extraDenom=0.0f;
    bool has;

    Vec3i NX=pos-Vec3i::Unit(0);
    Vec3i PX=pos+Vec3i::Unit(0);
    Vec3i NY=pos-Vec3i::Unit(1);
    Vec3i PY=pos+Vec3i::Unit(1);

    T2 extraNumTmp=Zero<T2>::value();
    T extraDenomTmp=0.0f;
    CHECK_DIR_EXTRA(x,NX,PX)
    CHECK_DIR_EXTRA(y,NY,PY)
    extraRef.get(pos)=extraNum/std::max<T>(extraDenom,ScalarUtil<T>::scalar_eps);
#undef CHECK_DIR_EXTRA
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::floodFill2D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags)
{
#define CHECK_DIR(dir)				\
if(tags.isSafeIndex(dir))			\
{									\
	if(tags.get(dir) == KNOWN)		\
	{tags.get(pos)=KNOWN;			\
	 from.get(pos)=from.get(dir);}	\
	else{id.push_back(dir);}		\
}

    while(!id.empty()) {
        const Vec3i pos=id.back();
        id.pop_back();
        CHECK_DIR(pos+Vec3i(-1,0,0))
        CHECK_DIR(pos+Vec3i(1,0,0))
        CHECK_DIR(pos+Vec3i(0,-1,0))
        CHECK_DIR(pos+Vec3i(0,1,0))
    }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::floodFill3D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags)
{
    while(!id.empty()) {
        const Vec3i pos=id.back();
        id.pop_back();
        CHECK_DIR(pos+Vec3i(-1,0,0))
        CHECK_DIR(pos+Vec3i(1,0,0))
        CHECK_DIR(pos+Vec3i(0,-1,0))
        CHECK_DIR(pos+Vec3i(0,1,0))
        CHECK_DIR(pos+Vec3i(0,0,-1))
        CHECK_DIR(pos+Vec3i(0,0,1))
    }

#undef CHECK_DIR
}

template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateClose2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close)
{
#define CHECK_DIR(DIR,N,P)											\
exist=false;														\
if(close.DIR() > 0 && tags.get(N) == KNOWN){						\
	B[off]=result.get(N);											\
	exist=true;														\
}																	\
if(close.DIR() < tags.getNrPoint().DIR()-1 && tags.get(P) == KNOWN)	\
{																	\
	tmp=result.get(P);												\
	if(exist){														\
		if(std::abs(tmp) < std::abs(B[off]))									\
			B[off]=tmp;												\
	}																\
	else{															\
		B[off]=tmp;													\
	}																\
	exist=true;														\
}																	\
if(exist)															\
	off++;

    //locals
    T B[3];
    sizeType off=0;
    T a,b,c,tmp,det;
    bool exist;
    typename GridType::ValueType invCellSz=tags.getInvCellSize();
    typename GridType::ValueType invCellSzSqr=(invCellSz.array()*invCellSz.array()).matrix();
    Vec3i NX=close-Vec3i::Unit(0);
    Vec3i PX=close+Vec3i::Unit(0);
    Vec3i NY=close-Vec3i::Unit(1);
    Vec3i PY=close+Vec3i::Unit(1);

    //check every direction
    CHECK_DIR(x,NX,PX)
    CHECK_DIR(y,NY,PY)

    //check same sign
    bool positive=B[0] > 0.0f;
    for(sizeType i=1; i<off; i++) {
        if(B[0]*B[i]<0.0f) {
            WARNING("In-Out Test Fail For Fast Marching,Mesh Not Watertight!")
            return false;
        }
    }

    //sort the minimal direction
    if(off == 2) {
        if(std::abs(B[1]) < std::abs(B[0]))
            swap(B[0],B[1]);
    }

    //try solving
    while(off>0) {
        //build a,b,c
        a=b=0.0f;
        c=-1.0f;
        for(sizeType i=0; i<off; i++) {
            a+=invCellSzSqr[i];
            b-=2.0f*B[i]*invCellSzSqr[i];
            c+=B[i]*B[i]*invCellSzSqr[i];
        }

        det=b*b-4.0f*a*c;
        if(det < 0.0f) {
            off--;
        } else {
            //solution to the quadratic
            a=max<T>(2.0f*a,EPS);
            result.get(close)=positive ? ((-b+sqrt(det))/a) : ((-b-sqrt(det))/a);
            result.get(close)=std::abs(result.get(close));
            if(!positive)
                result.get(close)*=-1.0f;

            //update heap
            if(heapOffsets.get(close) < 0)
                pushHeapAbs(result,heapOffsets,heap,close);
            else
                updateHeapAbs(result,heapOffsets,heap,close);

            break;
        }
    }

    if(off == 0) {
        WARNING("Numerically Impossible!")
        return false;
    }
    return true;
}

template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateClose3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close)
{
    //locals
    T B[3];
    sizeType off=0;
    T a,b,c,tmp,det;
    bool exist;
    typename GridType::ValueType invCellSz=tags.getInvCellSize();
    typename GridType::ValueType invCellSzSqr=(invCellSz.array()*invCellSz.array()).matrix();
    Vec3i NX=close-Vec3i::Unit(0);
    Vec3i PX=close+Vec3i::Unit(0);
    Vec3i NY=close-Vec3i::Unit(1);
    Vec3i PY=close+Vec3i::Unit(1);
    Vec3i NZ=close-Vec3i::Unit(2);
    Vec3i PZ=close+Vec3i::Unit(2);

    //check every direction
    CHECK_DIR(x,NX,PX)
    CHECK_DIR(y,NY,PY)
    CHECK_DIR(z,NZ,PZ)

    //check same sign
    bool positive=B[0] > 0.0f;
    for(sizeType i=1; i<off; i++) {
        if(B[0]*B[i]<0.0f) {
            WARNING("In-Out Test Fail For Fast Marching,Mesh Not Watertight!")
            return false;
        }
    }

    //sort the minimal direction
    if(off == 3) {
        if(std::abs(B[1]) < std::abs(B[0]))
            swap(B[0],B[1]);

        if(std::abs(B[2]) < std::abs(B[1]))
            swap(B[1],B[2]);

        if(std::abs(B[1]) < std::abs(B[0]))
            swap(B[0],B[1]);
    } else if(off == 2) {
        if(std::abs(B[1]) < std::abs(B[0]))
            swap(B[0],B[1]);
    }

    //try solving
    while(off>0) {
        //build a,b,c
        a=b=0.0f;
        c=-1.0f;
        for(sizeType i=0; i<off; i++) {
            a+=invCellSzSqr[i];
            b-=2.0f*B[i]*invCellSzSqr[i];
            c+=B[i]*B[i]*invCellSzSqr[i];
        }

        det=b*b-4.0f*a*c;
        if(det < 0.0f) {
            off--;
        } else {
            //solution to the quadratic
            a=max<T>(2.0f*a,EPS);
            result.get(close)=positive ? ((-b+sqrt(det))/a) : ((-b-sqrt(det))/a);
            result.get(close)=std::abs(result.get(close));
            if(!positive)
                result.get(close)*=-1.0f;

            //update heap
            if(heapOffsets.get(close) < 0)
                pushHeapAbs(result,heapOffsets,heap,close);
            else
                updateHeapAbs(result,heapOffsets,heap,close);

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

template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateKnown2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id)
{
    bool ret=true;
    if(id.x() > 0) {
        Vec3i nPos=Vec3i(id.x()-1,id.y(),0);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.x() < tags.getNrPoint().x()-1) {
        Vec3i nPos=Vec3i(id.x()+1,id.y(),0);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.y() > 0) {
        Vec3i nPos=Vec3i(id.x(),id.y()-1,0);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.y() < tags.getNrPoint().y()-1) {
        Vec3i nPos=Vec3i(id.x(),id.y()+1,0);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
        }
    }
    return ret;
}

template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateKnown3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id)
{
    bool ret=true;
    if(id.x() > 0) {
        Vec3i nPos=Vec3i(id.x()-1,id.y(),id.z());
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.x() < tags.getNrPoint().x()-1) {
        Vec3i nPos=Vec3i(id.x()+1,id.y(),id.z());
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.y() > 0) {
        Vec3i nPos=Vec3i(id.x(),id.y()-1,id.z());
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.y() < tags.getNrPoint().y()-1) {
        Vec3i nPos=Vec3i(id.x(),id.y()+1,id.z());
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.z() > 0) {
        Vec3i nPos=Vec3i(id.x(),id.y(),id.z()-1);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    if(id.z() < tags.getNrPoint().z()-1) {
        Vec3i nPos=Vec3i(id.x(),id.y(),id.z()+1);
        unsigned char& tag=tags.get(nPos);
        if(tag != KNOWN) {
            tag=CLOSE;
            ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
        }
    }
    return ret;
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::reinitialize(GridType& from)
{
    const T thres=from.getCellSize().maxCoeff()*2.0f;
    const Vec3i nrP=from.getNrPoint();
    solveEikonalPDE(from,thres*3.0f);
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > known;
    if(from.getDim() == 3) {
        for(sizeType x=0; x<nrP.x(); x++)
            for(sizeType y=0; y<nrP.y(); y++)
                for(sizeType z=0; z<nrP.z(); z++) {
                    const T val=from.get(Vec3i(x,y,z));
                    if(std::abs(val)<thres)
                        known.push_back(Vec3i(x,y,z));
                }
    } else {
        for(sizeType x=0; x<nrP.x(); x++)
            for(sizeType y=0; y<nrP.y(); y++) {
                const T val=from.get(Vec3i(x,y,0));
                if(std::abs(val)<thres)
                    known.push_back(Vec3i(x,y,0));
            }
    }
    solveEikonalFM<IndexType>(from,known);
}

//smooth
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth(GridType& from,const TagGridType* tagField,unsigned char tag)
{
    GridType tmp=from;
    if(from.getDim() == 3)
        smooth3D(from,tmp,tagField,tag);
    else if(from.getDim() == 2)
        smooth2D(from,tmp,tagField,tag);
    else ASSERT(false);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth3D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag)
{
    const scalar weight=1.0f;//0.6f;
    const Vec3i nrPoint=from.getNrPoint();

    OMP_PARALLEL_FOR_
    for(sizeType x=1; x<nrPoint.x()-1; x++)
        for(sizeType y=1; y<nrPoint.y()-1; y++)
            for(sizeType z=1; z<nrPoint.z()-1; z++) {
                if(!tagField || tagField->get(Vec3i(x,y,z)) == tag)
                    tmp.get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z))*(1.0f-weight)+
                                          (from.get(Vec3i(x-1,y,z))+from.get(Vec3i(x+1,y,z))+
                                           from.get(Vec3i(x,y-1,z))+from.get(Vec3i(x,y+1,z))+
                                           from.get(Vec3i(x,y,z-1))+from.get(Vec3i(x,y,z+1)))*weight/6.0f;
            }
    from.swap(tmp);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth2D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag)
{
    const scalar weight=1.0f;//0.6f;
    const Vec3i nrPoint=from.getNrPoint();

    OMP_PARALLEL_FOR_
    for(sizeType x=1; x<nrPoint.x()-1; x++)
        for(sizeType y=1; y<nrPoint.y()-1; y++) {
            if(!tagField || tagField->get(Vec3i(x,y,0)) == tag)
                tmp.get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0))*(1.0f-weight)+
                                      (from.get(Vec3i(x-1,y,0))+from.get(Vec3i(x+1,y,0))+
                                       from.get(Vec3i(x,y-1,0))+from.get(Vec3i(x,y+1,0)))*weight/4.0f;
        }
    from.swap(tmp);
}

//io
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridVTK(const string& path,const GridType& grd,bool moveVertex)
{
    ASSERT(grd.getDim() == 2)
    std::vector<IndexType> pos;
    std::vector<Vec3i> index;

    std::vector<T> pointData;
    const Vec3i nrPoint=grd.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            IndexType p=grd.getPt(Vec3i(x,y,0));
            if(moveVertex)
                p.z()=(TI)grd.get(Vec3i(x,y,0));
            pos.push_back(p);
            pointData.push_back(grd.get(Vec3i(x,y,0)));
        }

    for(sizeType x=0; x<nrPoint.x()-1; x++)
        for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
            index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
            index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
        }

    VTKWriter<TI> writer("2D Grid",path,true);
    writer.appendPoints(pos.begin(),pos.end());
    writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
    writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarBarChartVTK(const string& path,const GridType& grd,bool moveVertex,const TI& scale)
{
    //typedef Eigen::Matrix<sizeType,8,1> IDS;

    ASSERT(grd.getDim() == 2)
    IndexType cellSize=grd.getCellSize()*0.5f*scale;
    cellSize.z()=0.0f;

    std::vector<IndexType> pos;
    std::vector<TI> cellData;

    const Vec3i nrPoint=grd.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            T val=grd.get(Vec3i(x,y,0));
            IndexType zOff(0.0f,0.0f,(TI)val);
            if(!moveVertex)
                zOff.z()=1E-3f;
            IndexType p=grd.getPt(Vec3i(x,y,0));

            pos.push_back(p-cellSize);
            pos.push_back(p+cellSize+zOff);

            cellData.push_back((TI)val);
        }

    VTKWriter<TI> writer("2D Grid BarChart",path,true);
    writer.appendVoxels(pos.begin(),pos.end(),true);
    writer.appendCustomData("value",cellData.begin(),cellData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DMarkGridVTK(const string& path,const MarkGridType& grd,bool moveVertex)
{
    ASSERT(grd.getDim() == 2)
    std::vector<IndexType> pos;
    std::vector<Vec3i> index;

    std::vector<T> pointData;
    const Vec3i nrPoint=grd.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            IndexType p=grd.getPt(Vec3i(x,y,0));
            if(moveVertex)
                p.z()=(T)grd.get(Vec3i(x,y,0));
            pos.push_back(p);
            pointData.push_back((T)grd.get(Vec3i(x,y,0)));
        }

    for(sizeType x=0; x<nrPoint.x()-1; x++)
        for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
            index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
            index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
        }

    VTKWriter<TI> writer("2D Grid",path,true);
    writer.appendPoints(pos.begin(),pos.end());
    writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
    writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridGradVTK(const string& path,const GridType& grd,const sizeType& sampleInterval,const TI& len)
{
    ASSERT(grd.getDim() == 2)
    std::vector<IndexType> pos;
    std::vector<Vec3i> index;

    std::vector<TI> cellData;
    const Vec3i nrPoint=grd.nrPoint();
    for(sizeType x=0; x<nrPoint.x(); x+=sampleInterval)
        for(sizeType y=0; y<nrPoint.y(); y+=sampleInterval) {
            IndexType p=grd.getPt(Vec3i(x,y,0));
            p.z()=0.0f;

            ValueType delta=grd.sampleSafeGrad(p);
            if(delta.norm() < EPS)
                continue;

            delta.normalize();
            delta*=len;
            IndexType end(p.x()+delta.x(),p.y()+delta.y(),p.z()+delta.z());
            end.z()=0.0f;

            index.push_back(Vec3i(pos.size(),pos.size()+1,0));
            cellData.push_back(grd.get(Vec3i(x,y,0)));
            pos.push_back(p);
            pos.push_back(end);
        }

    VTKWriter<TI> writer("2D Grid Grad",path,false);
    writer.appendPoints(pos.begin(),pos.end());
    writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::LINE);
    writer.appendCustomData("value",cellData.begin(),cellData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridGradMagVTK(const string& path,const GridType& grd)
{
    ASSERT(grd.getDim() == 2)
    std::vector<IndexType> pos;
    std::vector<Vec3i> index;

    std::vector<T> pointData;
    const Vec3i nrPoint=grd.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            IndexType p=grd.getPt(Vec3i(x,y,0));
            pos.push_back(p);
            pointData.push_back(grd.sampleSafeGrad(grd.getPt(Vec3i(x,y,0))).norm());
        }

    for(sizeType x=0; x<nrPoint.x()-1; x++)
        for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
            index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
            index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
        }

    VTKWriter<TI> writer("2D Grid",path,true);
    writer.appendPoints(pos.begin(),pos.end());
    writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
    writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write3DScalarGridVTK(const string& path,const GridType& grd)
{
    ASSERT(grd.getDim() == 3)

    std::vector<T> pointData;
    const Vec3i nrPoint=grd.getNrPoint();
    for(sizeType z=0; z<nrPoint.z(); z++)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType x=0; x<nrPoint.x(); x++)
                pointData.push_back(grd.get(Vec3i(x,y,z)));

    VTKWriter<TI> writer("3D Grid",path,true,grd.getBB(),grd.getNrCell(),grd.isCenter());
    writer.appendDatas("value",pointData.begin(),pointData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DVectorGridVTK(const string& path,const Grid<IndexType,TI>& vel,const TI& len)
{
    ASSERT(vel.getDim() == 2)

    std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;

    {
        const Vec3i nrPoints=vel.getNrPoint();
        for(sizeType xx=0; xx<nrPoints.x(); xx++)
            for(sizeType yy=0; yy<nrPoints.y(); yy++) {
                IndexType pt=vel.getPt(Vec3i(xx,yy,0));
                IndexType velDir=vel.sampleSafe2D(pt);
                if(velDir.norm() < EPS)
                    continue;
                if(len > 0.0f) {
                    velDir.normalize();
                    velDir*=len;
                }

                indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
                vertices.push_back(pt);
                vertices.push_back(pt+velDir);
            }
    }

    VTKWriter<TI> writer("Vec",boost::filesystem::path(path),true);
    writer.appendPoints(vertices.begin(),vertices.end());
    writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DMACGridVTK(const string& path,const MACGridType& mac)
{
    ASSERT(mac.getDim() == 2)

    const IndexType cellSz=mac.getCellSize();
    const IndexType offX(cellSz.x(),0.0f,0.0f);
    const IndexType offY(0.0f,cellSz.y(),0.0f);

    std::vector<IndexType,Eigen::aligned_allocator<Vec3> > vertices;
    std::vector<TI> colorData;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;

    {
        const GridType& uu=mac.getGu();
        const Vec3i nrPoints=mac.getGu().getNrPoint();
        for(sizeType xx=0; xx<nrPoints.x(); xx++)
            for(sizeType yy=0; yy<nrPoints.y(); yy++) {
                IndexType pt=uu.getPt(Vec3i(xx,yy,0))-offY*0.5f;
                indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
                vertices.push_back(pt);
                vertices.push_back(pt+offY);
                colorData.push_back(uu.get(Vec3i(xx,yy,0)));
                colorData.push_back(uu.get(Vec3i(xx,yy,0)));
            }
    }

    {
        const GridType& vv=mac.getGv();
        const Vec3i nrPoints=mac.getGv().getNrPoint();
        for(sizeType xx=0; xx<nrPoints.x(); xx++)
            for(sizeType yy=0; yy<nrPoints.y(); yy++) {
                IndexType pt=vv.getPt(Vec3i(xx,yy,0))-offX*0.5f;
                indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
                vertices.push_back(pt);
                vertices.push_back(pt+offX);
                colorData.push_back(vv.get(Vec3i(xx,yy,0)));
                colorData.push_back(vv.get(Vec3i(xx,yy,0)));
            }
    }

    VTKWriter<TI> writer("Vel",boost::filesystem::path(path),true);
    writer.appendPoints(vertices.begin(),vertices.end());
    writer.appendCells(indices.begin(),indices.end(),VTKWriter<typename IndexType::Scalar>::LINE);
    writer.appendCustomPointData("gridData",colorData.begin(),colorData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DVelocityGridVTK(const string& path,const MACGridType& vel,const TI& len)
{
    ASSERT(vel.getDim() == 2)

    std::vector<IndexType,Eigen::aligned_allocator<Vec3> > vertices;
    std::vector<TI> colorData;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indicesPt;

    {
        const GridType& uu=vel.getGu();
        const Vec3i nrPoints=vel.getGu().getNrPoint();
        for(sizeType xx=0; xx<nrPoints.x(); xx++)
            for(sizeType yy=0; yy<nrPoints.y(); yy++) {
                IndexType pt=uu.getPt(Vec3i(xx,yy,0));
                IndexType velDir=vel.sampleSafe2D(pt);
                if(velDir.norm() < EPS)
                    continue;
                if(len > 0.0f) {
                    velDir.normalize();
                    velDir*=len;
                }

                indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
                indicesPt.push_back(Vec2i(vertices.size(),0));
                vertices.push_back(pt);
                vertices.push_back(pt+velDir);
                colorData.push_back(0.0f);
                colorData.push_back(1.0f);
            }
    }

    {
        const GridType& vv=vel.getGv();
        const Vec3i nrPoints=vel.getGv().getNrPoint();
        for(sizeType xx=0; xx<nrPoints.x(); xx++)
            for(sizeType yy=0; yy<nrPoints.y(); yy++) {
                IndexType pt=vv.getPt(Vec3i(xx,yy,0));
                IndexType velDir=vel.sampleSafe2D(pt);
                if(velDir.norm() < EPS)
                    continue;
                if(len > 0.0f) {
                    velDir.normalize();
                    velDir*=len;
                }

                indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
                indicesPt.push_back(Vec2i(vertices.size(),0));
                vertices.push_back(pt);
                vertices.push_back(pt+velDir);
                colorData.push_back(0.0f);
                colorData.push_back(1.0f);
            }
    }

    VTKWriter<TI> writer("Vel",boost::filesystem::path(path),true);
    writer.appendPoints(vertices.begin(),vertices.end());
    writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
    writer.appendCells(indicesPt.begin(),indicesPt.end(),VTKWriter<TI>::POINT);
    writer.appendCustomPointData("Color",colorData.begin(),colorData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write3DVelocityGridVTK(const string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter)
{
    std::vector<IndexType,Eigen::aligned_allocator<Vec3> > vertices;
    std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;
    std::vector<TI> colorData;

    for(T x=vel.getBB()._minC.x(); x<vel.getBB()._maxC.x(); x+=sample.x())
        for(T y=vel.getBB()._minC.y(); y<vel.getBB()._maxC.y(); y+=sample.y())
            for(T z=vel.getBB()._minC.z(); z<vel.getBB()._maxC.z(); z+=sample.z()) {
                IndexType pt(x,y,z);
                if(jitter) {
                    pt.x()+=(rand()/(TI)RAND_MAX-0.5f)*sample.x()*0.5f;
                    pt.y()+=(rand()/(TI)RAND_MAX-0.5f)*sample.y()*0.5f;
                    pt.z()+=(rand()/(TI)RAND_MAX-0.5f)*sample.z()*0.5f;
                }
                advectRK2(pt,time,vel,vertices,indices,colorData);
            }

    VTKWriter<TI> writer("Vel",boost::filesystem::path(path),true);
    writer.appendPoints(vertices.begin(),vertices.end());
    writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
    writer.appendCustomPointData("Color",colorData.begin(),colorData.end());
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::advectRK2(const IndexType& from,const TI& time,const MACGridType& vel,
                                std::vector<IndexType,Eigen::aligned_allocator<IndexType> >& vertices,
                                std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> >& indices,
                                std::vector<TI>& colorData)
{
#define SUB_D 20
    IndexType curr=from;
    const TI stepLen=time/(TI)SUB_D;
    const TI stepCol=1.0f/(TI)SUB_D;
    
    vertices.push_back(curr);
    colorData.push_back(stepCol);

    for(TI i=0; i<time; i+=stepLen) {
        ValueType v=vel.sampleSafe(curr);
        IndexType mid=curr+IndexType(v.x(),v.y(),v.z())*stepLen*0.5f;

        v=vel.sampleSafe(mid);
        curr+=IndexType(v.x(),v.y(),v.z())*stepLen;

        indices.push_back(Vec2i(vertices.size()-1,vertices.size()));
        vertices.push_back(curr);
        colorData.push_back(stepCol);
    }
#undef SUB_D
}

//miscellaneous
template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::copyVelFromFunc(MACGridType& vel,const VelFunc<T2>& func)
{
    if(vel.getDim() >= 1)
        copyVelFromFunc(vel.getGu(),func,0);
    if(vel.getDim() >= 2)
        copyVelFromFunc(vel.getGv(),func,1);
    if(vel.getDim() >= 3)
        copyVelFromFunc(vel.getGw(),func,2);
}

template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::copyVelFromFunc(ScalarField& vel,const VelFunc<T2>& func,const sizeType& a)
{
    const Vec3i nrPoint=vel.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++)
                vel.get(Vec3i(x,y,z))=func(vel.getPt(Vec3i(x,y,z)))[a];
}

template <typename T,typename TI,typename TV>
template <typename T2>
void GridOp<T,TI,TV>::copyFromImplictFunc(GridType& to,const ImplicitFunc<T2>& func)
{
    typedef typename ScalarUtil<T2>::ScalarVec3 Vec3;
    const Vec3i nrPoint=to.getNrPoint();
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++) {
        INFOV("%lu",x)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++) {
                IndexType pt=to.getPt(Vec3i(x,y,z));
                to.get(Vec3i(x,y,z))=(T)func(Vec3((T2)pt.x(),(T2)pt.y(),(T2)pt.z()));
            }
    }
}

template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGrid(GridType& to,const Grid<T2,TI2>& from)
{
    const Vec3i nrPoint=to.getNrPoint();
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++) {
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++)
                to.get(Vec3i(x,y,z))=(T)from.sampleSafe(to.getPt(Vec3i(x,y,z)));
    }
}

template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGridOfSameGeometry(GridType& to,const Grid<T2,TI2>& from)
{
    const Vec3i nrPoint=to.getNrPoint();
    ASSERT(nrPoint==from.getNrPoint());
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++) {
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++)
                to.get(Vec3i(x,y,z))=(T)from.get(Vec3i(x,y,z));
    }
}

template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGridOfSameGeometry(MACGridType& to,const MACGrid<T2,TI2>& from)
{
    if(to.getDim() >= 1)
        copyFromOtherGridOfSameGeometry(to.getGu(),from.getGu());
    if(to.getDim() >= 2)
        copyFromOtherGridOfSameGeometry(to.getGv(),from.getGv());
    if(to.getDim() >= 3)
        copyFromOtherGridOfSameGeometry(to.getGw(),from.getGw());
}

template <typename T,typename TI,typename TV>
template <typename COMPARE,typename ALLOC>
void GridOp<T,TI,TV>::getAllDifferentValues(const GridType& g,std::set<T,COMPARE,ALLOC>& vals)
{
    const Vec3i nrPoint=g.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++)
                vals.insert(g.get(Vec3i(x,y,z)));
}

template <typename T,typename TI,typename TV>
template <typename COMPARE,typename ALLOC>
void GridOp<T,TI,TV>::getAllDifferentValues(const MACGridType& g,std::set<T,COMPARE,ALLOC>& vals)
{
    if(g.getDim() >= 1)
        getAllDifferentValues(g.getGu(),vals);
    if(g.getDim() >= 2)
        getAllDifferentValues(g.getGv(),vals);
    if(g.getDim() >= 3)
        getAllDifferentValues(g.getGw(),vals);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter(const MACGridType& from,Grid<ValueType,T>& to)
{
    ASSERT(from.getNrCell() == to.getNrCell())
    ASSERT(to.isCenter())
    if(to.getDim() == 2)
        fromFaceToCenter2D(from,to);
    else
        fromFaceToCenter3D(from,to);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter3D(const MACGridType& from,Grid<ValueType,T>& to)
{
    const Vec3i nrCell=to.getNrCell();
    for(sizeType x=0; x<nrCell.x(); x++)
        for(sizeType y=0; y<nrCell.y(); y++)
            for(sizeType z=0; z<nrCell.z(); z++) {
                ValueType& val=to.get(Vec3i(x,y,z));
                val(0)=(from.getGu().get(Vec3i(x,y,z))+from.getGu().get(Vec3i(x+1,y,z)))*0.5f;
                val(1)=(from.getGv().get(Vec3i(x,y,z))+from.getGv().get(Vec3i(x,y+1,z)))*0.5f;
                val(2)=(from.getGw().get(Vec3i(x,y,z))+from.getGw().get(Vec3i(x,y,z+1)))*0.5f;
            }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter2D(const MACGridType& from,Grid<ValueType,T>& to)
{
    const Vec3i nrCell=to.getNrCell();
    for(sizeType x=0; x<nrCell.x(); x++)
        for(sizeType y=0; y<nrCell.y(); y++) {
            ValueType& val=to.get(Vec3i(x,y,0));
            val(0)=(from.getGu().get(Vec3i(x,y,0))+from.getGu().get(Vec3i(x+1,y,0)))*0.5f;
            val(1)=(from.getGv().get(Vec3i(x,y,0))+from.getGv().get(Vec3i(x,y+1,0)))*0.5f;
            val(2)=0.0f;
        }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace(const Grid<ValueType,T>& from,MACGridType& to)
{
    ASSERT(from.getNrCell() == to.getNrCell())
    ASSERT(from.isCenter())
    if(to.getDim() == 2)
        fromCenterToFace2D(from,to);
    else
        fromCenterToFace3D(from,to);
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace3D(const Grid<ValueType,T>& from,MACGridType& to)
{
    const Vec3i nrCell=to.getNrCell();
    for(sizeType x=0; x<nrCell.x(); x++)
        for(sizeType y=0; y<nrCell.y(); y++)
            for(sizeType z=0; z<nrCell.z(); z++) {
                //X
                if(x == 0)
                    to.getGu().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).x();
                else
                    to.getGu().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).x()+from.get(Vec3i(x-1,y,z)).x())*0.5f;
                if(x == nrCell.x()-1)
                    to.getGu().get(Vec3i(x+1,y,z))=from.get(Vec3i(x,y,z)).x();
                else
                    to.getGu().get(Vec3i(x+1,y,z))=(from.get(Vec3i(x,y,z)).x()+from.get(Vec3i(x+1,y,z)).x())*0.5f;

                //Y
                if(y == 0)
                    to.getGv().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).y();
                else
                    to.getGv().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).y()+from.get(Vec3i(x,y-1,z)).y())*0.5f;
                if(y == nrCell.y()-1)
                    to.getGv().get(Vec3i(x,y+1,z))=from.get(Vec3i(x,y,z)).y();
                else
                    to.getGv().get(Vec3i(x,y+1,z))=(from.get(Vec3i(x,y,z)).y()+from.get(Vec3i(x,y+1,z)).y())*0.5f;

                //Z
                if(z == 0)
                    to.getGw().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).z();
                else
                    to.getGw().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).z()+from.get(Vec3i(x,y-1,z)).z())*0.5f;
                if(z == nrCell.z()-1)
                    to.getGw().get(Vec3i(x,y+1,z))=from.get(Vec3i(x,y,z)).z();
                else
                    to.getGw().get(Vec3i(x,y+1,z))=(from.get(Vec3i(x,y,z)).z()+from.get(Vec3i(x,y+1,z)).z())*0.5f;
            }
}

template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace2D(const Grid<ValueType,T>& from,MACGridType& to)
{
    const Vec3i nrCell=to.getNrCell();
    for(sizeType x=0; x<nrCell.x(); x++)
        for(sizeType y=0; y<nrCell.y(); y++) {
            //X
            if(x == 0)
                to.getGu().get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0)).x();
            else
                to.getGu().get(Vec3i(x,y,0))=(from.get(Vec3i(x,y,0)).x()+from.get(Vec3i(x-1,y,0)).x())*0.5f;
            if(x == nrCell.x()-1)
                to.getGu().get(Vec3i(x+1,y,0))=from.get(Vec3i(x,y,0)).x();
            else
                to.getGu().get(Vec3i(x+1,y,0))=(from.get(Vec3i(x,y,0)).x()+from.get(Vec3i(x+1,y,0)).x())*0.5f;

            //Y
            if(y == 0)
                to.getGv().get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0)).y();
            else
                to.getGv().get(Vec3i(x,y,0))=(from.get(Vec3i(x,y,0)).y()+from.get(Vec3i(x,y-1,0)).y())*0.5f;
            if(y == nrCell.y()-1)
                to.getGv().get(Vec3i(x,y+1,0))=from.get(Vec3i(x,y,0)).y();
            else
                to.getGv().get(Vec3i(x,y+1,0))=(from.get(Vec3i(x,y,0)).y()+from.get(Vec3i(x,y+1,0)).y())*0.5f;
        }
}

PRJ_END

#endif
