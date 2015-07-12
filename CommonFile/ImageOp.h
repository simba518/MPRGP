#ifndef IMAGE_OP_H
#define IMAGE_OP_H

#include "GridOp.h"

PRJ_BEGIN

template <typename WEIGHT>
struct DisjointSetElem {
    sizeType _rank;	//rank based merging
    sizeType _p;	//pointer to parent
    sizeType _size;	//size of current subtree, size=0 is a tag to avoid joining
    WEIGHT _Int;	//another criterion for merging
};

template <typename WEIGHT,typename ELEM=DisjointSetElem<WEIGHT> >
class DisjointSet
{
public:
    DisjointSet():_nrSet(0) {}
    DisjointSet(sizeType elements) {
        resize(elements);
    }
    void resize(sizeType elements) {
        _elts.clear();
        _elts.resize(elements);
        _nrSet=elements;
        for(sizeType i=0; i < elements; i++) {
            _elts[i]._rank=0;
            _elts[i]._size=1;
            _elts[i]._p=i;
            _elts[i]._Int=(WEIGHT)0.0f;
        }
    }
    virtual ~DisjointSet() {}
    sizeType find(sizeType x) {
        sizeType y=x;
        while(y != _elts[y]._p)
            y=_elts[y]._p;
        _elts[x]._p=y;
        return y;
    }
    void join(sizeType x,sizeType y) {
        if(x == y)
            return;
        if(_elts[x]._rank > _elts[y]._rank) {
            _elts[y]._p=x;
            _elts[x]._size+=_elts[y]._size;
        } else {
            _elts[x]._p=y;
            _elts[y]._size+=_elts[x]._size;
            if(_elts[x]._rank == _elts[y]._rank)
                _elts[y]._rank++;
        }
        _nrSet--;
    }
    void joinSafe(sizeType x,sizeType y) {
        if(_elts[x]._size == 0 || _elts[y]._size == 0)
            return;
        join(find(x),find(y));
    }
    sizeType size(sizeType x) const {
        return _elts[x]._size;
    }
    sizeType numSets() const {
        return _nrSet;
    }
    std::vector<ELEM> _elts;
    sizeType _nrSet;
};

template <typename T,typename TI,typename WEIGHT>
class GraphBasedSegmentation
{
public:
    class Edge
    {
    public:
        sizeType _A;
        sizeType _B;
        WEIGHT _weight;
        Edge():_A(-1),_B(-1),_weight(0.0f) {}
        Edge(const sizeType& A,const sizeType& B,const WEIGHT& weight):_A(A),_B(B),_weight(weight) {}
        bool operator<(const Edge& b) const {
            return _weight < b._weight;
        }
    };
    GraphBasedSegmentation(Grid<sizeType,TI>& output,const Grid<T,TI>& input,const WEIGHT& k)
        :_tag(input.getNrPoint().prod()),_k(k) {
        const Vec3i nrPoint=input.getNrPoint();
        if(input.getDim() == 2) {
            sizeType nrEdge=(nrPoint.x()-1)*nrPoint.y()+
                            nrPoint.x()*(nrPoint.y()-1);
            _weights.resize(nrEdge);

            sizeType index=0;
            for(sizeType x=0; x<nrPoint.x()-1; x++)
                for(sizeType y=0; y<nrPoint.y(); y++) {
                    sizeType A=input.getIndexNoAlign(Vec3i(x,y,0));
                    sizeType B=input.getIndexNoAlign(Vec3i(x+1,y,0));
                    _weights[index++]=edge(A,B,diff(input.get(Vec3i(x,y,0)),input.get(Vec3i(x+1,y,0))));
                }
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y()-1; y++) {
                    sizeType A=input.getIndexNoAlign(Vec3i(x,y,0));
                    sizeType B=input.getIndexNoAlign(Vec3i(x,y+1,0));
                    _weights[index++]=edge(A,B,diff(input.get(Vec3i(x,y,0)),input.get(Vec3i(x,y+1,0))));
                }
            ASSERT(index == nrEdge)
        } else {
            sizeType nrEdge=(nrPoint.x()-1)*nrPoint.y()*nrPoint.z()+
                            nrPoint.x()*(nrPoint.y()-1)*nrPoint.z()+
                            nrPoint.x()*nrPoint.y()*(nrPoint.z()-1);
            _weights.resize(nrEdge);

            sizeType index=0;
            for(sizeType x=0; x<nrPoint.x()-1; x++)
                for(sizeType y=0; y<nrPoint.y(); y++)
                    for(sizeType z=0; z<nrPoint.z(); z++) {
                        sizeType A=input.getIndexNoAlign(Vec3i(x,y,z));
                        sizeType B=input.getIndexNoAlign(Vec3i(x+1,y,z));
                        _weights[index++]=edge(A,B,diff(input.get(Vec3i(x,y,z)),input.get(Vec3i(x+1,y,z))));
                    }
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y()-1; y++)
                    for(sizeType z=0; z<nrPoint.z(); z++) {
                        sizeType A=input.getIndexNoAlign(Vec3i(x,y,z));
                        sizeType B=input.getIndexNoAlign(Vec3i(x,y+1,z));
                        _weights[index++]=edge(A,B,diff(input.get(Vec3i(x,y,z)),input.get(Vec3i(x,y+1,z))));
                    }
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y(); y++)
                    for(sizeType z=0; z<nrPoint.z()-1; z++) {
                        sizeType A=input.getIndexNoAlign(Vec3i(x,y,z));
                        sizeType B=input.getIndexNoAlign(Vec3i(x,y,z+1));
                        _weights[index++]=edge(A,B,diff(input.get(Vec3i(x,y,z)),input.get(Vec3i(x,y,z+1))));
                    }
            ASSERT(index == nrEdge)
        }

        segment();

        output.makeSameGeometry(input);
        assemble(output);
    }
    virtual ~GraphBasedSegmentation() {}
    void segment() {
        adjustWeight();
        typedef typename std::vector<Edge>::const_iterator iter;
        std::sort(_weights.begin(),_weights.end());
        for(iter beg=_weights.begin(),end=_weights.end(); beg!=end; beg++) {
            const Edge& edg=*beg;
            const sizeType baseA=_tag.find(edg._A);
            const sizeType baseB=_tag.find(edg._B);
            if(baseA != baseB && edg._weight <= MInt(baseA,baseB)) {
                _tag.join(baseA,baseB);
                _tag._elts[_tag.find(baseA)]._Int=edg._weight;
            }
        }
    }
    void assemble(Grid<sizeType,TI>& output) {
        typedef typename std::vector<typename DisjointSet<WEIGHT>::Elem>::iterator Iter;

        const sizeType nrSet=_tag.numSets();
        sizeType id=0;
        for(Iter beg=_tag._elts.begin(),end=_tag._elts.end(); beg!=end; beg++)
            if(beg->_p == (sizeType)(beg-_tag._elts.begin()))
                beg->_rank=id++;
        ASSERT(id == nrSet)

        const Vec3i nrPoint=output.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const sizeType id=output.getIndexNoAlign(Vec3i(x,y,z));
                    output.get(Vec3i(x,y,z))=_tag._elts[_tag.find(id)]._rank;
                }
    }
    WEIGHT MInt(const sizeType& baseA,const sizeType& baseB) const {
        typedef typename DisjointSet<WEIGHT>::Elem ElemType;
        const ElemType& elemA=_tag._elts[baseA];
        const ElemType& elemB=_tag._elts[baseB];
        return std::min<WEIGHT>(elemA._Int+_k/elemA._size,elemB._Int+_k/elemB._size);
    }
    virtual WEIGHT diff(const T& a,const T& b) {
        return std::abs(a-b);
    }
    virtual void adjustWeight() {}
protected:
    DisjointSet<WEIGHT> _tag;
    std::vector<Edge> _weights;
    WEIGHT _k;
};

#define TEST_JOINT(id)								\
N=id;												\
if(tag.isSafeIndex(N) && tag.get(N) != -1)			\
{													\
	const sizeType nid=tag.getIndexNoAlign(N);		\
	if(nid != cid)									\
		set.joinSafe(cid,nid);						\
}
template<typename T,typename TI>
class ImageOp
{
    typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
    typedef typename ScalarUtil<T>::ScalarMat3 Mat3Type;
    typedef typename ScalarUtil<T>::ScalarMat4 Mat4Type;
public:
    //extract main feature
    static void extractMainFeatureLevelSet(Grid<T,TI>& levelSet,const T& featureSize,const Grid<T,TI>* solid) {
        const Vec3Type cellSz=levelSet.getCellSize();
        const Vec3i nrExpand=ceil(Vec3Type(std::abs(featureSize)*2.0f/std::max(cellSz.x(),EPS),
                                           std::abs(featureSize)*2.0f/std::max(cellSz.y(),EPS),
                                           std::abs(featureSize)*2.0f/std::max(cellSz.z(),EPS)));

        Grid<T,TI> levelSetExpand=levelSet;
        levelSetExpand.expand(nrExpand,std::abs(featureSize));
        GridOp<T,TI>::copyFromOtherGrid(levelSetExpand,levelSet);
        GridOp<T,TI>::reinitialize(levelSetExpand);

        if(solid)
            constraint(levelSetExpand,*solid);

        levelSetExpand.add(-featureSize);
        GridOp<T,TI>::reinitialize(levelSetExpand);
        levelSetExpand.add(featureSize);
        GridOp<T,TI>::copyFromOtherGrid(levelSet,levelSetExpand);
    }
    static void constraint(Grid<T,TI>& phi,const Grid<T,TI>& solid) {
        const Vec3i nrPoint=phi.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const Vec3i id(x,y,z);
                    scalar& phiVal=phi.get(id);
                    phiVal=std::max(phiVal,-solid.sampleSafe(phi.getPt(id)));
                }
    }
    //ridge based segmentation for level set
    static void segmentLevelSet(const Grid<T,TI>& levelSet,Grid<sizeType,TI>& tag,const T& k,const T& radius) {
        const Vec3i nrPoint=levelSet.getNrPoint();
        const sizeType nrErode=(sizeType)std::ceil(radius/levelSet.getCellSize().maxCoeff());

        //initial segment
        {
            Grid<T,TI> levelSetGrad;
            levelSetGrad.makeSameGeometry(levelSet);
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y(); y++)
                    for(sizeType z=0; z<nrPoint.z(); z++) {
                        const T gradMag=levelSet.sampleSafeGrad(levelSet.getPt(Vec3i(x,y,z))).norm();
                        levelSetGrad.get(Vec3i(x,y,z))=gradMag-1.0f;
                    }
            GraphBasedSegmentation<T,TI,T> segment(tag,levelSetGrad,k);
        }

        //identify important group
        std::set<sizeType> mainRegion;
        Grid<sizeType,TI> tagMain=tag;
        {
            for(sizeType i=0; i<nrErode; i++)
                erode(tagMain);
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y(); y++)
                    for(sizeType z=0; z<nrPoint.z(); z++) {
                        const sizeType& id=tagMain.get(Vec3i(x,y,z));
                        if(id != -1)
                            mainRegion.insert(id);
                    }
        }

        //identify weakly connected region and do more accurate segmentation
        Grid<unsigned char,TI> mainTag;
        {
            mainTag.makeSameGeometry(tag);
            mainTag.init(false);
            for(sizeType x=0; x<nrPoint.x(); x++)
                for(sizeType y=0; y<nrPoint.y(); y++)
                    for(sizeType z=0; z<nrPoint.z(); z++)
                        if(mainRegion.find(tag.get(Vec3i(x,y,z))) != mainRegion.end())
                            mainTag.get(Vec3i(x,y,z))=true;

            floodFill(tagMain);
            while(guardedDilate(tagMain,mainTag));
            tag.swap(tagMain);
        }

        //extend the main region to fill the whole domain
        while(dilate(tag));
    }
    //erode, return can still erode
    static bool erode(Grid<sizeType,TI>& tag) {
        bool hasFeature=false;
        Grid<sizeType,TI> tmp=tag;
        const Vec3i nrPoint=tag.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++)
                    hasFeature|=erodeCell(tmp,tag,Vec3i(x,y,z));
        tag.swap(tmp);
        return hasFeature;
    }
    static bool erodeCell(Grid<sizeType,TI>& tmp,const Grid<sizeType,TI>& tag,const Vec3i& id) {
        sizeType& val=tmp.get(id);
        if(val == -1)
            return false;

        for(sizeType xx=id.x()-1; xx<=id.x()+1; xx++)
            for(sizeType yy=id.y()-1; yy<=id.y()+1; yy++)
                for(sizeType zz=id.z()-1; zz<=id.z()+1; zz++)
                    if(tag.getSafe(Vec3i(xx,yy,zz)) != val) {
                        val=-1;
                        return false;
                    }

        return true;
    }
    //dilate, return can still dilate
    static bool dilate(Grid<sizeType,TI>& tag) {
        bool change=false;
        Grid<sizeType,TI> tmp=tag;
        const Vec3i nrPoint=tag.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++)
                    dilateCell(tmp,tag,Vec3i(x,y,z),change);
        tag.swap(tmp);
        return change;
    }
    static void dilateCell(Grid<sizeType,TI>& tmp,const Grid<sizeType,TI>& tag,const Vec3i& id,bool& change) {
        sizeType& val=tmp.get(id);
        if(val != -1)
            return;

        for(sizeType xx=id.x()-1; xx<=id.x()+1; xx++)
            for(sizeType yy=id.y()-1; yy<=id.y()+1; yy++)
                for(sizeType zz=id.z()-1; zz<=id.z()+1; zz++) {
                    const sizeType neighVal=tag.getSafe(Vec3i(xx,yy,zz));
                    if(neighVal != -1) {
                        val=neighVal;
                        change=true;
                        return;
                    }
                }
    }
    //flood fill
    static void floodFill(Grid<sizeType,TI>& tag,bool corner=true) {
        DisjointSet<sizeType> set(tag.getNrPoint().prod());

        //joint neighbor cells
        const Vec3i nrPoint=tag.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    if(tag.get(Vec3i(x,y,z)) == -1)
                        set._elts[tag.getIndexNoAlign(Vec3i(x,y,z))]._Int=-1;
                    else {
                        if(corner)
                            testJoint8(set,tag,Vec3i(x,y,z));
                        else
                            testJoint4(set,tag,Vec3i(x,y,z));
                    }
                }

        //assemble
        typedef typename std::vector<DisjointSetElem<sizeType> >::iterator Iter;
        //const sizeType nrSet=set.numSets();
        sizeType id=0;
        for(Iter beg=set._elts.begin(),end=set._elts.end(); beg!=end; beg++)
            if(beg->_Int == 0 && beg->_p == (sizeType)(beg-set._elts.begin()))
                beg->_rank=id++;

        //assign tags
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const sizeType id=tag.getIndexNoAlign(Vec3i(x,y,z));
                    if(tag.get(Vec3i(x,y,z)) != -1)
                        tag.get(Vec3i(x,y,z))=set._elts[set.find(id)]._rank;
                }
    }
    static void testJoint4(DisjointSet<sizeType> &set,const Grid<sizeType,TI>& tag,const Vec3i& id) {
        if(tag.get(id) == -1)
            return;

        const sizeType cid=tag.getIndexNoAlign(id);
        Vec3i N;
        //X
        TEST_JOINT(id-Vec3i(1,0,0))
        TEST_JOINT(id+Vec3i(1,0,0))
        //Y
        TEST_JOINT(id-Vec3i(0,1,0))
        TEST_JOINT(id+Vec3i(0,1,0))
        //Z
        if(tag.getDim() == 3) {
            TEST_JOINT(id-Vec3i(0,0,1))
            TEST_JOINT(id+Vec3i(0,0,1))
        }

    }
    static void testJoint8(DisjointSet<sizeType> &set,const Grid<sizeType,TI>& tag,const Vec3i& id) {
        if(tag.get(id) == -1)
            return;

        Vec3i N;
        const sizeType cid=tag.getIndexNoAlign(id);
        for(sizeType xx=id.x()-1; xx<=id.x()+1; xx++)
            for(sizeType yy=id.y()-1; yy<=id.y()+1; yy++)
                for(sizeType zz=id.z()-1; zz<=id.z()+1; zz++) {
                    TEST_JOINT(Vec3i(xx,yy,zz))
                }
    }
    //guarded retag
    static bool guardedDilate(Grid<sizeType,TI>& tag,const Grid<unsigned char,TI>& guard) {
        bool change=false;
        Grid<sizeType,TI> tmp=tag;
        const Vec3i nrPoint=tag.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++)
                    guardedDilateCell(tmp,tag,guard,Vec3i(x,y,z),change);
        tag.swap(tmp);
        return change;
    }
    static void guardedDilateCell(Grid<sizeType,TI>& tmp,const Grid<sizeType,TI>& tag,const Grid<unsigned char,TI>& guard,const Vec3i& id,bool& change) {
        sizeType& val=tmp.get(id);
        if(!guard.get(id) || val != -1)
            return;

        for(sizeType xx=id.x()-1; xx<=id.x()+1; xx++)
            for(sizeType yy=id.y()-1; yy<=id.y()+1; yy++)
                for(sizeType zz=id.z()-1; zz<=id.z()+1; zz++) {
                    const sizeType id=tag.getSafe(Vec3i(xx,yy,zz));
                    if(id != -1) {
                        val=id;
                        change=true;
                        return;
                    }
                }
    }
};
#undef TEST_JOINT

PRJ_END

#endif
