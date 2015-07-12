#ifndef GRID_BASIC_H
#define GRID_BASIC_H

#include "Config.h"
#include "MathBasic.h"
#include "IO.h"

PRJ_BEGIN

template <typename T,typename TI,typename TG=vector<T,Eigen::aligned_allocator<T> > >
struct Grid : public HasMagic {
    friend class DataStructureCL;
public:
    typedef T value_type;
    typedef typename Eigen::Matrix<T,3,1> ValueType;
    typedef typename Eigen::Matrix<TI,3,1> IndexType;
    typedef typename Eigen::Matrix<T,3,3> MatrixType;
    Grid():HasMagic(0xFCFCCFCFFCFCCFCF),_szPoint(0,0,0) {}
    virtual ~Grid() {}
    bool write(ostream &os) const {
        if(!HasMagic::writeMagic(os))
            return false;

        if(!writeBinaryData(_off,os).good())
            return false;
        if(!writeBinaryData(_szCell,os).good())
            return false;
        if(!writeBinaryData(_invSzCell,os).good())
            return false;
        if(!writeBinaryData(_minClamp,os).good())
            return false;
        if(!writeBinaryData(_maxClamp,os).good())
            return false;
        if(!writeBinaryData(_maxInterp,os).good())
            return false;
        if(!writeBinaryData(_szPoint,os).good())
            return false;
        if(!writeBinaryData(_maxIndex,os).good())
            return false;
        if(!writeBinaryData(_bb,os).good())
            return false;
        if(!writeBinaryData(_move,os).good())
            return false;
        if(!writeBinaryData(_szPointAligned,os).good())
            return false;
        if(!writeBinaryData(_stride,os).good())
            return false;
        if(!writeBinaryData(_align,os).good())
            return false;
        if(!writeVector(_grid,os).good())
            return false;
        if(!writeBinaryData(_dim,os).good())
            return false;

        return true;
    }
    bool read(istream &is) {
        if(!HasMagic::readMagic(is))
            return false;

        if(!readBinaryData(_off,is).good())
            return false;
        if(!readBinaryData(_szCell,is).good())
            return false;
        if(!readBinaryData(_invSzCell,is).good())
            return false;
        if(!readBinaryData(_minClamp,is).good())
            return false;
        if(!readBinaryData(_maxClamp,is).good())
            return false;
        if(!readBinaryData(_maxInterp,is).good())
            return false;
        if(!readBinaryData(_szPoint,is).good())
            return false;
        if(!readBinaryData(_maxIndex,is).good())
            return false;
        if(!readBinaryData(_bb,is).good())
            return false;
        if(!readBinaryData(_move,is).good())
            return false;
        if(!readBinaryData(_szPointAligned,is).good())
            return false;
        if(!readBinaryData(_stride,is).good())
            return false;
        if(!readBinaryData(_align,is).good())
            return false;
        if(!readVector(_grid,is).good())
            return false;
        if(!readBinaryData(_dim,is).good())
            return false;

        return true;
    }
    bool writeASCII1D(ostream& os) const {
        os << getNrPoint().x() << std::endl;
        for(sizeType x=0; x<getNrPoint().x(); x++)
            os << (T)get(Vec3i(x,0,0)) << " "  << std::endl;
        return os.good();
    }
    bool writeASCII2D(ostream& os) const {
        //y
        //|
        //|
        //|
        //0--------x
        os << getNrPoint().x() << " " << getNrPoint().y() << std::endl;
        for(sizeType y=getNrPoint().y()-1; y>=0; y--) {
            for(sizeType x=0; x<getNrPoint().x(); x++)
                os << (T)get(Vec3i(x,y,0)) << " ";
            os << std::endl;
        }
        return os.good();
    }
    void reset(const Vec3i& nrCell,const BBox<TI>& bb,const T& val,bool center=true,sizeType align=4,bool shadow=false,bool ZFirst=true,const IndexType& move=IndexType::Zero()) {
        const IndexType extent=bb.getExtent();
        _dim=3;
        if(extent.z() == 0.0f)_dim--;
        if(extent.y() == 0.0f)_dim--;

        _bb=bb;
        _move=move;
        _szPoint=center ? nrCell : (nrCell+Vec3i::Constant(1));
        //fill alignment
        _szPointAligned=_szPoint;
        if(_dim == 3)
            _szPointAligned[2]=((_szPoint[2]+align-1)/align)*align;
        else if(_dim == 2)
            _szPointAligned[1]=((_szPoint[1]+align-1)/align)*align;
        else
            _szPointAligned[0]=((_szPoint[0]+align-1)/align)*align;
        //set all pending dimension to 1
        for(sizeType j=_dim; j<3; j++)
            _szPoint(j)=_szPointAligned(j)=1;
        if((_ZFirst=ZFirst))
            _stride=Vec3i(_szPointAligned.y()*_szPointAligned.z(),_szPointAligned.z(),1);
        else
            _stride=Vec3i(1,_szPointAligned.x(),_szPointAligned.x()*_szPointAligned.y());
        _maxIndex=_szPoint-Vec3i::Constant(1);
        _align=align;

        _off=center ? 0.5f : 0.0f;
        _szCell=IndexType(extent.x()/nrCell.x(),
                          (_dim < 2) ? 0.0f : extent.y()/nrCell.y(),
                          (_dim < 3) ? 0.0f : extent.z()/nrCell.z());
        _invSzCell=IndexType(1.0f/_szCell.x(),
                             (_dim < 2) ? ScalarUtil<TI>::scalar_max : 1.0f/_szCell.y(),
                             (_dim < 3) ? ScalarUtil<TI>::scalar_max : 1.0f/_szCell.z());

        _minClamp=IndexType::Zero();
        _maxClamp=IndexType((TI)_szPoint.x()-1.0f,(TI)_szPoint.y()-1.0f,(TI)_szPoint.z()-1.0f);
        _maxInterp=_maxIndex;
        if(_dim >= 3)
            _maxInterp.z()--;
        if(_dim >= 2)
            _maxInterp.y()--;
        if(_dim >= 1)
            _maxInterp.x()--;

        if(!shadow)
            init(val);
        else {
            TG tmp;
            _grid.swap(tmp);
        }
    }
    template <typename T2,typename TI2>
    void makeSameGeometry(const Grid<T2,TI2>& other,bool shadow=false,bool ZFirst=true,sizeType align=-1) {
        BBox<TI> bb;
        bb.copy(other.getRawBB());
        if(align == -1)
            align=other.getAlign();
        reset(other.getNrCell(),bb,T(),other.isCenter(),align,shadow,ZFirst,other.getMove());
    }
    void expand(const Vec3i& nr,const T& val) {
        IndexType nrD(_szCell.x()*(scalar)nr.x(),
                      _szCell.y()*(scalar)nr.y(),
                      _szCell.y()*(scalar)nr.z());

        //expand cell size
        Vec3i nrCell=getNrCell()+nr*2;

        //expand bounding box size
        BBox<scalar> bb=getRawBB();
        bb._minC-=nrD;
        bb._maxC+=nrD;
        for(sizeType j=2; j>=_dim; j--)
            bb._maxC[j]=bb._minC[j];

        reset(nrCell,bb,val,isCenter(),getAlign(),false,_ZFirst,getMove());
    }
    void decimate(const Grid<T,TI>& child,const sizeType& decimate,const T& val,bool uniformScale,sizeType align=4) {
        _align=align;

        //we only support node grid decimate
        //ASSERT(child._off == 0.0f);
        ASSERT(decimate > 1);

        if(uniformScale) {
            //cell size
            IndexType cellSz=child.getCellSize()*(TI)decimate;
            Vec3i nrCell=ceil((IndexType)((child._bb.getExtent().array()/cellSz.array()).matrix()));
            nrCell=compMax(nrCell,Vec3i(2,2,2));	//ensure correct interpolation
            ASSERT(nrCell.x()*nrCell.y()*nrCell.z() > 1)

            BBox<TI> bb;
            bb._minC=child._bb._minC;
            bb._maxC=(cellSz.array()*IndexType((TI)nrCell.x(),(TI)nrCell.y(),(TI)nrCell.z()).array()).matrix()+bb._minC;
            for(sizeType j=child.getDim(); j<3; j++)
                bb._maxC(j)=bb._minC(j);
            reset(nrCell,bb,val,child.isCenter(),align,child.getMove());
        } else {
            Vec3i nrCell=child.getNrCell()/decimate;
            nrCell=compMax(nrCell,Vec3i(2,2,2));
            ASSERT(nrCell.x()*nrCell.y()*nrCell.z() > 1)

            reset(nrCell,child._bb,val,child.isCenter(),align,child.getMove());
        }
    }
public:
    const sizeType getAlign() const {
        return _align;
    }
    const sizeType getIndexNoAlign(const Vec3i& index) const {
        return Vec3i(_szPoint.y()*_szPoint.z(),_szPoint.z(),1).dot(index);
    }
    const sizeType getIndex(const Vec3i& index) const {
        return _stride.dot(index);
    }
    const Vec3i getIndex(const sizeType& index) const {
        return Vec3i(index/_stride.x(),(index%_stride.x())/_stride.y(),index%_stride.y());
    }
    FORCE_INLINE const T& operator[](const Vec3i& index) const {
        return get(index);
    }
    FORCE_INLINE const T& operator[](const sizeType& index) const {
        return _grid[index];
    }
    const T& get(const Vec3i& index) const {
        return _grid[getIndex(index)];
    }
    const T& getSafe(const Vec3i& index) const {
        return _grid[getIndex(compMax(compMin(index,_maxIndex),Vec3i(0,0,0)))];
    }
    bool isSafeIndex(const Vec3i& id) const {
        return compGE(id,Vec3i(0,0,0)) && compLE(id,_maxIndex);
    }
    const Vec3i& getNrPoint() const {
        return _szPoint;
    }
    const Vec3i& getMaxInterp() const {
        return _maxInterp;
    }
    Vec3i getNrCell() const {
        if(_dim == 1)
            return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i(1,0,0);
        else if(_dim == 2)
            return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i(1,1,0);
        else
            return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i::Constant(1);
    }
    BBox<TI> getBB() const {
        return BBox<TI>(_bb._minC+_move,_bb._maxC+_move);
    }
    BBox<TI> getRawBB() const{
        return _bb;
    }
    const IndexType& getMove() const{return _move;}
    FORCE_INLINE T& operator[](const Vec3i& index) {
        return get(index);
    }
    FORCE_INLINE T& operator[](const sizeType& index) {
        return _grid[index];
    }
    T& get(const Vec3i& index) {
        return _grid[getIndex(index)];
    }
    T& getSafe(const Vec3i& index) {
        return _grid[getIndex(compMax(compMin(index,_maxIndex),Vec3i(0,0,0)))];
    }
    const IndexType& getInvCellSize() const {
        return _invSzCell;
    }
    const IndexType& getCellSize() const {
        return _szCell;
    }
    IndexType getIndexFrac(const IndexType& pos) const {
        IndexType ret=((pos-_move-_bb._minC).array()*getInvCellSize().array()).matrix()-IndexType::Constant(_off);
        if(_dim < 2)
            ret.y()=0.0f;
        if(_dim < 3)
            ret.z()=0.0f;
        return ret;
    }
    IndexType getIndexFracSafe(const IndexType& pos) const {
        return compMax(_minClamp,compMin(getIndexFrac(pos),_maxClamp));
    }
    IndexType getPt(const Vec3i& cell) const {
        return _move+_bb._minC+(getCellSize().array()*IndexType(TI(cell.x())+_off,TI(cell.y())+_off,TI(cell.z())+_off).array()).matrix();
    }
    IndexType getCellCtr(const Vec3i& cell) const {
        return _move+_bb._minC+(getCellSize().array()*IndexType(TI(cell.x())+0.5f,TI(cell.y())+0.5f,TI(cell.z())+0.5f).array()).matrix();
    }
    IndexType getCellCorner(const Vec3i& cell) const {
        return _move+_bb._minC+(getCellSize().array()*IndexType(TI(cell.x())     ,TI(cell.y())     ,TI(cell.z())).array()).matrix();
    }
    const sizeType getSzLinear() const {
        return _szPointAligned.prod();
    }
    const sizeType getSzLinearNoAlign() const {
        return _szPoint.x()*_szPoint.y()*_szPoint.z();
    }
#define DETERMINE_VALUE3D(FRAC_FUNC)    \
ASSERT(_dim == 3)   \
ASSERT(_szPoint.x() > 1)    \
ASSERT(_szPoint.y() > 1)    \
ASSERT(_szPoint.z() > 1);    \
IndexType frac=FRAC_FUNC(pos);    \
const Vec3i base=compMin(floor(frac),_maxInterp);    \
frac-=IndexType((TI)base.x(),(TI)base.y(),(TI)base.z());    \
const sizeType id000=getIndex(base);    \
const sizeType id100=id000+_stride.x();    \
const sizeType id010=id000+_stride.y();    \
const sizeType id110=id100+_stride.y();    \
const sizeType id001=id000+_stride.z();    \
const sizeType id101=id001+_stride.x();    \
const sizeType id011=id001+_stride.y();    \
const sizeType id111=id101+_stride.y();

#define DETERMINE_VALUE2D(FRAC_FUNC)    \
ASSERT(_dim == 2)    \
ASSERT(_szPoint.x() > 1)    \
ASSERT(_szPoint.y() > 1);    \
IndexType frac=FRAC_FUNC(pos);    \
const Vec3i base=compMin(floor(frac),_maxInterp);    \
frac-=IndexType((TI)base.x(),(TI)base.y(),0.0f);    \
const sizeType id000=getIndex(base);    \
const sizeType id100=id000+_stride.x();    \
const sizeType id010=id000+_stride.y();    \
const sizeType id110=id100+_stride.y();

#define DETERMINE_VALUE1D(FRAC_FUNC)    \
ASSERT(_dim == 1)    \
ASSERT(_szPoint.x() > 1);    \
IndexType frac=FRAC_FUNC(pos);    \
const Vec3i base=compMin(floor(frac),_maxInterp);    \
frac-=IndexType((TI)base.x(),0.0f,0.0f);    \
const sizeType id000=getIndex(base);    \
const sizeType id100=id000+_stride.x();
    //sample value
    T sample(const IndexType& pos) const {
        return _dim == 1 ? sample1D(pos) : _dim == 2 ? sample2D(pos) : sample3D(pos);
    }
    T sampleSafe(const IndexType& pos) const {
        return _dim == 1 ? sampleSafe1D(pos) : _dim == 2 ? sampleSafe2D(pos) : sampleSafe3D(pos);
    }
    T sample3D(const IndexType& pos) const {
        DETERMINE_VALUE3D(getIndexFrac)
        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z());//higher z plane
    }
    T sampleSafe3D(const IndexType& pos) const {
        DETERMINE_VALUE3D(getIndexFracSafe)
        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z());//higher z plane
    }
    T sample2D(const IndexType& pos) const {
        DETERMINE_VALUE2D(getIndexFrac)
        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y());	//higher z plane
    }
    T sampleSafe2D(const IndexType& pos) const {
        DETERMINE_VALUE2D(getIndexFracSafe)
        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y());
    }
    T sample1D(const IndexType& pos) const {
        DETERMINE_VALUE1D(getIndexFrac)
        return interp1D(_grid[id000],_grid[id100],frac.x());	//higher z plane
    }
    T sampleSafe1D(const IndexType& pos) const {
        DETERMINE_VALUE1D(getIndexFracSafe)
        return interp1D(_grid[id000],_grid[id100],frac.x());
    }
    //sample stencil
    T getSampleStencil(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        return _dim == 1 ? getSampleStencil1D(pos,coefs,pts) : 
               _dim == 2 ? getSampleStencil2D(pos,coefs,pts) : 
                           getSampleStencil3D(pos,coefs,pts);
    }
    T getSampleStencilSafe(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        return _dim == 1 ? getSampleStencilSafe1D(pos,coefs,pts) : 
               _dim == 2 ? getSampleStencilSafe2D(pos,coefs,pts) : 
                           getSampleStencilSafe3D(pos,coefs,pts);
    }
    T getSampleStencil3D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE3D(getIndexFrac)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        pts[2]=base+Vec3i(0,1,0);
        pts[3]=base+Vec3i(1,1,0);
        pts[4]=base+Vec3i(0,0,1);
        pts[5]=base+Vec3i(1,0,1);
        pts[6]=base+Vec3i(0,1,1);
        pts[7]=base+Vec3i(1,1,1);
        return stencil3D(coefs,frac.x(),frac.y(),frac.z());//higher z plane
    }
    T getSampleStencilSafe3D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE3D(getIndexFracSafe)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        pts[2]=base+Vec3i(0,1,0);
        pts[3]=base+Vec3i(1,1,0);
        pts[4]=base+Vec3i(0,0,1);
        pts[5]=base+Vec3i(1,0,1);
        pts[6]=base+Vec3i(0,1,1);
        pts[7]=base+Vec3i(1,1,1);
        return stencil3D(coefs,frac.x(),frac.y(),frac.z());//higher z plane
    }
    T getSampleStencil2D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE2D(getIndexFrac)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        pts[2]=base+Vec3i(0,1,0);
        pts[3]=base+Vec3i(1,1,0);
        return stencil2D(coefs,frac.x(),frac.y());//higher z plane
    }
    T getSampleStencilSafe2D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE2D(getIndexFracSafe)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        pts[2]=base+Vec3i(0,1,0);
        pts[3]=base+Vec3i(1,1,0);
        return stencil2D(coefs,frac.x(),frac.y());//higher z plane
    }
    T getSampleStencil1D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE1D(getIndexFrac)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        return stencil1D(coefs,frac.x());//higher z plane
    }
    T getSampleStencilSafe1D(const IndexType& pos,TI* coefs,Vec3i* pts) const {
        DETERMINE_VALUE1D(getIndexFracSafe)
        pts[0]=base+Vec3i(0,0,0);
        pts[1]=base+Vec3i(1,0,0);
        return stencil1D(coefs,frac.x());//higher z plane
    }
    //sample value with minmax
    T sample3D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE3D(getIndexFrac)
        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z(),minV,maxV);//higher z plane
    }
    T sampleSafe3D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE3D(getIndexFracSafe)
        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z(),minV,maxV);//higher z plane
    }
    T sample2D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE2D(getIndexFrac)
        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y(),minV,maxV);	//higher z plane
    }
    T sampleSafe2D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE2D(getIndexFracSafe)
        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y(),minV,maxV);
    }
    T sample1D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE1D(getIndexFrac)
        return interp1D(_grid[id000],_grid[id100],frac.x(),minV,maxV);	//higher z plane
    }
    T sampleSafe1D(const IndexType& pos,T& minV,T& maxV) const {
        DETERMINE_VALUE1D(getIndexFracSafe)
        return interp1D(_grid[id000],_grid[id100],frac.x(),minV,maxV);
    }
    //sample value with default
    T sample3D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE3D(getIndexFrac)
        if(valid[id000] == 0 || valid[id100] == 0 ||
           valid[id010] == 0 || valid[id110] == 0 ||
           valid[id001] == 0 || valid[id101] == 0 ||
           valid[id011] == 0 || valid[id111] == 0)
            return def;

        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z());//higher z plane
    }
    T sampleSafe3D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE3D(getIndexFracSafe)
        if(valid[id000] == 0 || valid[id100] == 0 ||
           valid[id010] == 0 || valid[id110] == 0 ||
           valid[id001] == 0 || valid[id101] == 0 ||
           valid[id011] == 0 || valid[id111] == 0)
            return def;

        return interp3D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],	 //lower z plane
                        _grid[id001],_grid[id101],
                        _grid[id011],_grid[id111],
                        frac.x(),frac.y(),frac.z());//higher z plane
    }
    T sample2D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE2D(getIndexFrac)
        if(valid[id000] == 0 || valid[id100] == 0 ||
           valid[id010] == 0 || valid[id110] == 0)
            return def;

        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y());	//higher z plane
    }
    T sampleSafe2D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE2D(getIndexFracSafe)
        if(valid[id000] == 0 || valid[id100] == 0 ||
           valid[id010] == 0 || valid[id110] == 0)
            return def;

        return interp2D(_grid[id000],_grid[id100],
                        _grid[id010],_grid[id110],
                        frac.x(),frac.y());
    }
    T sample1D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE1D(getIndexFrac)
        if(valid[id000] == 0 || valid[id100] == 0)
            return def;

        return interp1D(_grid[id000],_grid[id100],frac.x());	//higher z plane
    }
    T sampleSafe1D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def) const {
        DETERMINE_VALUE1D(getIndexFracSafe)
        if(valid[id000] == 0 || valid[id100] == 0)
            return def;

        return interp1D(_grid[id000],_grid[id100],frac.x());
    }
    //sample gradient
    ValueType sampleGrad(const IndexType& pos) const {
        return _dim == 1 ? sample1DGrad(pos) : _dim == 2 ? sample2DGrad(pos) : sample3DGrad(pos);
    }
    ValueType sampleSafeGrad(const IndexType& pos) const {
        return _dim == 1 ? sampleSafe1DGrad(pos) : _dim == 2 ? sampleSafe2DGrad(pos) : sampleSafe3DGrad(pos);
    }
    ValueType sample3DGrad(const IndexType& pos) const {
        DETERMINE_VALUE3D(getIndexFrac)
        return ValueType(interp2D(_grid[id100]-_grid[id000],
                                  _grid[id110]-_grid[id010],
                                  _grid[id101]-_grid[id001],
                                  _grid[id111]-_grid[id011],
                                  frac.y(),frac.z())*_invSzCell.x(),

                         interp2D(_grid[id010]-_grid[id000],
                                  _grid[id110]-_grid[id100],
                                  _grid[id011]-_grid[id001],
                                  _grid[id111]-_grid[id101],
                                  frac.x(),frac.z())*_invSzCell.y(),

                         interp2D(_grid[id001]-_grid[id000],
                                  _grid[id101]-_grid[id100],
                                  _grid[id011]-_grid[id010],
                                  _grid[id111]-_grid[id110],
                                  frac.x(),frac.y())*_invSzCell.z() );//higher z plane
    }
    ValueType sampleSafe3DGrad(const IndexType& pos) const {
        DETERMINE_VALUE3D(getIndexFracSafe)
		ValueType ret	(interp2D(_grid[id100]-_grid[id000],
                                  _grid[id110]-_grid[id010],
                                  _grid[id101]-_grid[id001],
                                  _grid[id111]-_grid[id011],
                                  frac.y(),frac.z())*_invSzCell.x(),

                         interp2D(_grid[id010]-_grid[id000],
                                  _grid[id110]-_grid[id100],
                                  _grid[id011]-_grid[id001],
                                  _grid[id111]-_grid[id101],
                                  frac.x(),frac.z())*_invSzCell.y(),

                         interp2D(_grid[id001]-_grid[id000],
                                  _grid[id101]-_grid[id100],
                                  _grid[id011]-_grid[id010],
                                  _grid[id111]-_grid[id110],
                                  frac.x(),frac.y())*_invSzCell.z() );
        if(pos.x() <= _bb._minC.x() || pos.x() >= _bb._maxC.x())
            ret(0)=0.0f;
        if(pos.y() <= _bb._minC.y() || pos.y() >= _bb._maxC.y())
            ret(1)=0.0f;
        if(pos.z() <= _bb._minC.z() || pos.z() >= _bb._maxC.z())
            ret(2)=0.0f;
        return ret;
    }
    ValueType sample2DGrad(const IndexType& pos) const {
        DETERMINE_VALUE2D(getIndexFrac)
        return ValueType(interp1D(_grid[id100]-_grid[id000],_grid[id110]-_grid[id010],frac.y())*_invSzCell.x(),
                         interp1D(_grid[id010]-_grid[id000],_grid[id110]-_grid[id100],frac.x())*_invSzCell.y(),
                         0.0f);
    }
    ValueType sampleSafe2DGrad(const IndexType& pos) const {
		DETERMINE_VALUE2D(getIndexFracSafe)
        ValueType ret	(interp1D(_grid[id100]-_grid[id000],_grid[id110]-_grid[id010],frac.y())*_invSzCell.x(),
                         interp1D(_grid[id010]-_grid[id000],_grid[id110]-_grid[id100],frac.x())*_invSzCell.y(),
                         0.0f);
        if(pos.x() <= _bb._minC.x() || pos.x() >= _bb._maxC.x())
            ret(0)=0.0f;
        if(pos.y() <= _bb._minC.y() || pos.y() >= _bb._maxC.y())
            ret(1)=0.0f;
        return ret;
    }
    ValueType sample1DGrad(const IndexType& pos) const {
        DETERMINE_VALUE1D(getIndexFrac)
        return ValueType((_grid[id100]-_grid[id000])*_invSzCell.x(),
                         0.0f,
                         0.0f);	//higher z plane
    }
    ValueType sampleSafe1DGrad(const IndexType& pos) const {
        DETERMINE_VALUE1D(getIndexFracSafe)
		ValueType ret	((_grid[id100]-_grid[id000])*_invSzCell.x(),
                         0.0f,
                         0.0f);
        if(pos.x() <= _bb._minC.x() || pos.x() >= _bb._maxC.x())
            ret(0)=0.0f;
        return ret;
    }
    //sample laplace
    void sampleLaplace3D(const IndexType& pos,MatrixType& lap) const {
        DETERMINE_VALUE3D(getIndexFrac)

        lap.row(0)=ValueType(0.0f,
                             interp1D((_grid[id110]-_grid[id010])-(_grid[id100]-_grid[id000]),
                                      (_grid[id111]-_grid[id011])-(_grid[id101]-_grid[id001]),frac.z()),
                             interp1D((_grid[id101]-_grid[id001])-(_grid[id100]-_grid[id000]),
                                      (_grid[id111]-_grid[id011])-(_grid[id110]-_grid[id010]),frac.y()))*
                   (_invSzCell.y()*_invSzCell.z());

        lap.row(1)=ValueType(interp1D((_grid[id110]-_grid[id100])-(_grid[id010]-_grid[id000]),
                                      (_grid[id111]-_grid[id101])-(_grid[id011]-_grid[id001]),frac.z()),
                             0.0f,
                             interp1D((_grid[id011]-_grid[id001])-(_grid[id010]-_grid[id000]),
                                      (_grid[id111]-_grid[id101])-(_grid[id110]-_grid[id100]),frac.x()))*
                   (_invSzCell.x()*_invSzCell.z());

        lap.row(2)=ValueType(interp1D((_grid[id101]-_grid[id100])-(_grid[id001]-_grid[id000]),
                                      (_grid[id111]-_grid[id110])-(_grid[id011]-_grid[id010]),frac.y()),
                             interp1D((_grid[id011]-_grid[id010])-(_grid[id001]-_grid[id000]),
                                      (_grid[id111]-_grid[id110])-(_grid[id101]-_grid[id100]),frac.x()),
                             0.0f)*
                   (_invSzCell.x()*_invSzCell.y());
    }
    void sampleLaplaceSafe3D(const IndexType& pos,MatrixType& lap) const {
        DETERMINE_VALUE3D(getIndexFracSafe)

        lap.row(0)=ValueType(0.0f,
                             interp1D((_grid[id110]-_grid[id010])-(_grid[id100]-_grid[id000]),
                                      (_grid[id111]-_grid[id011])-(_grid[id101]-_grid[id001]),frac.z()),
                             interp1D((_grid[id101]-_grid[id001])-(_grid[id100]-_grid[id000]),
                                      (_grid[id111]-_grid[id011])-(_grid[id110]-_grid[id010]),frac.y()))*
                   (_invSzCell.y()*_invSzCell.z());

        lap.row(1)=ValueType(interp1D((_grid[id110]-_grid[id100])-(_grid[id010]-_grid[id000]),
                                      (_grid[id111]-_grid[id101])-(_grid[id011]-_grid[id001]),frac.z()),
                             0.0f,
                             interp1D((_grid[id011]-_grid[id001])-(_grid[id010]-_grid[id000]),
                                      (_grid[id111]-_grid[id101])-(_grid[id110]-_grid[id100]),frac.x()))*
                   (_invSzCell.x()*_invSzCell.z());

        lap.row(2)=ValueType(interp1D((_grid[id101]-_grid[id100])-(_grid[id001]-_grid[id000]),
                                      (_grid[id111]-_grid[id110])-(_grid[id011]-_grid[id010]),frac.y()),
                             interp1D((_grid[id011]-_grid[id010])-(_grid[id001]-_grid[id000]),
                                      (_grid[id111]-_grid[id110])-(_grid[id101]-_grid[id100]),frac.x()),
                             0.0f)*
                   (_invSzCell.x()*_invSzCell.y());
    }
    //simple operation
    T dot(const Grid& other) const {
        const sizeType n=(sizeType)_grid.size();
        T ret=0.0f;
        for(sizeType i=0; i<n; i++) {
            ret+=_grid[i]*other._grid[i];
        }
        return ret;
    }
    T getAbsMax() const {
        const sizeType n=(sizeType)_grid.size();
        T maxV=-ScalarUtil<T>::scalar_max;
        for(sizeType i=0; i<n; i++) {
            const T& val=_grid[i];
            if(std::abs(val) > maxV)
                maxV=std::abs(val);
        }
        return maxV;
    }
    T getAbsMaxVec3() const {
        const sizeType n=(sizeType)_grid.size();
        T maxV=T::Constant(-ScalarUtil<typename T::Scalar>::scalar_max);
        for(sizeType i=0; i<n; i++) {
            const T& val=_grid[i];
            if(std::abs(val[0]) > maxV[0])
                maxV[0]=std::abs(val[0]);
            if(std::abs(val[1]) > maxV[1])
                maxV[1]=std::abs(val[1]);
            if(std::abs(val[2]) > maxV[2])
                maxV[2]=std::abs(val[2]);
        }
        return maxV;
    }
    T sum() const {
        T ret=0.0f;
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            ret+=_grid[i];
        return ret;
    }
    void minMax(T& minV,T& maxV) const {
        const sizeType n=(sizeType)_grid.size();
        minV=ScalarUtil<T>::scalar_max;
        maxV=-minV;
        for(sizeType i=0; i<n; i++) {
            const T& val=_grid[i];
            minV=mmin<T>(val,minV);
            maxV=mmax<T>(val,maxV);
        }
    }
    void add(const Grid& other) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]+=other._grid[i];
    }
    void add(const T& coef) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]+=coef;
    }
    void sub(const Grid& other) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]-=other._grid[i];
    }
	void min(const Grid& other) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
			_grid[i]=std::min(_grid[i],other._grid[i]);
    }
    void sub(const T& coef) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]-=coef;
    }
    void addScaled(const Grid& other,const T& coef) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]+=other._grid[i]*coef;
    }
    void mul(const T& coef) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++)
            _grid[i]*=coef;
    }
    void clamp(const T& maxVal) {
        const sizeType n=(sizeType)_grid.size();
        for(sizeType i=0; i<n; i++) {
            T val=std::abs(_grid[i]);
            if(val > maxVal)
                _grid[i]*=maxVal/val;
        }
    }
    //dimension reduction
    void getSlice(Grid<T,TI>& slice,const sizeType& dim0,const sizeType& dim1,const sizeType& dim2,const T& x) const {
        Vec3i nrCell(getNrCell()[dim0],getNrCell()[dim1],0);
        BBox<TI> bb;
        bb._minC.x()=_bb._minC[dim0];
        bb._maxC.x()=_bb._maxC[dim0];
        bb._minC.y()=_bb._minC[dim1];
        bb._maxC.y()=_bb._maxC[dim1];
        bb._minC.z()=0.0f;
        bb._maxC.z()=0.0f;
        slice.reset(nrCell,bb,get(Vec3i(0,0,0)),isCenter(),getAlign(),false,_ZFirst,getMove());

        sizeType sliceId=(sizeType)std::floor((x-_bb._minC[dim2])/getCellSize()[dim2]);

        Vec3i nrPoint=slice.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++) {
                Vec3i id;
                id[dim0]=x;
                id[dim1]=y;
                id[dim2]=sliceId;
                slice.get(Vec3i(x,y,0))=get(id);
            }
    }
    void getSliceYZ(Grid<T,TI>& slice,const T& x) const {
        getSlice(slice,1,2,0,x);
    }
    void getSliceXZ(Grid<T,TI>& slice,const T& y) const {
        getSlice(slice,0,2,1,y);
    }
    void getSliceXY(Grid<T,TI>& slice,const T& z) const {
        getSlice(slice,0,1,2,z);
    }
    Grid<T,TI> getSliceYZ(const T& x) const {
        Grid<T,TI> ret;
        getSlice(ret,1,2,0,x);
        return ret;
    }
    Grid<T,TI> getSliceXZ(const T& y) const {
        Grid<T,TI> ret;
        getSlice(ret,0,2,1,y);
        return ret;
    }
    Grid<T,TI> getSliceXY(const T& z) const {
        Grid<T,TI> ret;
        getSlice(ret,0,1,2,z);
        return ret;
    }
    //move the grid
    IndexType moveBy(Vec3i nrGrid,bool copyData=true)
    {
        Vec3i minC=Vec3i::Zero();
        Vec3i maxC=_szPoint;
        Vec3i dir=Vec3i::Ones();
        
        IndexType rel=IndexType::Zero();
        if(_dim>=1){
            if(nrGrid[0] < 0){
                std::swap(minC[0],maxC[0]);
                minC[0]--;maxC[0]--;
                dir[0]=-1;
            }
        }else{
            maxC[0]=1;
            nrGrid[0]=0;
        }rel[0]=(TI)nrGrid[0]*_szCell[0];

        if(_dim>=2){
            if(nrGrid[1] < 0){
                std::swap(minC[1],maxC[1]);
                minC[1]--;maxC[1]--;
                dir[1]=-1;
            }
        }else{
            maxC[1]=1;
            nrGrid[1]=0;
        }rel[1]=(TI)nrGrid[1]*_szCell[1];

        if(_dim>=3){
            if(nrGrid[2] < 0){
                std::swap(minC[2],maxC[2]);
                minC[2]--;maxC[2]--;
                dir[2]=-1;
            }
        }else{
            maxC[2]=1;
            nrGrid[2]=0;
        }rel[2]=(TI)nrGrid[2]*_szCell[2];

        move(rel);
        if(copyData)
        for(sizeType x=minC[0];x!=maxC[0];x+=dir[0])
        for(sizeType y=minC[1];y!=maxC[1];y+=dir[1])
        for(sizeType z=minC[2];z!=maxC[2];z+=dir[2])
        {
            Vec3i id(x,y,z);
            get(id)=getSafe(id+nrGrid);
        }
        return rel;
    }
    void resetMove(){_move.setZero();}
    void move(const IndexType& delta) {
        _move+=delta;
    }
    //dangerous methods
    const TG& data() const {
        return _grid;
    }
    void setData(const Grid<T,TI>& other){
        for(sizeType x=0;x<_szPoint.x();x++)
        for(sizeType y=0;y<_szPoint.y();y++)
        for(sizeType z=0;z<_szPoint.z();z++)
            get(Vec3i(x,y,z))=other.get(Vec3i(x,y,z));
    }
    const T& get(const sizeType& index) const {
        return _grid[index];
    }
    T& get(const sizeType& index) {
        return _grid[index];
    }
    const Vec3i& getStride() const {
        return _stride;
    }
    bool isCenter() const {
        return _off == 0.5f;
    }
    void init(const T& val) {
        {
            TG tmp;
            _grid.swap(tmp);
        }
        _grid.assign(_szPointAligned.x()*_szPointAligned.y()*_szPointAligned.z(),val);
    }
    const sizeType& getDim() const {
        return _dim;
    }
    bool getZFirst() const {
        return _ZFirst;
    }
    void swap(Grid& other) {
        ASSERT(_grid.size() == other._grid.size())
        std::swap(_off,other._off);
        std::swap(_szCell,other._szCell);
        std::swap(_invSzCell,other._invSzCell);
        std::swap(_minClamp,other._minClamp);
        std::swap(_maxClamp,other._maxClamp);
        std::swap(_maxInterp,other._maxInterp);
        std::swap(_szPoint,other._szPoint);
        std::swap(_maxIndex,other._maxIndex);
        std::swap(_bb,other._bb);
        std::swap(_szPointAligned,other._szPointAligned);
        std::swap(_stride,other._stride);
        std::swap(_align,other._align);
        _grid.swap(other._grid);
        std::swap(_dim,other._dim);
    }
    void setBB(const BBox<TI>& bb) {
        ASSERT((bb.getExtent()-_bb.getExtent()).norm() < ScalarUtil<TI>::scalar_eps);
        _bb._minC=bb._minC;
        _bb._maxC=bb._maxC;
        _move=IndexType::Zero();
    }
protected:
    //param
    TI _off;
    IndexType _szCell;
    IndexType _invSzCell;
    IndexType _minClamp;
    IndexType _maxClamp;
    Vec3i _maxInterp;
    Vec3i _szPoint;
    Vec3i _maxIndex;
    BBox<TI> _bb;
    IndexType _move;
    //memory
    Vec3i _szPointAligned;
    Vec3i _stride;
    sizeType _align;
    TG _grid;
    sizeType _dim;
    bool _ZFirst;
};

template <typename T,typename TI,typename TG=vector<T,Eigen::aligned_allocator<T> > >
struct MACGrid : public HasMagic {
    friend class dataStructureCL;
public:
    typedef typename Eigen::Matrix<T,3,1> ValueType;
    typedef typename Eigen::Matrix<TI,3,1> IndexType;
    MACGrid():HasMagic(0xAAAABBBBAAAABBBB) {}
    ~MACGrid() {}
    bool write(ostream &os) const {
        if(!HasMagic::writeMagic(os))
            return false;

        if(!writeBinaryData(_bb,os).good())
            return false;
        if(!writeBinaryData(_bbSafe,os).good())
            return false;
        if(!writeBinaryData(_move,os).good())
            return false;
        if(!writeBinaryData(_cellSz,os).good())
            return false;
        if(!writeBinaryData(_nrCell,os).good())
            return false;
        if(!writeBinaryData(_maxIndex,os).good())
            return false;
        if(!writeBinaryData(_dim,os).good())
            return false;

        bool ret=true;
        if(_dim >= 1 && ret)
            ret=ret && _gu.write(os);
        if(_dim >= 2 && ret)
            ret=ret && _gv.write(os);
        if(_dim >= 3 && ret)
            ret=ret && _gw.write(os);
        return ret;
    }
    bool read(istream &is) {
        if(!HasMagic::readMagic(is))
            return false;

        if(!readBinaryData(_bb,is).good())
            return false;
        if(!readBinaryData(_bbSafe,is).good())
            return false;
        if(!readBinaryData(_move,is).good())
            return false;
        if(!readBinaryData(_cellSz,is).good())
            return false;
        if(!readBinaryData(_nrCell,is).good())
            return false;
        if(!readBinaryData(_maxIndex,is).good())
            return false;
        if(!readBinaryData(_dim,is).good())
            return false;

        bool ret=true;
        if(_dim >= 1 && ret)
            ret=ret && _gu.read(is);
        if(_dim >= 2 && ret)
            ret=ret && _gv.read(is);
        if(_dim >= 3 && ret)
            ret=ret && _gw.read(is);
        return ret;
    }
    template <typename T2,typename TI2,typename TG2>
    void reset(const Grid<T2,TI2,TG2>& ref,bool shadow=false,bool edge=false) {
        //ASSERT(ref.isCenter())
        _cellSz=IndexType((TI)ref.getCellSize().x(),(TI)ref.getCellSize().y(),(TI)ref.getCellSize().z());
        _nrCell=ref.getNrCell();
        _maxIndex=compMax((Vec3i)(_nrCell-Vec3i::Ones()),Vec3i::Zero());

        _dim=ref.getDim();
        _bb.copy(ref.getBB());
        
        typename MACGrid<T2,TI2,TG2>::IndexType refM=ref.getMove();
        _move=IndexType((TI)refM[0],(TI)refM[1],(TI)refM[2]);

        if(ref.getDim() >= 1)
            resetGu(ref,shadow,edge);
        if(ref.getDim() >= 2)
            resetGv(ref,shadow,edge);
        if(ref.getDim() >= 3)
            resetGw(ref,shadow,edge);
        resetBB();
    }
    void resetBB()
    {
        _bbSafe=_bb;
        if(getGu().getDim() >= 1) {
            _bbSafe._minC.x()+=_cellSz.x()*0.75f;
            _bbSafe._maxC.x()-=_cellSz.x()*0.75f;
        }
        if(getGu().getDim() >= 2) {
            _bbSafe._minC.y()+=_cellSz.y()*0.75f;
            _bbSafe._maxC.y()-=_cellSz.y()*0.75f;
        }
        if(getGu().getDim() >= 3) {
            _bbSafe._minC.z()+=_cellSz.z()*0.75f;
            _bbSafe._maxC.z()-=_cellSz.z()*0.75f;
        }
    }
    template <typename T2,typename TI2,typename TG2>
    void makeSameGeometry(const MACGrid<T2,TI2,TG2>& other) {
        _bb.copy(other.getBBRaw());
        _bbSafe.copy(other.getBBSafe());

        typename MACGrid<T2,TI2,TG2>::IndexType otherM=other.getMove();
        _move=IndexType((TI)otherM[0],(TI)otherM[1],(TI)otherM[2]);

        typename MACGrid<T2,TI2,TG2>::IndexType cellSz=other.getCellSize();
        _cellSz=IndexType(cellSz.x(),cellSz.y(),cellSz.z());
        _nrCell=other.getNrCell();
        _maxIndex=other.getMaxIndex();
        _dim=other.getDim();

        if(_dim >= 1)
            getGu().makeSameGeometry(other.getGu());
        if(_dim >= 2)
            getGv().makeSameGeometry(other.getGv());
        if(_dim >= 3)
            getGw().makeSameGeometry(other.getGw());
    }
public:
    const sizeType getAlign() const {
        return _gu.getAlign();
    }
    const Vec3i& getNrCell() const {
        return _nrCell;
    }
    const Vec3i& getMaxIndex() const {
        return _maxIndex;
    }
    BBox<TI> getBBRaw() const {
        return _bb;
    }
    BBox<TI> getBB() const {
        return BBox<TI>(_bb._minC+_move,_bb._maxC+_move);
    }
    BBox<TI> getBBSafe() const{
        return _bbSafe;
    }
    const IndexType& getMove() const{return _move;}
    const IndexType& getInvCellSize() const {
        return _gu.getInvCellSize();
    }
    const IndexType& getCellSize() const {
        return _gu.getCellSize();
    }
    ValueType get(const Vec3i& index) const {
        if(getDim() == 1)
            return get1D(index);
        else if(getDim() == 2)
            return get2D(index);
        else
            return get3D(index);
    }
    ValueType getSafe(const Vec3i& index) const {
        return get(compMin(compMax(index,Vec3i::Zero()),_maxIndex));
    }
    ValueType get1D(const Vec3i& index) const {
        ASSERT(getDim() == 1)
        return ValueType((_gu.get(index)+_gu.get(index+Vec3i(1,0,0)))*0.5f,
                         0.0f,0.0f);
    }
    ValueType getSafe1D(const Vec3i& index) const {
        return get1D(compMin(compMax(index,Vec3i::Zero()),_maxIndex));
    }
    ValueType get2D(const Vec3i& index) const {
        ASSERT(getDim() == 2)
        return ValueType((_gu.get(index)+_gu.get(index+Vec3i(1,0,0)))*0.5f,
                         (_gv.get(index)+_gv.get(index+Vec3i(0,1,0)))*0.5f,
                         0.0f);
    }
    ValueType getSafe2D(const Vec3i& index) const {
        return get2D(compMin(compMax(index,Vec3i::Zero()),_maxIndex));
    }
    ValueType get3D(const Vec3i& index) const {
        ASSERT(getDim() == 3)
        return ValueType((_gu.get(index)+_gu.get(index+Vec3i(1,0,0)))*0.5f,
                         (_gv.get(index)+_gv.get(index+Vec3i(0,1,0)))*0.5f,
                         (_gw.get(index)+_gw.get(index+Vec3i(0,0,1)))*0.5f);
    }
    ValueType getSafe3D(const Vec3i& index) const {
        return get2D(compMin(compMax(index,Vec3i::Zero()),_maxIndex));
    }
    ValueType sample(const IndexType& pos) const {
        if(getDim() == 1)
            return sample1D(pos);
        else if(getDim() == 2)
            return sample2D(pos);
        else return sample3D(pos);
    }
    //safe version
    ValueType sampleSafe(const IndexType& pos) const {
        return sample(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move));
    }
    ValueType sample3D(const IndexType& pos) const {
        ASSERT(getDim() == 3);
        return ValueType(_gu.sample3D(pos),_gv.sample3D(pos),_gw.sample3D(pos));
    }
    ValueType sampleSafe3D(const IndexType& pos) const {
        return sample3D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move));
    }
    ValueType sample2D(const IndexType& pos) const {
        ASSERT(getDim() == 2);
        return ValueType(_gu.sample2D(pos),_gv.sample2D(pos),0.0f);
    }
    ValueType sampleSafe2D(const IndexType& pos) const {
        return sample2D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move));
    }
    ValueType sample1D(const IndexType& pos) const {
        ASSERT(getDim() == 1);
        return ValueType(_gu.sample1D(pos),0.0f,0.0f);
    }
    ValueType sampleSafe1D(const IndexType& pos) const {
        return sample1D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move));
    }
    //safe version with default value
    ValueType sample3D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        ASSERT(getDim() == 3);
        return ValueType(_gu.sample3D(pos,valid.getGu(),def.x()),_gv.sample3D(pos,valid.getGv(),def.y()),_gw.sample3D(pos,valid.getGw(),def.z()));
    }
    ValueType sampleSafe3D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        return sample3D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move),valid,def);
    }
    ValueType sample2D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        ASSERT(getDim() == 2);
        return ValueType(_gu.sample2D(pos,valid.getGu(),def.x()),_gv.sample2D(pos,valid.getGv(),def.y()),0.0f);
    }
    ValueType sampleSafe2D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        return sample2D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move),valid,def);
    }
    ValueType sample1D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        ASSERT(getDim() == 1);
        return ValueType(_gu.sample1D(pos,valid.getGu(),def.x()),0.0f,0.0f);
    }
    ValueType sampleSafe1D(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def) const {
        return sample1D(compMin(compMax(pos,_bbSafe._minC+_move),_bbSafe._maxC+_move),valid,def);
    }
    //simple operation
    T dot(const MACGrid& other) const {
        T ret=0.0f;
        if(getDim() >= 1)
            ret+=_gu.dot(other._gu);
        if(getDim() >= 2)
            ret+=_gv.dot(other._gv);
        if(getDim() >= 3)
            ret+=_gw.dot(other._gw);
        return ret;
    }
    ValueType getAbsMax() const {
        ValueType ret=ValueType::Zero();
        if(getDim() >= 1)
            ret.x()=_gu.getAbsMax();
        if(getDim() >= 2)
            ret.y()=_gv.getAbsMax();
        if(getDim() >= 3)
            ret.z()=_gw.getAbsMax();
        return ret;
    }
    void minMax(ValueType& minV,ValueType& maxV) const {
        minV=maxV=ValueType::Zero();
        if(getDim() >= 1)
            _gu.minMax(minV.x(),maxV.x());
        if(getDim() >= 2)
            _gv.minMax(minV.y(),maxV.y());
        if(getDim() >= 3)
            _gw.minMax(minV.z(),maxV.z());
    }
    void add(const MACGrid& other) {
        if(getDim() >= 1)
            _gu.add(other._gu);
        if(getDim() >= 2)
            _gv.add(other._gv);
        if(getDim() >= 3)
            _gw.add(other._gw);
    }
    void sub(const MACGrid& other) {
        if(getDim() >= 1)
            _gu.sub(other._gu);
        if(getDim() >= 2)
            _gv.sub(other._gv);
        if(getDim() >= 3)
            _gw.sub(other._gw);
    }
	void min(const MACGrid& other) {
        if(getDim() >= 1)
            _gu.min(other._gu);
        if(getDim() >= 2)
            _gv.min(other._gv);
        if(getDim() >= 3)
            _gw.min(other._gw);
    }
    void mul(const T& other) {
        if(getDim() >= 1)
            _gu.mul(other);
        if(getDim() >= 2)
            _gv.mul(other);
        if(getDim() >= 3)
            _gw.mul(other);
    }
    void addScaled(const MACGrid& other,const T& coef) {
        if(getDim() >= 1)
            _gu.addScaled(other._gu,coef);
        if(getDim() >= 2)
            _gv.addScaled(other._gv,coef);
        if(getDim() >= 3)
            _gw.addScaled(other._gw,coef);
    }
    void clamp(const T& maxVal) {
        if(getDim() >= 1)
            _gu.clamp(maxVal);
        if(getDim() >= 2)
            _gv.clamp(maxVal);
        if(getDim() >= 3)
            _gw.clamp(maxVal);
    }
     //move the grid
    void moveBy(const Vec3i& nrGrid,bool copyData=true)
    {
        IndexType rel;
        if(_dim>=1)rel=getGu().moveBy(nrGrid,copyData);
        if(_dim>=2)rel=getGv().moveBy(nrGrid,copyData);
        if(_dim>=3)rel=getGw().moveBy(nrGrid,copyData);
        _move=rel;
    }
    void resetMove()
    {
        if(_dim>=1)getGu().resetMove();
        if(_dim>=2)getGv().resetMove();
        if(_dim>=3)getGw().resetMove();
        _move.setZero();
    }
    //dangerous methods
    void setData(const MACGrid<T,TI>& other)
    {
        if(_dim>=1)getGu().setData(other.getGu());
        if(_dim>=2)getGv().setData(other.getGv());
        if(_dim>=2)getGw().setData(other.getGw());
    }
    void init(const ValueType& val) {
        if(getDim() >= 1)
            _gu.init(val.x());
        if(getDim() >= 2)
            _gv.init(val.y());
        if(getDim() >= 3)
            _gw.init(val.z());
    }
    const sizeType& getDim() const {
        return _dim;
    }
    Grid<T,TI,TG>& getGu() {
        return _gu;
    }
    Grid<T,TI,TG>& getGv() {
        return _gv;
    }
    Grid<T,TI,TG>& getGw() {
        return _gw;
    }
    Grid<T,TI,TG>& getComp(const sizeType& i) {
        return i==0 ? _gu : i==1 ? _gv : _gw;
    }
    const Grid<T,TI,TG>& getGu() const {
        return _gu;
    }
    const Grid<T,TI,TG>& getGv() const {
        return _gv;
    }
    const Grid<T,TI,TG>& getGw() const {
        return _gw;
    }
    const Grid<T,TI,TG>& getComp(const sizeType& i) const {
        return i==0 ? _gu : i==1 ? _gv : _gw;
    }
    void swap(MACGrid& other) {
        std::swap(_bb,other._bb);
        std::swap(_bbSafe,other._bbSafe);
        std::swap(_cellSz,other._cellSz);
        std::swap(_nrCell,other._nrCell);
        std::swap(_maxIndex,other._maxIndex);
        std::swap(_dim,other._dim);
        _gu.swap(other._gu);
        _gv.swap(other._gv);
        _gw.swap(other._gw);
    }
protected:
    template <typename T2,typename TI2,typename TG2>
    void resetGu(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge) {
        Vec3i nrCell=ref.getNrCell();
        BBox<TI> bb;
        bb.copy(ref.getBB());
        if(edge) {
            ASSERT(ref.getDim() > 1)
            if(ref.getDim() == 3) {
                bb._minC.x()+=_cellSz.x()*0.5f;
                bb._maxC.x()-=_cellSz.x()*0.5f;
                nrCell.x()--;
            }
        } else {
            if(ref.getDim() >= 2) {
                bb._minC.y()+=_cellSz.y()*0.5f;
                bb._maxC.y()-=_cellSz.y()*0.5f;
                nrCell.y()--;
            }
            if(ref.getDim() >= 3) {
                bb._minC.z()+=_cellSz.z()*0.5f;
                bb._maxC.z()-=_cellSz.z()*0.5f;
                nrCell.z()--;
            }
        }

        _gu.reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true,_move);
    }
    template <typename T2,typename TI2,typename TG2>
    void resetGv(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge) {
        Vec3i nrCell=ref.getNrCell();
        BBox<TI> bb;
        bb.copy(ref.getBB());
        if(edge) {
            ASSERT(ref.getDim() > 1)
            if(ref.getDim() == 3) {
                bb._minC.y()+=_cellSz.y()*0.5f;
                bb._maxC.y()-=_cellSz.y()*0.5f;
                nrCell.y()--;
            }
        } else {
            if(ref.getDim() >= 1) {
                bb._minC.x()+=_cellSz.x()*0.5f;
                bb._maxC.x()-=_cellSz.x()*0.5f;
                nrCell.x()--;
            }
            if(ref.getDim() >= 3) {
                bb._minC.z()+=_cellSz.z()*0.5f;
                bb._maxC.z()-=_cellSz.z()*0.5f;
                nrCell.z()--;
            }
        }

        _gv.reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true,_move);
    }
    template <typename T2,typename TI2,typename TG2>
    void resetGw(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge) {
        Vec3i nrCell=ref.getNrCell();
        BBox<TI> bb;
        bb.copy(ref.getRawBB());
        if(edge) {
            ASSERT(ref.getDim() > 1)
            bb._minC.z()+=_cellSz.z()*0.5f;
            bb._maxC.z()-=_cellSz.z()*0.5f;
            nrCell.z()--;
        } else {
            if(ref.getDim() >= 1) {
                bb._minC.x()+=_cellSz.x()*0.5f;
                bb._maxC.x()-=_cellSz.x()*0.5f;
                nrCell.x()--;
            }
            if(ref.getDim() >= 2) {
                bb._minC.y()+=_cellSz.y()*0.5f;
                bb._maxC.y()-=_cellSz.y()*0.5f;
                nrCell.y()--;
            }
        }

        _gw.reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true,_move);
    }
    //params
    BBox<TI> _bb;
    BBox<TI> _bbSafe;
    IndexType _move;
    IndexType _cellSz;
    Vec3i _nrCell;
    Vec3i _maxIndex;
    sizeType _dim;
    //data
    Grid<T,TI,TG> _gu;
    Grid<T,TI,TG> _gv;
    Grid<T,TI,TG> _gw;
};

typedef Grid<unsigned char,scalarF> TagFieldF;
typedef Grid<scalarF,scalarF> ScalarFieldF;
typedef Grid<Vec3f,scalarF> VectorFieldF;
typedef MACGrid<unsigned char,scalarF> MACTagFieldF;
typedef MACGrid<scalarF,scalarF> MACVelocityFieldF;

typedef Grid<unsigned char,scalarD> TagFieldD;
typedef Grid<scalarD,scalarD> ScalarFieldD;
typedef Grid<Vec3d,scalarD> VectorFieldD;
typedef MACGrid<unsigned char,scalarD> MACTagFieldD;
typedef MACGrid<scalarD,scalarD> MACVelocityFieldD;

typedef Grid<unsigned char,scalar> TagField;
typedef Grid<scalar,scalar> ScalarField;
typedef Grid<Vec3,scalar> VectorField;
typedef MACGrid<unsigned char,scalar> MACTagField;
typedef MACGrid<scalar,scalar> MACVelocityField;

PRJ_END

#endif
