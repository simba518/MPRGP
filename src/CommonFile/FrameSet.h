#ifndef FRAME_SET_IO
#define FRAME_SET_IO

#include "ParticleSet.h"
#include "ParticleCD.h"

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>

PRJ_BEGIN

template <typename T>
struct HasId : boost::false_type {};
template <typename T>
struct HasType : boost::false_type {};

template <typename A,typename B,bool has>
struct CopyId {
    static FORCE_INLINE void copy(A& a,const B& b) {
        a._id=b._id;
    }
};
template <typename A,typename B>
struct CopyId<A,B,false> {
    static FORCE_INLINE void copy(A& a,const B& b) {}
};

template <typename A,typename B,bool has>
struct CopyType {
    static FORCE_INLINE void copy(A& a,const B& b) {
        a._type=b._type;
    }
};
template <typename A,typename B>
struct CopyType<A,B,false> {
    static FORCE_INLINE void copy(A& a,const B& b) {}
};

enum PARTICLE_TYPE {
    FREE_SURFACE=1,
    FREE_SURFACE_SAMPLED=2,
    SOLID_BOUNDARY=4,
    SOLID_BOUNDARY_SAMPLED=8,
    MEDIAL_AXIS=16,
    MEDIAL_AXIS_SAMPLED=32,
    OTHER=64,

    VALID=128,
};
template<typename T>
struct ParticleSetTplWithTag : public ParticleSetTpl<T> {
public:
    void setNrFS(const sizeType& nr) {
        _nrParticleFS=nr;
    }
    sizeType getBegin(const sizeType& type) const {
        return _typeMap[type].first;
    }
    sizeType getEnd(const sizeType& type) {
        return _typeMap[type].second;
    }
    ParticleSetTplWithTag<T> getSubset(const sizeType& type) const {
        ParticleSetTplWithTag<T> ret;
        for(sizeType i=0; i<size(); i++)
            if(get(i)._type&type)
                ret.addParticle(get(i));
        return ret;
    }
    void debug(const string& path,const sizeType& i) const {
        {
            ostringstream oss;
            oss << path << "/PSET_" << i << ".vtk";
            writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> fs=getSubset(FREE_SURFACE);
            ostringstream oss;
            oss << path << "/FREE_SURFACE_" << i << ".vtk";
            fs.writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> fss=getSubset(FREE_SURFACE_SAMPLED);
            ostringstream oss;
            oss << path << "/FREE_SURFACE_SAMPLED_" << i << ".vtk";
            fss.writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> b=getSubset(SOLID_BOUNDARY);
            ostringstream oss;
            oss << path << "/SOLID_BOUNDARY_" << i << ".vtk";
            b.writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> bs=getSubset(SOLID_BOUNDARY_SAMPLED);
            ostringstream oss;
            oss << path << "/SOLID_BOUNDARY_SAMPLED_" << i << ".vtk";
            bs.writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> ma=getSubset(MEDIAL_AXIS);
            ostringstream oss;
            oss << path << "/MEDIAL_AXIS_" << i << ".vtk";
            ma.writeVTK(oss.str());
        }

        {
            ParticleSetTplWithTag<T> mas=getSubset(MEDIAL_AXIS_SAMPLED);
            ostringstream oss;
            oss << path << "/MEDIAL_AXIS_SAMPLED_ " << i << ".vtk";
            mas.writeVTK(oss.str());
        }
    }
};
template <typename T>
struct ParticleWithTag : public ParticleBase<T> {
    ALIGN_16 sizeType _id;
    ALIGN_16 sizeType _type;
    void set(sizeType type) {
        _type|=type;
    }
    void clear(sizeType type) {
        _type&=~type;
    }
    ParticleWithTag():_id(0),_type(0) {}
    template<typename OTHER>
    ParticleWithTag& copy(const OTHER& p) {
        ParticleBase::copy(p);
        CopyId<ParticleWithTag,OTHER,HasId<OTHER>::value>::copy(*this,p);
        CopyType<ParticleWithTag,OTHER,HasType<OTHER>::value>::copy(*this,p);
        return *this;
    }
    virtual bool read(istream& is) {
        return ParticleBase::read(is) && readBinaryData(_id,is).good() && readBinaryData(_type,is).good();
    }
    virtual bool write(ostream& os) const {
        return ParticleBase::write(os) && writeBinaryData(_id,os).good() && writeBinaryData(_type,os).good();
    }
};
template <typename T>
struct ParticleWithVelTag : public Particle<T> {
    ALIGN_16 sizeType _id;
    ALIGN_16 sizeType _type;
    void set(sizeType type) {
        _type|=type;
    }
    void clear(sizeType type) {
        _type&=~type;
    }
    ParticleWithVelTag():_id(0),_type(0) {}
    template<typename OTHER>
    ParticleWithVelTag& copy(const OTHER& p) {
        ParticleBase::copy(p);
        CopyId<ParticleWithVelTag,OTHER,HasId<OTHER>::value>::copy(*this,p);
        CopyType<ParticleWithVelTag,OTHER,HasType<OTHER>::value>::copy(*this,p);
        return *this;
    }
    virtual bool read(istream& is) {
        return Particle::read(is) && readBinaryData(_id,is).good() && readBinaryData(_type,is).good();
    }
    virtual bool write(ostream& os) const {
        return Particle::write(os) && writeBinaryData(_id,os).good() && writeBinaryData(_type,os).good();
    }
};

template <typename T,typename TI,typename TV=vector<T,Eigen::aligned_allocator<T> > >
class LevelSetFrameData
{
public:
    typedef Grid<T,TI,TV> GridType;
    typedef MACGrid<T,TI,TV> MACGridType;
    typedef T scalarType;
    bool read(istream& is) {
        return _phi.read(is) && _vel.read(is);
    }
    bool write(ostream& os) const {
        return _phi.write(os) && _vel.write(os);
    }
    GridType _phi;
    MACGridType _vel;
};
template <typename PS_TYPE>
struct FrameSetHeader {
    FrameSetHeader() {
        _nFrames=0;
        _idFrom=-1;
        _idTo=-1;

        _vTimeFrom=0.0f;
        _vTimeTo=0.0f;
    }
    sizeType _nFrames;
    sizeType _idFrom;
    sizeType _idTo;
    typename PS_TYPE::scalarType _vTimeFrom;
    typename PS_TYPE::scalarType _vTimeTo;
};
template <typename PS_TYPE>
class DefaultWriter
{
public:
    bool read(PS_TYPE& ps,istream& is) const {
        return ps.read(is);
    }
    bool write(const PS_TYPE& ps,ostream& os) const {
        return ps.write(os);
    }
};
template <typename PS_TYPE,bool hasVel=true>
class ParticleSetCompressor
{
    //a particle set quantization lib
public:
    typedef typename PS_TYPE::ParticleType P_TYPE;
    typedef typename P_TYPE::scalarType T;
    typedef typename ScalarUtil<T>::ScalarVec3 PT3;
    template <bool HAS_VEL_TRAIT> struct VelOp {
        static void writeVel(const P_TYPE& p,unsigned int& dat,ostream& os,const T& maxVel) {
            ParticleSetCompressor::compressVel(p._vel,dat,maxVel);
            writeBinaryData(dat,os);
        }
        static void readVel(P_TYPE& p,unsigned int& dat,istream& is,const T& maxMag,bool is2D) {
            readBinaryData(dat,is);
            ParticleSetCompressor::uncompressVel(dat,p._vel,maxMag);
            if(is2D)
                p._vel.z()=0.0f;
        }
        static PT3 getVel(const P_TYPE& p) {
            return p._vel;
        }
    };
    template <> struct VelOp<false> {
        static void writeVel(const P_TYPE& p,unsigned int& dat,ostream& os,const T& maxVel) {}
        static void readVel(P_TYPE& p,unsigned int& dat,istream& is,const T& maxMag,bool is2D) {}
        static PT3 getVel(const P_TYPE& p) {
            return PT3();
        }
    };
public:
    void reset(const BBox<T>& bb,const PT3& cellSz) {
        _cd.reset(new CollisionInterface<T,PS_TYPE>());
        _cd->reset(bb,cellSz);
    }
    static void compressVel(const PT3& vel,unsigned int& dat,const T& maxMag) {
        const T angX=atan2(vel.y(),vel.x())+M_PI;
        const T angY=atan2(vel.z(),vel.block<2,1>(0,0).norm())+M_PI;
        T mag=vel.norm()/maxMag;
        mag=std::min<T>(std::max<T>(mag,0.0f),1.0f);
        dat=(((unsigned int)(mag*65535)) << 16)+
            (((unsigned int)((angX/(M_PI*2.0f))*255)) << 8)+
            (((unsigned int)((angY/(M_PI*2.0f))*255)) << 0);
    }
    static void uncompressVel(const unsigned int& dat,PT3& vel,const T& maxMag) {
        const T mag=(((dat>>16)&65535)/65535.0f)*maxMag;
        const T angX=(((dat>>8)&255)/255.0f)*(M_PI*2.0f)-M_PI;
        const T angY=(((dat>>0)&255)/255.0f)*(M_PI*2.0f)-M_PI;
        vel.z()=sin(angY);
        vel.x()=cos(angY);
        vel.y()=vel.x()*sin(angX);
        vel.x()*=cos(angX);
        vel*=mag;
    }
    static void compressPos(const PT3& pos,unsigned int& dat,const PT3& base,const PT3& cellSz) {
        const unsigned int res=(1<<10)-1;
        PT3 tmp=((pos-base).array()/cellSz.array()).matrix();
        tmp.x()=std::min<T>(std::max<T>(tmp.x(),0.0f),1.0f);
        tmp.y()=std::min<T>(std::max<T>(tmp.y(),0.0f),1.0f);
        tmp.z()=std::min<T>(std::max<T>(tmp.z(),0.0f),1.0f);
        dat=(((unsigned int)(tmp.x()*res)) << (unsigned int)20)+
            (((unsigned int)(tmp.y()*res)) << (unsigned int)10)+
            (((unsigned int)(tmp.z()*res)) << (unsigned int)0);
    }
    static void uncompressPos(const unsigned int& dat,PT3& pos,const PT3& base,const PT3& cellSz) {
        const unsigned int res=(1<<10)-1;
        pos.x()=((dat>>(unsigned int)20)&res)*cellSz.x()/(T)res+base.x();
        pos.y()=((dat>>(unsigned int)10)&res)*cellSz.y()/(T)res+base.y();
        pos.z()=((dat>>(unsigned int)0)&res)*cellSz.z()/(T)res+base.z();
    }
    bool read(PS_TYPE& ps,istream& is) const {
        if(!ps.HasMagic::readMagic(is))
            return false;

        //resize particle count
        sizeType n;
        readBinaryData(n,is);
        ps.resize(n);
        //read maximum velocity length
        T maxMag;
        readBinaryData(maxMag,is);
        //read grid count in XYZ
        Vec3i nrGrid;
        readBinaryData(nrGrid,is);
        //read cell edge length
        PT3 cellSz;
        readBinaryData(cellSz,is);
        //read left-bottom-front corner
        PT3 base;
        readBinaryData(base,is);
        //reconstruct
        //store all the tags for velocity
        sizeType tagId=0;
        std::vector<unsigned char> tags(nrGrid.prod());
        //now read all position
        sizeType total=0;
        for(sizeType x=0; x<nrGrid.x(); x++)
            for(sizeType y=0; y<nrGrid.y(); y++) {
                for(sizeType z=0; z<nrGrid.z(); z++) {
                    //read how many
                    unsigned char nr;
                    readBinaryData(nr,is);
                    tags[tagId++]=nr;
                    //read particles
                    PT3 baseCell=base+PT3(cellSz.x()*(T)x,cellSz.y()*(T)y,cellSz.z()*(T)z);
                    unsigned int dat;
                    for(sizeType i=0; i<nr; i++,total++) {
                        P_TYPE& p=ps.get(total);
                        readBinaryData(dat,is);
                        uncompressPos(dat,p._pos,baseCell,cellSz);
                        VelOp<hasVel>::readVel(p,dat,is,maxMag,nrGrid.z() <= 1);
                        if(nrGrid.z() <= 1)
                            p._pos.z()=0.0f;
                    }
                }
            }
        ps.resize(total);
        return is.good();
    }
    bool write(const PS_TYPE& ps,ostream& os) const {
        if(!ps.HasMagic::writeMagic(os))
            return false;

        //prepare collision detection
        _cd->prepare(ps);
        const Grid<pair<sizeType,sizeType>,T>& ssGrid=_cd->getStartEndGrid();
        const std::vector<sizeType,Eigen::aligned_allocator<sizeType> >& ssIndex=_cd->getIndex();

        //collect max velocity
        T maxVel=0.0f;
        for(sizeType i=0; i<ps.size(); i++)
            maxVel=std::max<T>(maxVel,VelOp<hasVel>::getVel(ps.get(i)).squaredNorm());
        maxVel=sqrt(maxVel+ScalarUtil<T>::scalar_eps);

        //write header
        //resize particle count
        sizeType n=ps.size();
        writeBinaryData(n,os);
        //read maximum velocity length
        writeBinaryData(maxVel,os);
        //read grid count in XYZ
        Vec3i nrGrid=ssGrid.getNrCell();
        writeBinaryData(nrGrid,os);
        //read cell edge length
        PT3 cellSz=ssGrid.getCellSize();
        if(cellSz.z() < ScalarUtil<T>::scalar_eps)
            cellSz.z()=1.0f;
        writeBinaryData(cellSz,os);
        //read left-bottom-front corner
        PT3 base=ssGrid.getPt(Vec3i::Zero());
        writeBinaryData(base,os);

        //write all positions
        sizeType nrTotal=0;
        std::vector<unsigned short> tags;
        for(sizeType x=0; x<nrGrid.x(); x++)
            for(sizeType y=0; y<nrGrid.y(); y++) {
                for(sizeType z=0; z<nrGrid.z(); z++) {
                    const std::pair<sizeType,sizeType>& ss=ssGrid.get(Vec3i(x,y,z));
                    unsigned char nrP=(unsigned char)std::min<sizeType>(ss.second-ss.first,255);
                    writeBinaryData(nrP,os);
                    nrTotal+=nrP;

                    PT3 baseC=base+PT3(cellSz.x()*(T)x,cellSz.y()*(T)y,cellSz.z()*(T)z);
                    unsigned int dat;
                    for(sizeType i=0; i<nrP; i++) {
                        const P_TYPE& p=ps.get(ssIndex[ss.first+i]);
                        compressPos(p._pos,dat,baseC,cellSz);
                        writeBinaryData(dat,os);
                        VelOp<hasVel>::writeVel(p,dat,os,maxVel);
                    }
                }
            }
        INFOV("PSet Size: %d, Total Written: %d",ps.size(),nrTotal)
        return os.good();
    }
protected:
    boost::shared_ptr<CollisionInterface<T,PS_TYPE> > _cd;
};
template <typename PS_TYPE,typename WRITER=DefaultWriter<PS_TYPE> >
struct FrameSetTpl : public HasMagic {
#define USE_ZLIB

    typedef typename PS_TYPE::scalarType scalarType;
    struct Frame {
        Frame():_isProxy(false) {}
        //data
        PS_TYPE _pSet;
        scalarType _vTime;
        sizeType _id;
        //does this frame contain data
        bool _isProxy;
    };
    FrameSetTpl():HasMagic(0xABCDABCDABCDABCD) {}
    virtual ~FrameSetTpl() {}
    //getter
    bool empty() const {
        return _frames.empty() && _framesFastAccess.empty();
    }
    sizeType size() const {
        return _header._idTo-_header._idFrom+1;
        //return max(_frames.size(),_framesFastAccess.size());
    }
    scalarType from() const {
        return _header._vTimeFrom;
    }
    scalarType to() const {
        return _header._vTimeTo;
    }
    sizeType idFrom() const {
        return _header._idFrom;
    }
    sizeType idTo() const {
        return _header._idTo;
    }
    //setter
    void appendFrame(Frame frm,bool headerOnly=false) {
        ASSERT(frm._vTime >= _header._vTimeTo);
        if(idFrom() != -1)
            ASSERT(frm._id == _header._idTo+1);

        optimizeAccess(false);
        _frames.push_back(frm);
        if(headerOnly) {
            _frames.back()._pSet.clear();
            _frames.back()._isProxy=true;
        } else {
            ASSERT(frm._isProxy == false)
        }
        syncHeader();
    }
    void reset() {
        _frames.swap(std::list<Frame>());
        _framesFastAccess.swap(std::vector<Frame>());
        syncHeader();
    }
    void discardSince(const sizeType& id) {
        if(id < idFrom()) {
            reset();
            return;
        }
        if(id > idTo()) {
            return;
        }

        sizeType offset=id-idFrom();
        if(_frames.empty()) {
            _framesFastAccess.erase(_framesFastAccess.begin()+offset,_framesFastAccess.end());
        } else {
            std::list<Frame>::const_iterator iter=_frames.begin();
            for(sizeType j=0; j<offset; j++)
                iter++;
            _frames.erase(iter,_frames.end());
        }

        syncHeader();
    }
    bool getSubset(FrameSetTpl& ret,const sizeType& fFrom,const sizeType& fTo,bool headerOnly=false) {
        if(fFrom < _header._idFrom || fTo > _header._idTo)
            return false;
        if(fTo < _header._idFrom || fTo > _header._idTo)
            return false;
        if(fTo < fFrom)
            return false;

        optimizeAccess(true);

        ret.reset();
        for(sizeType i=fFrom-idFrom(); i<=fTo-idFrom(); i++)
            ret.appendFrame(getFrame(i),headerOnly);
        return true;
    }
    //io
    template <typename PS_TYPE2,typename WRITER2>
    bool readBinary(istream& is,const sizeType& interval=1,FrameSetTpl<PS_TYPE2,WRITER2>* to=NULL,bool useWriter=true) {
        if(!HasMagic::readMagic(is))
            return false;

        {
            reset();

            readBinaryData(_header._nFrames,is);
            readBinaryData(_header._idFrom,is);
            readBinaryData(_header._idTo,is);
            readBinaryData(_header._vTimeFrom,is);
            readBinaryData(_header._vTimeTo,is);

            if(to) {
                to->reset();
                to->_header._nFrames=_header._nFrames;
                to->_header._idFrom=_header._idFrom;
                to->_header._idTo=_header._idTo;
                to->_header._vTimeFrom=_header._vTimeFrom;
                to->_header._vTimeTo=_header._vTimeTo;
            }
        }

        {
            PS_TYPE tmp;
            if(to)
                to->_framesFastAccess.resize(_header._nFrames);
            else
                _framesFastAccess.resize(_header._nFrames);
            for(sizeType i=0; i<(sizeType)_header._nFrames; i++) {
                if(i%50 == 0)
                    INFOV("Read %d Frames",i)

                    if(!is.good())
                        return false;

                Frame frm;
                readBinaryData(frm._vTime,is);
                readBinaryData(frm._id,is);
                if(i%interval == 0) {
                    if(useWriter) _writer.read(frm._pSet,is);
                    else frm._pSet.read(is);
                    frm._isProxy=false;
                } else {
                    if(useWriter) _writer.read(tmp,is);
                    else tmp.read(is);
                    frm._pSet.clear();
                    frm._isProxy=true;
                }

                if(to) {
                    to->_framesFastAccess[i]._id=frm._id;
                    to->_framesFastAccess[i]._vTime=frm._vTime;
                    to->_framesFastAccess[i]._isProxy=frm._isProxy;
                    to->_framesFastAccess[i]._pSet.copy(frm._pSet);
                } else {
                    _framesFastAccess[i]=frm;
                }
            }
        }

        if(to)
            to->syncHeader();
        else
            syncHeader();
        return true;
    }
    bool writeBinary(ostream& os,bool proxy=false,bool useWriter=true) const {
        typedef std::list<Frame>::const_iterator iterL;
        typedef std::vector<Frame>::const_iterator iterV;

        //check no proxy
        if(!proxy) {
            if(!_frames.empty()) {
                for(iterL beg=_frames.begin(),end=_frames.end(); beg!=end; beg++)
                    if(beg->_isProxy) {
                        INFO("You are writing to a frameSet with proxy!")
                        return false;
                    }
            } else {
                for(iterV beg=_framesFastAccess.begin(),end=_framesFastAccess.end(); beg!=end; beg++)
                    if(beg->_isProxy) {
                        INFO("You are writing to a frameSet with proxy!")
                        return false;
                    }
            }
        }

        //write header
        {
            if(!HasMagic::writeMagic(os))
                return false;

            writeBinaryData(_header._nFrames,os);
            writeBinaryData(_header._idFrom,os);
            writeBinaryData(_header._idTo,os);
            writeBinaryData(_header._vTimeFrom,os);
            writeBinaryData(_header._vTimeTo,os);
        }

        //write frameSet data
        {
            PS_TYPE emptyPSet;
            sizeType id=0;
            if(!_frames.empty()) {
                for(iterL beg=_frames.begin(),end=_frames.end(); beg!=end; beg++) {
                    writeBinaryData(beg->_vTime,os);
                    writeBinaryData(beg->_id,os);
                    if(!proxy) _writer.write(beg->_pSet,os);
                    else emptyPSet.write(os);
                }
            } else {
                for(iterV beg=_framesFastAccess.begin(),end=_framesFastAccess.end(); beg!=end; beg++) {
                    writeBinaryData(beg->_vTime,os);
                    writeBinaryData(beg->_id,os);
                    if(!proxy) _writer.write(beg->_pSet,os);
                    else emptyPSet.write(os);
                }
            }
        }

        //syncHeader();
        return true;
    }
    template <typename PS_TYPE2,typename WRITER2>
    bool readBinarySepFile(const boost::filesystem::path& path,const sizeType& interval=1,FrameSetTpl<PS_TYPE2,WRITER2>* to=NULL,bool useWriter=true) {
        if(!boost::filesystem::exists(path) || !boost::filesystem::is_directory(path))
            return false;

        {
            ostringstream oss;
            oss << "header.dat";
            boost::filesystem::ifstream is(path/oss.str(),ios::binary);
            if(!HasMagic::readMagic(is))
                return false;

            reset();

            readBinaryData(_header._nFrames,is);
            readBinaryData(_header._idFrom,is);
            readBinaryData(_header._idTo,is);
            readBinaryData(_header._vTimeFrom,is);
            readBinaryData(_header._vTimeTo,is);
            ASSERT(_header._nFrames == _header._idTo-_header._idFrom+1)

            if(to) {
                to->reset();
                to->_header._nFrames=_header._nFrames;
                to->_header._idFrom=_header._idFrom;
                to->_header._idTo=_header._idTo;
                to->_header._vTimeFrom=_header._vTimeFrom;
                to->_header._vTimeTo=_header._vTimeTo;
            }
        }

        {
            if(to)
                to->_framesFastAccess.resize(_header._nFrames);
            else
                _framesFastAccess.resize(_header._nFrames);
            for(sizeType i=0; i<(sizeType)_header._nFrames; i++) {
                if(i%50 == 0)
                    INFOV("Read %d Frames",i)

                    ostringstream oss;
                oss << "f" << (i+_header._idFrom) << ".dat";

                Frame frm;
                if(!readFrameSepFile(frm,(path/oss.str()).string(),i%interval != 0))
                    return false;

                if(to) {
                    to->_framesFastAccess[i]._id=frm._id;
                    to->_framesFastAccess[i]._vTime=frm._vTime;
                    to->_framesFastAccess[i]._isProxy=frm._isProxy;
                    to->_framesFastAccess[i]._pSet.copy(frm._pSet);
                } else {
                    _framesFastAccess[i]=frm;
                }
            }

            if(to)
                to->syncHeader();
            else
                syncHeader();
        }
        return true;
    }
    bool writeBinarySepFile(const boost::filesystem::path& path,bool proxy=false,bool useWriter=true) const {
        typedef std::list<Frame>::const_iterator iterL;
        typedef std::vector<Frame>::const_iterator iterV;

        //check path
        if(!boost::filesystem::exists(path) || !boost::filesystem::is_directory(path))
            return false;

        //check no proxy
        if(!proxy) {
            if(!_frames.empty()) {
                for(iterL beg=_frames.begin(),end=_frames.end(); beg!=end; beg++)
                    if(beg->_isProxy) {
                        INFO("You are writing to a frameSet with proxy!")
                        return false;
                    }
            } else {
                for(iterV beg=_framesFastAccess.begin(),end=_framesFastAccess.end(); beg!=end; beg++)
                    if(beg->_isProxy) {
                        INFO("You are writing to a frameSet with proxy!")
                        return false;
                    }
            }
        }

        //write header
        {
            ostringstream oss;
            oss << "header.dat";
            boost::filesystem::ofstream os(path/oss.str(),ios::binary);
            if(!HasMagic::writeMagic(os))
                return false;

            writeBinaryData(_header._nFrames,os);
            writeBinaryData(_header._idFrom,os);
            writeBinaryData(_header._idTo,os);
            writeBinaryData(_header._vTimeFrom,os);
            writeBinaryData(_header._vTimeTo,os);
        }

#define CREATE_OS boost::filesystem::ofstream os(path/oss.str(),ios::binary);
#define CREATE_OS_ZLIB														\
boost::iostreams::filtering_ostream os;										\
os.push(boost::iostreams::bzip2_compressor());								\
os.push(boost::iostreams::file_sink((path/oss.str()).string(),ios::binary));
#define WRITE_FRM															\
if(!HasMagic(0xABCDABCDABCDAABB).writeMagic(os))							\
	return false;															\
writeBinaryData(beg->_vTime,os);											\
writeBinaryData(beg->_id,os);												\
if(!proxy){																	\
	if(useWriter) _writer.write(beg->_pSet,os);								\
	else beg->_pSet.write(os);												\
}else{																		\
	if(useWriter) _writer.write(emptyPSet,os);								\
	else emptyPSet.write(os);												\
}
#define WRITE_FRAME					\
{									\
	ostringstream oss;				\
	oss << "f" << id++ << ".dat";	\
	CREATE_OS						\
	WRITE_FRM						\
}
#define WRITE_FRAME_ZLIB			\
{									\
	ostringstream oss;				\
	oss << "f" << id++ << ".dat";	\
	CREATE_OS_ZLIB					\
	WRITE_FRM						\
}

        //write data
        PS_TYPE emptyPSet;
        sizeType id=_header._idFrom;
        if(!_frames.empty()) {
            for(iterL beg=_frames.begin(),end=_frames.end(); beg!=end; beg++)
#ifdef USE_ZLIB
                WRITE_FRAME_ZLIB
#else
                WRITE_FRAME
#endif
            } else {
            for(iterV beg=_framesFastAccess.begin(),end=_framesFastAccess.end(); beg!=end; beg++)
#ifdef USE_ZLIB
                WRITE_FRAME_ZLIB
#else
                WRITE_FRAME
#endif
            }

        //syncHeader();
        return true;
#undef CREATE_OS
#undef CREATE_OS_ZLIB
#undef WRITE_FRM
#undef WRITE_FRAME
#undef WRITE_FRAME_ZLIB
    }
    bool readFrameSepFile(Frame& frm,const std::string& path,bool proxy=false,bool useWriter=true) {
#ifdef USE_ZLIB
        boost::iostreams::filtering_istream is;
        is.push(boost::iostreams::bzip2_decompressor());
        is.push(boost::iostreams::file_source(path,ios::binary));
#else
        boost::filesystem::ifstream is(path,ios::binary);
#endif
        if(!HasMagic(0xABCDABCDABCDAABB).readMagic(is))
            return false;

        readBinaryData(frm._vTime,is);
        readBinaryData(frm._id,is);
        if(!proxy) {
            if(useWriter)
                _writer.read(frm._pSet,is);
            else
                frm._pSet.read(is);
            frm._isProxy=false;
        } else {
            //since this is the end of file,
            //we don't need to read this.
            //tmp.read(is);
            frm._pSet.clear();
            frm._isProxy=true;
        }

        return true;
    }
    void writeVTK(const string& path) const {
        boost::filesystem::create_directory(path);
        if(_frames.empty()) {
            sizeType n=size();
            for(sizeType i=0; i<n; i++) {
                ostringstream oss;
                oss << path << "frame_" << i << ".vtk";
                _framesFastAccess[i]._pSet.writeVTK(oss.str());
            }
        } else {
            typedef std::list<Frame>::const_iterator iter;
            sizeType i=0;
            for(iter beg=_frames.begin(),end=_frames.end(); beg!=end; beg++) {
                ostringstream oss;
                oss << path << "frame_" << i << ".vtk";
                beg->_pSet.writeVTK(oss.str());
                i++;
            }
        }
    }
    void writeVTKSampled(const string& path) const {
        boost::filesystem::create_directory(path);
        if(_frames.empty()) {
            sizeType n=size();
            for(sizeType i=0; i<n; i++)
                _framesFastAccess[i]._pSet.debug(path,i);
        } else {
            typedef std::list<Frame>::const_iterator iter;
            sizeType i=0;
            for(iter beg=_frames.begin(),end=_frames.end(); beg!=end; beg++)
                beg->_pSet.debug(path,i);
        }
    }
    //memory access
    Frame& getFrame(const sizeType& id) {
        ASSERT(id >= idFrom() && id <= idTo())

        sizeType offset=id-idFrom();
        if(!_framesFastAccess.empty())
            return _framesFastAccess[offset];
        else {
            std::list<Frame>::iterator iter=_frames.begin();
            for(sizeType j=0; j<offset; j++)
                iter++;
            return *iter;
        }
    }
    Frame& getFrameSafe(const sizeType& id) {
        ASSERT(id >= idFrom() && id <= idTo())

        sizeType offset=id-idFrom();
        if(!_framesFastAccess.empty()) {
            ASSERT(!_framesFastAccess[0]._isProxy)
            while(_framesFastAccess[offset]._isProxy)
                offset--;
            return _framesFastAccess[offset];
        } else {
            ASSERT(!_frames.front()._isProxy)
            std::list<Frame>::iterator iter=_frames.begin();
            for(sizeType j=0; j<offset; j++)
                iter++;
            while(iter->_isProxy)
                iter--;
            return *iter;
        }
    }
    const Frame& getFrame(const sizeType& id) const {
        return (const Frame&)(const_cast<FrameSetTpl<PS_TYPE>&>(*this).getFrame(id));
    }
    void optimizeAccess(bool fastAccess) {
        if(fastAccess) {
            if(_framesFastAccess.empty()) {
                std::insert_iterator<std::vector<Frame> > inserter(_framesFastAccess,_framesFastAccess.end());
                std::copy(_frames.begin(),_frames.end(),inserter);
                _frames.swap(std::list<Frame>());
            }
        } else {
            if(_frames.empty()) {
                std::insert_iterator<std::list<Frame> > inserter(_frames,_frames.end());
                std::copy(_framesFastAccess.begin(),_framesFastAccess.end(),inserter);
                _framesFastAccess.swap(std::vector<Frame>());
            }
        }
    }
    void syncHeader() {
        if(!_frames.empty()) {
            _header._nFrames=_frames.size();
            if(_header._nFrames > 0) {
                _header._idFrom=_frames.front()._id;
                _header._idTo=_frames.back()._id;
                _header._vTimeFrom=_frames.front()._vTime;
                _header._vTimeTo=_frames.back()._vTime;
            } else {
                _header=FrameSetHeader<PS_TYPE>();
            }
        } else {
            _header._nFrames=_framesFastAccess.size();
            if(_header._nFrames > 0) {
                _header._idFrom=_framesFastAccess.front()._id;
                _header._idTo=_framesFastAccess.back()._id;
                _header._vTimeFrom=_framesFastAccess.front()._vTime;
                _header._vTimeTo=_framesFastAccess.back()._vTime;
            } else {
                _header=FrameSetHeader<PS_TYPE>();
            }
        }
    }
    FrameSetHeader<PS_TYPE> getHeader() const {
        return _header;
    }
    template <typename PS_TYPE_IN>
    void copy(FrameSetTpl<PS_TYPE_IN>& other) {
        reset();

        _header._nFrames=other._header._nFrames;
        _header._idFrom=other._header._idFrom;
        _header._idTo=other._header._idTo;
        _header._vTimeFrom=(PS_TYPE::scalarType)other._header._vTimeFrom;
        _header._vTimeTo=(PS_TYPE::scalarType)other._header._vTimeTo;

        other.optimizeAccess(true);
        _framesFastAccess.resize(other._framesFastAccess.size());
        for(sizeType i=0; i<(sizeType)_framesFastAccess.size(); i++) {
            _framesFastAccess[i]._pSet.copy(other._framesFastAccess[i]._pSet);
            _framesFastAccess[i]._vTime=(PS_TYPE::scalarType)other._framesFastAccess[i]._vTime;
            _framesFastAccess[i]._id=other._framesFastAccess[i]._id;
        }
    }
    //random access and modify
    bool readHeaderR(const boost::filesystem::path& path,FrameSetHeader<PS_TYPE>& header) const {
        if(!boost::filesystem::exists(path) || !boost::filesystem::is_directory(path))
            return false;

        ostringstream oss;
        oss << "header.dat";
        boost::filesystem::ifstream is(path/oss.str(),ios::binary);
        if(!HasMagic(0xABCDABCDABCDABCD).readMagic(is))
            return false;

        readBinaryData(header._nFrames,is);
        readBinaryData(header._idFrom,is);
        readBinaryData(header._idTo,is);
        readBinaryData(header._vTimeFrom,is);
        readBinaryData(header._vTimeTo,is);
        return true;
    }
    bool readFrameR(const boost::filesystem::path& path,Frame& frm,bool useWriter=true) const {
        ostringstream oss;
        oss << "f" << frm._id << ".dat";
#ifdef USE_ZLIB
        boost::iostreams::filtering_istream is;
        is.push(boost::iostreams::bzip2_decompressor());
        is.push(boost::iostreams::file_source((path/oss.str()).string(),ios::binary));
#else
        boost::filesystem::ifstream is(path/oss.str(),ios::binary);
#endif
        return readFrm(frm,is,useWriter);
    }
    bool setHeaderR(const boost::filesystem::path& path,const FrameSetHeader<PS_TYPE>& header) const {
        if(!boost::filesystem::exists(path) || !boost::filesystem::is_directory(path))
            return false;

        ASSERT(header._nFrames >= 0)
        if(header._nFrames == 0)
            ASSERT(header._idFrom == -1 && header._idTo == -1)
            else
                ASSERT(header._nFrames == header._idTo-header._idFrom+1)

            {
                ostringstream oss;
                oss << "header.dat";
                boost::filesystem::ofstream os(path/oss.str(),ios::binary|ios::trunc);
                if(!HasMagic::writeMagic(os))
                    return false;

                writeBinaryData(header._nFrames,os);
                writeBinaryData(header._idFrom,os);
                writeBinaryData(header._idTo,os);
                writeBinaryData(header._vTimeFrom,os);
                writeBinaryData(header._vTimeTo,os);
            }
        return true;
    }
    bool appendFrameR(const boost::filesystem::path& path,const Frame& frm,bool useWriter=true,bool checkHeader=true) const {
        ASSERT(frm._isProxy == false)

        if(checkHeader) {
            FrameSetHeader<PS_TYPE> header;
            if(!readHeaderR(path,header))
                return false;

            ASSERT(frm._id == header._idTo+1 || header._idTo  == -1)
            if(header._idTo == -1) {
                header._nFrames++;
                header._idFrom=header._idTo=frm._id;
                header._vTimeFrom=header._vTimeTo=frm._vTime;
            } else {
                header._nFrames++;
                header._idTo=frm._id;
                header._vTimeTo=frm._vTime;
            }
            ASSERT(header._nFrames == header._idTo-header._idFrom+1)
            if(!setHeaderR(path,header))
                return false;
        }

        {
            ostringstream oss;
            oss << "f" << frm._id << ".dat";
#ifdef USE_ZLIB
            boost::iostreams::filtering_ostream os;
            os.push(boost::iostreams::bzip2_compressor());
            os.push(boost::iostreams::file_sink((path/oss.str()).string(),ios::binary));
#else
            boost::filesystem::ofstream os(path/oss.str(),ios::binary);
#endif
            return writeFrm(frm,os,useWriter);
        }
    }
    bool discardSinceR(const boost::filesystem::path& path,const sizeType& id) const {
        FrameSetHeader<PS_TYPE> header;
        if(!readHeaderR(path,header))
            return false;

        if(id <= header._idFrom) {
            header=FrameSetHeader<PS_TYPE>();
            return setHeaderR(path,header);
        } else if(id > header._idTo) {
            return true;
        } else {
            Frame frm;
            frm._id=id;
            if(!readFrameR(path,frm))
                return false;

            header._idTo=frm._id-1;
            header._vTimeTo=frm._vTime;
            header._nFrames=header._idTo-header._idFrom+1;
            return setHeaderR(path,header);
        }
    }
    bool createR(const boost::filesystem::path& path) const {
        boost::filesystem::create_directory(path);
        return setHeaderR(path,FrameSetHeader<PS_TYPE>());
    }
    bool changeDataR(const boost::filesystem::path& path,const Frame& frm,bool useWriter=true,bool checkHeader=true) const {
        ASSERT(frm._isProxy == false)

        if(checkHeader) {
            FrameSetHeader<PS_TYPE> header;
            if(!readHeaderR(path,header))
                return false;

            if(frm._id < header._idFrom || frm._id > header._idTo)
                return false;
        }

        {
            ostringstream oss;
            oss << "f" << frm._id << ".dat";
#ifdef USE_ZLIB
            boost::iostreams::filtering_ostream os;
            os.push(boost::iostreams::bzip2_compressor());
            os.push(boost::iostreams::file_sink((path/oss.str()).string(),ios::binary|ios::trunc));
#else
            boost::filesystem::ofstream os(path/oss.str(),ios::binary|ios::trunc);
#endif
            return writeFrm(frm,os,useWriter);
        }
    }
    bool writeFrm(const Frame& frm,ostream& os,bool useWriter) const {
        if(!HasMagic(0xABCDABCDABCDAABB).writeMagic(os))
            return false;

        writeBinaryData(frm._vTime,os);
        writeBinaryData(frm._id,os);
        if(useWriter) _writer.write(frm._pSet,os);
        else frm._pSet.write(os);
        return true;
    }
    bool readFrm(Frame& frm,istream& is,bool useWriter) const {
        if(!HasMagic(0xABCDABCDABCDAABB).readMagic(is))
            return false;

        readBinaryData(frm._vTime,is);
        readBinaryData(frm._id,is);
        if(useWriter) _writer.read(frm._pSet,is);
        else frm._pSet.read(is);
        return true;
    }

    FrameSetHeader<PS_TYPE> _header;
    std::list<Frame> _frames;
    std::vector<Frame> _framesFastAccess;
    WRITER _writer;

#undef USE_ZLIB
};

typedef ParticleSetTplWithTag<ParticleWithTag<scalarF> > ParticleSetTF;
typedef ParticleSetTplWithTag<ParticleWithVelTag<scalarF> > ParticleSetVTF;

typedef ParticleSetTplWithTag<ParticleWithTag<scalarD> > ParticleSetTD;
typedef ParticleSetTplWithTag<ParticleWithVelTag<scalarD> > ParticleSetVTD;

typedef ParticleSetTplWithTag<ParticleWithTag<scalar> > ParticleSetT;
typedef ParticleSetTplWithTag<ParticleWithVelTag<scalar> > ParticleSetVT;

typedef FrameSetTpl<ParticleSetB> FrameSetRawB;
typedef FrameSetTpl<ParticleSetBD> FrameSetRawBD;
typedef FrameSetTpl<ParticleSet> FrameSetRaw;
typedef FrameSetTpl<ParticleSetD> FrameSetRawD;
typedef FrameSetTpl<ParticleSetN> FrameSetRawN;
typedef FrameSetTpl<ParticleSetND> FrameSetRawND;
typedef FrameSetTpl<ParticleSetT> FrameSetT;
typedef FrameSetTpl<ParticleSetVT> FrameSetVT;

template <typename T>
struct HasPos<ParticleWithTag<T> > : public boost::true_type {};
template <typename T>
struct HasId<ParticleWithTag<T> > : public boost::true_type {};
template <typename T>
struct HasType<ParticleWithTag<T> > : public boost::true_type {};

template <typename T>
struct HasPos<ParticleWithVelTag<T> > : public boost::true_type {};
template <typename T>
struct HasVel<ParticleWithVelTag<T> > : public boost::true_type {};
template <typename T>
struct HasId<ParticleWithVelTag<T> > : public boost::true_type {};
template <typename T>
struct HasType<ParticleWithVelTag<T> > : public boost::true_type {};

PRJ_END

#endif