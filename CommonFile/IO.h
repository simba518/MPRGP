#ifndef IO_H
#define IO_H

#include "Config.h"
#include "MathBasic.h"
#include <Eigen/Sparse>

#include <iostream>
#include <map>

#include <boost/static_assert.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/interprocess/streams/vectorstream.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

using namespace std;

struct IOData;
struct Serializable
{
	Serializable(sizeType type):_type(type){}
	virtual bool read(istream& is){ASSERT_MSG(false,"Not Implemented!");return false;}
	virtual bool write(ostream& os) const{ASSERT_MSG(false,"Not Implemented!");return false;}
	virtual bool read(istream& is,IOData* dat){return read(is);}
	virtual bool write(ostream& os,IOData* dat) const{return write(os);}
	virtual boost::shared_ptr<Serializable> copy() const{ASSERT_MSG(false,"Not Implemented!");return boost::shared_ptr<Serializable>(new Serializable(_type));}
	sizeType type() const{return _type;}
private:
	sizeType _type;
};

//io for basic type
#define IO_BASIC(T)	\
FORCE_INLINE ostream& writeBinaryData(const T& val,ostream& os,IOData* dat=NULL){os.write((char*)&val,sizeof(T));return os;}	\
FORCE_INLINE istream& readBinaryData(T& val,istream& is,IOData* dat=NULL){is.read((char*)&val,sizeof(T));return is;}
IO_BASIC(char)
IO_BASIC(unsigned char)
IO_BASIC(short)
IO_BASIC(unsigned short)
IO_BASIC(int)
IO_BASIC(unsigned int)
IO_BASIC(scalarD)
IO_BASIC(bool)
IO_BASIC(sizeType)

//io for float is double
FORCE_INLINE ostream& writeBinaryData(const scalarF& val,ostream& os,IOData* dat=NULL)
{
    scalarD valD=val;
    os.write((char*)&valD,sizeof(scalarD));
    return os;
}
FORCE_INLINE istream& readBinaryData(scalarF& val,istream& is,IOData* dat=NULL)
{
    scalarD valD;
    is.read((char*)&valD,sizeof(scalarD));
    val=(scalarF)valD;
    return is;
}

//io for shared_ptr
struct IOData
{
public:
	typedef boost::shared_ptr<Serializable> SPTR;
	typedef boost::fast_pool_allocator<std::pair<const SPTR,sizeType> > ALLOC_MAP_WR;
	typedef boost::fast_pool_allocator<std::pair<const sizeType,SPTR> > ALLOC_MAP_RD;
	typedef boost::unordered_map<SPTR,sizeType,boost::hash<SPTR>,std::equal_to<SPTR>,ALLOC_MAP_WR> MAP_WR;
	typedef boost::unordered_map<sizeType,SPTR,boost::hash<sizeType>,std::equal_to<sizeType>,ALLOC_MAP_RD> MAP_RD;
	typedef boost::unordered_map<sizeType,SPTR,boost::hash<sizeType>,std::equal_to<sizeType>,ALLOC_MAP_RD> TYPE_SET;
public:
	IOData():_index(0){}
	sizeType getIndex(){return _index++;}
	void registerType(boost::shared_ptr<Serializable> type)
	{
		ASSERT_MSG(type->type() >= 0,"Given type doesn't support shared_ptr serialization!");
        ASSERT_MSGV(_tSet.find(type->type()) == _tSet.end(),"Conflicit type id: %ld",type->type());
		_tSet[type->type()]=type;
	}
	template <typename T>
	void registerType(){registerType(boost::shared_ptr<Serializable>(new T));}
	template <typename T>void createNew(istream& is,boost::shared_ptr<T>& val) const
	{
		sizeType type;
		readBinaryData(type,is);
		ASSERT_MSG(type != -1,"Type not found!")
		for(TYPE_SET::const_iterator beg=_tSet.begin(),end=_tSet.end();beg!=end;beg++)
			if(beg->first == type)
			{
                val=boost::dynamic_pointer_cast<T>(beg->second->copy());
				return;
			}
		ASSERT_MSG(false,"Cannot find compatible type!")
	}
	MAP_WR _ptrMapWr;
	MAP_RD _ptrMapRd;
private:
	sizeType _index;
	TYPE_SET _tSet;
};
template <typename T>FORCE_INLINE ostream& writeBinaryData(boost::shared_ptr<T> val,ostream& os,IOData* dat=NULL)
{
    ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
	if(val.get() == NULL)
	{
		writeBinaryData((sizeType)-1,os);
		return os;
	}
    boost::shared_ptr<Serializable> ptrS=boost::dynamic_pointer_cast<Serializable>(val);
	ASSERT_MSG(ptrS,"Not serializable type!")
	IOData::MAP_WR::const_iterator it=dat->_ptrMapWr.find(ptrS);
	if(it == dat->_ptrMapWr.end()){
		sizeType id=dat->getIndex();
		writeBinaryData(id,os,dat);
		writeBinaryData(val->type(),os,dat);
		dat->_ptrMapWr[val]=id;
		writeBinaryData(*val,os,dat);
	}else{
		writeBinaryData(it->second,os,dat);
	}
	return os;
}
template <typename T>FORCE_INLINE istream& readBinaryData(boost::shared_ptr<T>& val,istream& is,IOData* dat=NULL)
{
    ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
	sizeType id;
	readBinaryData(id,is,dat);
	if(id == -1)
	{
		val.reset((T*)NULL);
		return is;
	}
	IOData::MAP_RD::const_iterator it=dat->_ptrMapRd.find(id);
	if(it == dat->_ptrMapRd.end()){
		dat->createNew(is,val);
		dat->_ptrMapRd[id]=val;
		readBinaryData(*val,is,dat);
	}else{
        val=boost::dynamic_pointer_cast<T>(dat->_ptrMapRd[id]);
	}
    return is;
}
FORCE_INLINE ostream& writeBinaryData(const Serializable& val,ostream& os,IOData* dat=NULL){val.write(os,dat);return os;}
FORCE_INLINE istream& readBinaryData(Serializable& val,istream& is,IOData* dat=NULL){val.read(is,dat);return is;}

//io for fixed matrix
#define IO_FIXED(NAME,TO_TYPE)																	\
FORCE_INLINE ostream& writeBinaryData(const NAME& v,ostream& os,IOData* dat=NULL)				\
{																								\
	sizeType d0=NAME::RowsAtCompileTime;os.write((char*)&d0,sizeof(sizeType));					\
	sizeType d1=NAME::ColsAtCompileTime;os.write((char*)&d1,sizeof(sizeType));					\
	for(sizeType r=0;r<NAME::RowsAtCompileTime;r++)												\
	for(sizeType c=0;c<NAME::ColsAtCompileTime;c++)												\
	{																							\
		TO_TYPE val=(TO_TYPE)v(r,c);															\
		os.write((char*)&val,sizeof(TO_TYPE));													\
	}																							\
	return os;																					\
}																								\
FORCE_INLINE istream& readBinaryData(NAME& v,istream& is,IOData* dat=NULL)						\
{																								\
	sizeType d0;is.read((char*)&d0,sizeof(sizeType));ASSERT(d0 == NAME::RowsAtCompileTime)		\
	sizeType d1;is.read((char*)&d1,sizeof(sizeType));ASSERT(d1 == NAME::ColsAtCompileTime)		\
	for(sizeType r=0;r<NAME::RowsAtCompileTime;r++)												\
	for(sizeType c=0;c<NAME::ColsAtCompileTime;c++)												\
	{																							\
		TO_TYPE val;																			\
		is.read((char*)&val,sizeof(TO_TYPE));													\
		v(r,c)=(NAME::Scalar)val;																\
	}																							\
	return is;																					\
}

//io for non-fixed matrix
#define IO_NON_FIXED(NAME,TO_TYPE)																\
FORCE_INLINE ostream& writeBinaryData(const NAME& v,ostream& os,IOData* dat=NULL)				\
{																								\
	sizeType d0=v.rows();																		\
	sizeType d1=v.cols();																		\
	os.write((char*)&d0,sizeof(sizeType));														\
	os.write((char*)&d1,sizeof(sizeType));														\
    std::vector<TO_TYPE> data(d0*d1);                                                           \
    sizeType index=0;                                                                           \
	for(sizeType r=0;r<d0;r++)																	\
	for(sizeType c=0;c<d1;c++)																	\
		data[index++]=(TO_TYPE)v(r,c);															\
	if(d0*d1 > 0)																				\
    os.write((char*)&(data[0]),sizeof(TO_TYPE)*d0*d1);											\
	return os;																					\
}																								\
FORCE_INLINE istream& readBinaryData(NAME& v,istream& is,IOData* dat=NULL)						\
{																								\
	sizeType d0;is.read((char*)&d0,sizeof(sizeType));											\
	sizeType d1;is.read((char*)&d1,sizeof(sizeType));											\
	v.resize(d0,d1);																			\
    std::vector<TO_TYPE> data(d0*d1);                                                           \
	if(d0*d1 > 0)																				\
    is.read((char*)&(data[0]),sizeof(TO_TYPE)*d0*d1);											\
    sizeType index=0;                                                                           \
	for(sizeType r=0;r<d0;r++)																	\
	for(sizeType c=0;c<d1;c++)																	\
		v(r,c)=(NAME::Scalar)data[index++];														\
	return is;																					\
}

#define IO_FIXED_QUAT(NAME,TO_TYPE)																\
FORCE_INLINE ostream& writeBinaryData(const NAME& v,ostream& os,IOData* dat=NULL)				\
{																								\
	TO_TYPE val=(TO_TYPE)v.w();																	\
	os.write((char*)&val,sizeof(TO_TYPE));														\
	val=(TO_TYPE)v.x();																			\
	os.write((char*)&val,sizeof(TO_TYPE));														\
	val=(TO_TYPE)v.y();																			\
	os.write((char*)&val,sizeof(TO_TYPE));														\
	val=(TO_TYPE)v.z();																			\
	os.write((char*)&val,sizeof(TO_TYPE));														\
	return os;																					\
}																								\
FORCE_INLINE istream& readBinaryData(NAME& v,istream& is,IOData* dat=NULL)						\
{																								\
	TO_TYPE val;																				\
	is.read((char*)&val,sizeof(TO_TYPE));														\
	v.w()=(NAME::Scalar)val;																	\
	is.read((char*)&val,sizeof(TO_TYPE));														\
	v.x()=(NAME::Scalar)val;																	\
	is.read((char*)&val,sizeof(TO_TYPE));														\
	v.y()=(NAME::Scalar)val;																	\
	is.read((char*)&val,sizeof(TO_TYPE));														\
	v.z()=(NAME::Scalar)val;																	\
	return is;																					\
}

IO_FIXED(Vec2d,scalarD)
IO_FIXED(Vec3d,scalarD)
IO_FIXED(Vec4d,scalarD)
IO_FIXED(Mat2d,scalarD)
IO_FIXED(Mat3d,scalarD)
IO_FIXED(Mat4d,scalarD)
IO_FIXED_QUAT(Quatd,scalarD)

IO_FIXED(Vec2f,scalarD)
IO_FIXED(Vec3f,scalarD)
IO_FIXED(Vec4f,scalarD)
IO_FIXED(Mat2f,scalarD)
IO_FIXED(Mat3f,scalarD)
IO_FIXED(Mat4f,scalarD)
IO_FIXED_QUAT(Quatf,scalarD)

IO_FIXED(Vec2i,sizeType)
IO_FIXED(Vec3i,sizeType)
IO_FIXED(Vec4i,sizeType)

IO_FIXED(Vec2c,char)
IO_FIXED(Vec3c,char)
IO_FIXED(Vec4c,char)

IO_NON_FIXED(Rowd,scalarD)
IO_NON_FIXED(Cold,scalarD)
IO_NON_FIXED(Matd,scalarD)

IO_NON_FIXED(Rowf,scalarD)
IO_NON_FIXED(Colf,scalarD)
IO_NON_FIXED(Matf,scalarD)

IO_NON_FIXED(Rowi,sizeType)
IO_NON_FIXED(Coli,sizeType)
IO_NON_FIXED(Mati,sizeType)

template <typename T1,typename T2>
FORCE_INLINE ostream& writeBinaryData(const std::pair<T1,T2>& m,ostream& os,IOData* dat=NULL)
{
	writeBinaryData(m.first,os,dat);
	writeBinaryData(m.second,os,dat);
	return os;
}
template <typename T1,typename T2>
FORCE_INLINE istream& readBinaryData(std::pair<T1,T2>& m,istream& is,IOData* dat=NULL)
{
	readBinaryData(m.first,is,dat);
	readBinaryData(m.second,is,dat);
	return is;
}

#define IO_BB(NAME)																				\
FORCE_INLINE ostream& writeBinaryData(const BBox<NAME>& b,ostream& os,IOData* dat=NULL)			\
{																								\
	writeBinaryData(b._minC,os,dat);															\
	writeBinaryData(b._maxC,os,dat);															\
	return os;																					\
}																								\
FORCE_INLINE istream& readBinaryData(BBox<NAME>& b,istream& is,IOData* dat=NULL)				\
{																								\
	readBinaryData(b._minC,is,dat);																\
	readBinaryData(b._maxC,is,dat);																\
	return is;																					\
}

IO_BB(scalarD)
IO_BB(scalarF)

FORCE_INLINE istream& readBinaryData(string& str,istream& is,IOData* dat=NULL)
{
    sizeType len;
    readBinaryData(len,is,dat);
    str.assign(len,' ');
    is.read(&(str[0]),len);
    return is;
}
FORCE_INLINE ostream& writeBinaryData(const string& str,ostream& os,IOData* dat=NULL)
{
    const sizeType len=(sizeType)str.length();
    writeBinaryData(len,os,dat);
    return os.write(str.c_str(),len);
}

template <typename T,typename ALLOC>FORCE_INLINE istream& readVector(vector<T,ALLOC>& v,istream& is,IOData* dat=NULL)
{
    sizeType size;
    readBinaryData(size,is,dat);
    v.resize(size);
    for(sizeType i=0; i<(sizeType)v.size(); i++)
        readBinaryData(v[i],is,dat);
    return is;
}
template <typename T,typename ALLOC>FORCE_INLINE ostream& writeVector(const vector<T,ALLOC>& v,ostream& os,IOData* dat=NULL)
{
    sizeType size=(sizeType)v.size();
    writeBinaryData(size,os,dat);
    for(sizeType i=0; i<(sizeType)v.size(); i++)
        writeBinaryData(v[i],os,dat);
    return os;
}

FORCE_INLINE istream& readVector(vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,istream& is,IOData* dat=NULL)
{
    sizeType size;
    readBinaryData(size,is,dat);
    v.resize(size);
    if(size > 0)
        is.read((char*)&(v[0]),sizeof(scalarD)*size);
    return is;
}
FORCE_INLINE ostream& writeVector(const vector<scalarD,Eigen::aligned_allocator<scalarD> >& v,ostream& os,IOData* dat=NULL)
{
    sizeType size=(sizeType)v.size();
    writeBinaryData(size,os,dat);
    if(size > 0)
        os.write((char*)&(v[0]),sizeof(scalarD)*size);
    return os;
}

FORCE_INLINE istream& readVector(vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,istream& is,IOData* dat=NULL)
{
    std::vector<scalarD,Eigen::aligned_allocator<scalarD> > tmp;
    readVector(tmp,is,dat);

    v.resize(tmp.size());
    for(sizeType i=0; i<(sizeType)v.size(); i++)
        v[i]=(scalarF)tmp[i];
    return is;
}
FORCE_INLINE ostream& writeVector(const vector<scalarF,Eigen::aligned_allocator<scalarF> >& v,ostream& os,IOData* dat=NULL)
{
    std::vector<scalarD,Eigen::aligned_allocator<scalarD> > tmp(v.size());
    for(sizeType i=0; i<(sizeType)v.size(); i++)
        tmp[i]=v[i];

    writeVector(tmp,os,dat);
    return os;
}

#define IO_VEC_FIXED(NAME)																									\
FORCE_INLINE istream& readVector(vector<NAME,Eigen::aligned_allocator<NAME> >& v,istream& is,IOData* dat=NULL)				\
{																															\
	sizeType size;																											\
	readBinaryData(size,is,dat);																							\
	v.resize(size);																											\
	for(sizeType i=0;i<(sizeType)v.size();i++)																				\
		readBinaryData(v[i],is,dat);																						\
	return is;																												\
}																															\
FORCE_INLINE ostream& writeVector(const vector<NAME,Eigen::aligned_allocator<NAME> >& v,ostream& os,IOData* dat=NULL)		\
{																															\
	sizeType size=(sizeType)v.size();																						\
	writeBinaryData(size,os,dat);																							\
	for(sizeType i=0;i<(sizeType)v.size();i++)																				\
		writeBinaryData(v[i],os,dat);																						\
	return os;																												\
}

IO_VEC_FIXED(Vec2)
IO_VEC_FIXED(Vec3)
IO_VEC_FIXED(Vec4)
IO_VEC_FIXED(Mat2)
IO_VEC_FIXED(Mat3)
IO_VEC_FIXED(Mat4)
IO_VEC_FIXED(Vec2i)
IO_VEC_FIXED(Vec3i)
IO_VEC_FIXED(Vec4i)

class HasMagic
{
public:
    HasMagic(const sizeType& MAGIC):_MAGIC(MAGIC) {}
    virtual ~HasMagic() {}
    bool readMagic(istream& is) {
        if(!is.good())
            return false;

        sizeType MAGIC_Test;
        readBinaryData(MAGIC_Test,is);
        if(MAGIC_Test != _MAGIC)
            return false;

        return true;
    }
    bool writeMagic(ostream& os) const {
        if(!os.good())
            return false;

        writeBinaryData(_MAGIC,os);
        return true;
    }
    HasMagic& operator=(const HasMagic& other) {
        ASSERT(other._MAGIC == _MAGIC)return *this;
    }
    const sizeType _MAGIC;
};

FORCE_INLINE bool beginsWith(const string& name,const string& regex){return name.size() >= regex.size() && name.substr(0,regex.size()) == regex;}
FORCE_INLINE bool endsWith(const string& name,const string& regex){return name.size() >= regex.size() && name.substr(name.size()-regex.size(),regex.size()) == regex;}
FORCE_INLINE string replaceDot(const string& name)
{
    string ret=name;
    for(sizeType i=0,nr=(sizeType)ret.size(); i<nr; i++)
        if(ret[i]=='.')
            ret[i]='_';
    return ret;
}

class Endianness
{
public:
    static bool isLittleEndian() {
        union u {
            unsigned long l;
            unsigned char c[sizeof(unsigned long)];
        };
        u dummy;
        dummy.l = 1;
        return dummy.c[0] == 1;
    }
    static void swap2Bytes(unsigned char* &ptr) {
        unsigned char tmp;
        tmp = ptr[0];
        ptr[0] = ptr[1];
        ptr[1] = tmp;
    }
    static void swap4Bytes(unsigned char* &ptr) {
        unsigned char tmp;
        tmp = ptr[0];
        ptr[0] = ptr[3];
        ptr[3] = tmp;
        tmp = ptr[1];
        ptr[1] = ptr[2];
        ptr[2] = tmp;
    }
    static void swap8Bytes(unsigned char* &ptr) {
        unsigned char tmp;
        tmp = ptr[0];
        ptr[0] = ptr[7];
        ptr[7] = tmp;
        tmp = ptr[1];
        ptr[1] = ptr[6];
        ptr[6] = tmp;
        tmp = ptr[2];
        ptr[2] = ptr[5];
        ptr[5] = tmp;
        tmp = ptr[3];
        ptr[3] = ptr[4];
        ptr[4] = tmp;
    }
};
template<class T2>
FORCE_INLINE void byteSwap(T2& val)
{
    sizeType n=sizeof(T2);
    unsigned char *p=(unsigned char*)&val;
    switch( n ) {
    case 1:
        return;
    case 2:
        Endianness::swap2Bytes(p);
        break;
    case 4:
        Endianness::swap4Bytes(p);
        break;
    case 8:
        Endianness::swap8Bytes(p);
        break;
    default:
        break;
    }
}
template <typename T2>
FORCE_INLINE void vtkWrite(ostream& oss,T2 val)
{
    if(Endianness::isLittleEndian())
        byteSwap(val);
    oss.write((const char*)&val,sizeof(T2));
}
template<typename T>
struct VTKWriter {
    enum VTK_DATA_TYPE {
        UNSTRUCTURED_GRID,
        STRUCTURED_POINTS,
    };
    enum VTK_CELL_TYPE {
        POINT=1,
        LINE=3,
        TRIANGLE=5,
        TETRA=10,
        VOXEL=11,
        HEX=12,
        POLYLINE=4,
        QUAD=9,
		QUADRATIC_LINE=21,
    };
public:
	struct Data
	{
		Data():_nr(0){}
		string _str;
		sizeType _nr;
	};
    template <typename ITER>
    struct ValueTraits {
        typedef typename ITER::value_type value_type;
    };
    template <typename POINTED_TYPE>
    struct ValueTraits<POINTED_TYPE*> {
        typedef POINTED_TYPE value_type;
    };
    template <typename ID>
    struct IteratorIndex {
        typedef ID value_type;
        IteratorIndex(const sizeType& id,const sizeType stride,const sizeType& off)
            :_id(id),_stride(stride),_off(off) {}
        void operator++() {
            _id++;
        }
        bool operator!=(const IteratorIndex& other) const {
            return _id != other._id;
        }
        virtual ID operator*() const {
            ID ret;
            for(sizeType i=0; i<ret.size(); i++)ret(i)=(_stride == 0) ? _id+_off*i : _id*_stride+i;
            return ret;
        }
        sizeType _id;
        sizeType _stride;
        sizeType _off;
    };
    template <typename ITER>
    struct IteratorAdd {
        typedef typename ValueTraits<ITER>::value_type value_type;
        IteratorAdd(ITER beg0,ITER beg1):_beg0(beg0),_beg1(beg1) {}
        void operator++() {
            _beg0++;
            _beg1++;
        }
        bool operator!=(const IteratorAdd& other) const {
            return _beg0 != other._beg0;
        }
        virtual value_type operator*() const {
            return (*_beg0)+(*_beg1);
        }
        ITER _beg0,_beg1;
    };
    template <typename ITER,typename SCALAR>
    struct IteratorAddMult {
        typedef typename ValueTraits<ITER>::value_type value_type;
        IteratorAddMult(ITER beg0,ITER beg1,SCALAR mult):_beg0(beg0),_beg1(beg1),_mult(mult) {}
        void operator++() {
            _beg0++;
            _beg1++;
        }
        bool operator!=(const IteratorAddMult& other) const {
            return _beg0 != other._beg0;
        }
        virtual value_type operator*() const {
            return (*_beg0)+(*_beg1)*_mult;
        }
        ITER _beg0,_beg1;
        SCALAR _mult;
    };
public:
    VTKWriter(const string& name,const boost::filesystem::path& path,bool binary)
        :_os(path,binary ? ios_base::binary : ios_base::out),
         _points(binary ? ios_base::binary : ios_base::out),
         _cells(binary ? ios_base::binary : ios_base::out),
         _cellTypes(binary ? ios_base::binary : ios_base::out),
         _cellDatas(binary ? ios_base::binary : ios_base::out),
         _nrPoint(0),_nrCell(0),_nrIndex(0),_nrData(0),_vdt(UNSTRUCTURED_GRID),
         _binary(binary) 
    {
        _os << "# vtk DataFile Version 1.0" << endl;
        _os << name << endl;
        _os << (binary ? "BINARY" : "ASCII") << endl;
        _os << "DATASET " << "UNSTRUCTURED_GRID" << endl;
    }
    VTKWriter(const string& name,const boost::filesystem::path& path,bool binary,const BBox<T>& bb,const Vec3i& nrCell,bool center)
        :_os(path,binary ? ios_base::binary : ios_base::out),
         _points(binary ? ios_base::binary : ios_base::out),
         _cells(binary ? ios_base::binary : ios_base::out),
         _cellTypes(binary ? ios_base::binary : ios_base::out),
         _cellDatas(binary ? ios_base::binary : ios_base::out),
         _nrPoint(0),_nrCell(0),_nrIndex(0),_nrData(0),_vdt(STRUCTURED_POINTS),
         _binary(binary) 
    {
        typename BBox<T>::PT ext=bb.getExtent();
        typename BBox<T>::PT spacing(ext.x()/nrCell.x(),ext.y()/nrCell.y(),ext.z()/nrCell.z());
        _os << "# vtk DataFile Version 1.0" << endl;
        _os << name << endl;
        _os << (binary ? "BINARY" : "ASCII") << endl;
        _os << "DATASET " << "STRUCTURED_POINTS" << endl;
        if(center) {
            typename BBox<T>::PT origin=bb._minC+spacing*0.5f;
            _os << "DIMENSIONS " << nrCell.x() << " " << nrCell.y() << " " << nrCell.z() << endl;
            _os << "ORIGIN " << origin.x() << " " << origin.y() << " " << origin.z() << endl;
            _os << "SPACING " << spacing.x() << " " << spacing.y() << " " << spacing.z() << endl;
        } else {
            typename BBox<T>::PT origin=bb._minC;
            _os << "DIMENSIONS " << (nrCell.x()+1) << " " << (nrCell.y()+1) << " " << (nrCell.z()+1) << endl;
            _os << "ORIGIN " << origin.x() << " " << origin.y() << " " << origin.z() << endl;
            _os << "SPACING " << spacing.x() << " " << spacing.y() << " " << spacing.z() << endl;
        }
    }
    virtual ~VTKWriter() {
		bool first;
        switch(_vdt) {
        case UNSTRUCTURED_GRID:
            _os << "POINTS " << _nrPoint << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
            _os << _points.str();
            _os << "CELLS " << _nrCell << " " << _nrIndex << endl;
            _os << _cells.str();
            _os << "CELL_TYPES " << _nrCell << endl;
            _os << _cellTypes.str();
			first=false;
			for(typename std::map<string,Data>::const_iterator beg=_customData.begin(),end=_customData.end();beg!=end;beg++)
            {
				if(!first)
					_os << "CELL_DATA " << beg->second._nr << endl;
				first=true;
				_os << beg->second._str << endl;
			}
			first=false;
            for(typename std::map<string,Data>::const_iterator beg=_customPointData.begin(),end=_customPointData.end();beg!=end;beg++)
            {
				if(!first)
					_os << "POINT_DATA " << beg->second._nr << endl;
				first=true;
                _os << beg->second._str << endl;
			}
            break;
        case STRUCTURED_POINTS:
            ASSERT(_nrData == 1)
            _os << "POINT_DATA " << _nrPoint << std::endl;
            _os << "SCALARS data " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
            _os << "LOOKUP_TABLE default" << endl;
            _os << _cellDatas.str();
            break;
        default:
            ASSERT_MSG(false,"Unsupported!")
        }
    }
    template <typename ITER> VTKWriter& appendPoints(ITER beg,ITER end) {
        typedef typename ValueTraits<ITER>::value_type value_type;
        sizeType nr=0;
        if(_binary) {
            for(; beg != end; ++beg) {
                const value_type& val=*beg;
                vtkWrite<T>(_points,val(0));
                vtkWrite<T>(_points,val(1));
                vtkWrite<T>(_points,val(2));
                nr++;
            }
        } else {
            for(; beg != end; ++beg) {
                const value_type& val=*beg;
                _points << (T)val(0) << " " << (T)val(1) << " " << (T)val(2) << endl;
                nr++;
            }
        }
        _nrPoint+=nr;
        return *this;
    }
    template <typename ITER> VTKWriter& appendVoxels(ITER beg,ITER end,bool hex) {
        typedef typename Eigen::Matrix<sizeType,8,1> IDS;
        typedef typename ValueTraits<ITER>::value_type value_type;
        std::vector<value_type,Eigen::aligned_allocator<value_type> > points;
        std::vector<IDS,Eigen::aligned_allocator<IDS> > cells;
        for(; beg!=end;) {
            IDS ids;
            sizeType base=(sizeType)points.size();

            value_type minC=*beg++;
            value_type maxC=*beg++;
            value_type ext=maxC-minC;

            points.push_back(minC+value_type(0.0f   ,0.0f   ,0.0f   ));
            points.push_back(minC+value_type(ext.x(),   0.0f,0.0f   ));
            points.push_back(minC+value_type(0.0f   ,ext.y(),0.0f   ));
            points.push_back(minC+value_type(ext.x(),ext.y(),0.0f   ));

            points.push_back(minC+value_type(0.0f   ,0.0f   ,ext.z()));
            points.push_back(minC+value_type(ext.x(),   0.0f,ext.z()));
            points.push_back(minC+value_type(0.0f   ,ext.y(),ext.z()));
            points.push_back(minC+value_type(ext.x(),ext.y(),ext.z()));

            if(hex)
                ids << base  ,base+1,base+3,base+2,
                       base+4,base+5,base+7,base+6;
            else
                ids << base  ,base+1,base+2,base+3,
                       base+4,base+5,base+6,base+7;
            cells.push_back(ids);
        }
        appendPoints(points.begin(),points.end());
        appendCells(cells.begin(),cells.end(),hex ? HEX : VOXEL);
        return *this;
    }
    template <typename ITER> VTKWriter& appendCells(ITER beg,ITER end,VTK_CELL_TYPE ct,bool relativeIndex=false) {
        if(relativeIndex)
            ASSERT(_relativeCellIndex >= -1)
        int base=relativeIndex ? (int)_relativeCellIndex : 0;
        
        typedef typename ValueTraits<ITER>::value_type value_type;
        sizeType nr=0;
        sizeType nrIndex=0;
        if(_binary) {
            for(; beg != end; ++beg) {
                const value_type& val=*beg;
                switch(ct) {
                case POINT:
                    nrIndex+=2;
                    vtkWrite<int>(_cells,1);
                    vtkWrite<int>(_cells,base+(int)val(0));
                    break;
                case LINE:
                    nrIndex+=3;
                    vtkWrite<int>(_cells,2);
                    vtkWrite<int>(_cells,base+(int)val(0));
                    vtkWrite<int>(_cells,base+(int)val(1));
                    break;
                case TRIANGLE:
				case QUADRATIC_LINE:
                    nrIndex+=4;
                    vtkWrite<int>(_cells,3);
                    vtkWrite<int>(_cells,base+(int)val(0));
                    vtkWrite<int>(_cells,base+(int)val(1));
                    vtkWrite<int>(_cells,base+(int)val(2));
                    break;
                case TETRA:
                case QUAD:
                    nrIndex+=5;
                    vtkWrite<int>(_cells,4);
                    vtkWrite<int>(_cells,base+(int)val(0));
                    vtkWrite<int>(_cells,base+(int)val(1));
                    vtkWrite<int>(_cells,base+(int)val(2));
                    vtkWrite<int>(_cells,base+(int)val(3));
                    break;
                case VOXEL:
                case HEX:
                    nrIndex+=9;
                    vtkWrite<int>(_cells,8);
                    vtkWrite<int>(_cells,base+(int)val(0));
                    vtkWrite<int>(_cells,base+(int)val(1));
                    vtkWrite<int>(_cells,base+(int)val(2));
                    vtkWrite<int>(_cells,base+(int)val(3));
                    vtkWrite<int>(_cells,base+(int)val(4));
                    vtkWrite<int>(_cells,base+(int)val(5));
                    vtkWrite<int>(_cells,base+(int)val(6));
                    vtkWrite<int>(_cells,base+(int)val(7));
                    break;
                case POLYLINE:
                    nrIndex+=val.rows()+1;
                    vtkWrite<int>(_cells,(int)val.rows());
                    for(sizeType i=0;i<(int)val.rows();i++)
                        vtkWrite<int>(_cells,base+(int)val[i]);
                    break;
                }
                vtkWrite<int>(_cellTypes,(int)ct);
                nr++;
            }
        } else {
            for(; beg != end; ++beg) {
                const value_type& val=*beg;
                switch(ct) {
                case POINT:
                    nrIndex+=2;
                    _cells << "1 " << (base+(int)val(0)) << endl;
                    break;
                case LINE:
                    nrIndex+=3;
                    _cells << "2 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << endl;
                    break;
                case TRIANGLE:
				case QUADRATIC_LINE:
                    nrIndex+=4;
                    _cells << "3 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << endl;
                    break;
                case TETRA:
                case QUAD:
                    nrIndex+=5;
                    _cells << "4 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << " " << (base+(int)val(3)) << endl;
                    break;
                case VOXEL:
                case HEX:
                    nrIndex+=9;
                    _cells << "8 " << (base+(int)val(0)) << " " << (base+(int)val(1)) << " " << (base+(int)val(2)) << " " << (base+(int)val(3)) << " "
                                   << (base+(int)val(4)) << " " << (base+(int)val(5)) << " " << (base+(int)val(6)) << " " << (base+(int)val(7)) << endl;
                    break;
                case POLYLINE:
                    nrIndex+=val.rows()+1;
                    _cells << val.rows() << " ";
                    for(sizeType i=0;i<(int)val.rows();i++)
                        _cells << (base+(int)val[i]) << " ";
                    _cells << endl;
                    break;
                }
                _cellTypes << ct << endl;
                nr++;
            }
        }
        _nrCell+=nr;
        _nrIndex+=nrIndex;
        return *this;
    }
    template <typename ITER> VTKWriter& appendDatas(const string name,ITER beg,ITER end) {
        if(_binary)
            for(; beg != end; ++beg,++_nrPoint)
                vtkWrite<T>(_cellDatas,*beg);
        else
            for(; beg != end; ++beg,++_nrPoint)
                _cellDatas << (T)*beg << endl;
        _nrData++;
        return *this;
    }
    template <typename ITER> VTKWriter& appendCustomData(const string name,ITER beg,ITER end) 
	{    
		ostringstream os;
		if(_customData.find(name) == _customData.end())
		{
			os << "SCALARS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
			os << "LOOKUP_TABLE default" << endl;
		}
		
		Data& dat=_customData[name];
        if(_binary)
            for(; beg != end; ++beg,++dat._nr)
                vtkWrite<T>(os,(T)*beg);
        else
            for(; beg != end; ++beg,++dat._nr)
                os << (T)*beg << endl;
		dat._str+=os.str();
        ASSERT(dat._nr == _nrCell)
        return *this;
    }
    template <typename ITER> VTKWriter& appendCustomPointData(const string name,ITER beg,ITER end) 
	{
		ostringstream os;
        if(_customPointData.find(name) == _customPointData.end())
		{
			Data& dat=_customPointData[name];
			os << "SCALARS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
			os << "LOOKUP_TABLE default" << endl;
		}
		
		Data& dat=_customPointData[name];
        if(_binary)
            for(; beg != end; ++beg,++dat._nr)
                vtkWrite<T>(os,(T)*beg);
        else
            for(; beg != end; ++beg,++dat._nr)
                os << (T)*beg << endl;
		dat._str+=os.str();
        ASSERT(dat._nr == _nrPoint)
        return *this;
    }
	template <typename ITER> VTKWriter& appendCustomVectorData(const string name,ITER beg,ITER end)
	{
		ostringstream os;
		if(_customData.find(name) == _customData.end())
		{
			os << "VECTORS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
		}
		
		Data& dat=_customData[name];
        if(_binary)
            for(; beg != end; ++beg,++dat._nr){
                vtkWrite<T>(os,(T)(*beg)[0]);
                vtkWrite<T>(os,(T)(*beg)[1]);
                vtkWrite<T>(os,(T)(*beg)[2]);
			}
        else
            for(; beg != end; ++beg,++dat._nr)
                os << (T)(*beg)[0] << " " << (T)(*beg)[1] << " " << (T)(*beg)[2] << endl;
		dat._str+=os.str();
        ASSERT(dat._nr == _nrCell)
        return *this;
	}
	template <typename ITER> VTKWriter& appendCustomPointVectorData(const string name,ITER beg,ITER end)
	{
		ostringstream os;
        if(_customPointData.find(name) == _customPointData.end())
		{
			Data& dat=_customPointData[name];
			os << "VECTORS " << name << " " << (sizeof(T) == sizeof(float) ? "float" : "double") << endl;
		}
		
		Data& dat=_customPointData[name];
        if(_binary)
            for(; beg != end; ++beg,++dat._nr){
                vtkWrite<T>(os,(T)(*beg)[0]);
                vtkWrite<T>(os,(T)(*beg)[1]);
                vtkWrite<T>(os,(T)(*beg)[2]);
			}
        else
            for(; beg != end; ++beg,++dat._nr)
                os << (T)(*beg)[0] << " " << (T)(*beg)[1] << " " << (T)(*beg)[2] << endl;
		dat._str+=os.str();
        ASSERT(dat._nr == _nrPoint)
        return *this;
	}
    template <typename ITER> VTKWriter& appendPointsByAdd(ITER beg0,ITER beg1,ITER end0) {
        appendPoints(IteratorAdd<ITER>(beg0,beg1),IteratorAdd<ITER>(end0,end0));
        return *this;
    }
    void setRelativeIndex(sizeType rel=-1)
    {
        if(rel == -1)
            _relativeCellIndex=_nrPoint;
        else _relativeCellIndex=rel;
    }
private:
    boost::filesystem::ofstream _os;
    ostringstream _points,_cells,_cellTypes,_cellDatas;
    std::map<string,Data> _customData;
    std::map<string,Data> _customPointData;
    sizeType _nrPoint,_nrCell,_nrIndex,_nrData;
    VTK_DATA_TYPE _vdt;
    bool _binary;
    sizeType _relativeCellIndex;
};

PRJ_END

#endif
