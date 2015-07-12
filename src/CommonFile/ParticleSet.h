#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include "Config.h"
#include "MathBasic.h"
#include "IO.h"
#include <boost/type_traits/integral_constant.hpp>

PRJ_BEGIN

template <typename T>
struct HasPos : boost::false_type {};
template <typename T>
struct HasVel : boost::false_type {};
template <typename T>
struct HasNormal : boost::false_type {};

template <typename A,typename B,bool has>
struct CopyPos {
    static FORCE_INLINE void copy(A& a,const B& b) {
        a._pos.x()=(typename A::scalarType)b._pos.x();
        a._pos.y()=(typename A::scalarType)b._pos.y();
        a._pos.z()=(typename A::scalarType)b._pos.z();
    }
};
template <typename A,typename B>
struct CopyPos<A,B,false> {
    static FORCE_INLINE void copy(A& a,const B& b) {}
};

template <typename A,typename B,bool has>
struct CopyVel {
    static FORCE_INLINE void copy(A& a,const B& b) {
        a._vel.x()=(typename A::scalarType)b._vel.x();
        a._vel.y()=(typename A::scalarType)b._vel.y();
        a._vel.z()=(typename A::scalarType)b._vel.z();
    }
};
template <typename A,typename B>
struct CopyVel<A,B,false> {
    static FORCE_INLINE void copy(A& a,const B& b) {}
};

template <typename A,typename B,bool has>
struct CopyNormal {
    static FORCE_INLINE void copy(A& a,const B& b) {
        a._normal.x()=(typename A::scalarType)b._normal.x();
        a._normal.y()=(typename A::scalarType)b._normal.y();
        a._normal.z()=(typename A::scalarType)b._normal.z();
    }
};
template <typename A,typename B>
struct CopyNormal<A,B,false> {
    static FORCE_INLINE void copy(A& a,const B& b) {}
};

//============================================================================
// Flip is one of the most efficient surface tracker
template <typename P_TYPE>
struct ParticleSetTpl : public HasMagic {
    typedef typename P_TYPE::scalarType scalarType;
    typedef typename ScalarUtil<scalarType>::ScalarVec3 Vec3Type;
    typedef P_TYPE ParticleType;
    typedef P_TYPE value_type;
    friend class DataStructureCL;
public:
    ParticleSetTpl():HasMagic(0xFFFFFFFFCCCCCCCC) {}
    virtual ~ParticleSetTpl() {}
    virtual bool read(istream& is) {
        if(!HasMagic::readMagic(is))
            return false;
        readVector(_pSet,is);
        return is.good();
    }
    virtual bool write(ostream& os) const {
        if(!HasMagic::writeMagic(os))
            return false;
		writeVector(_pSet,os);
        return os.good();
    }
    bool readFrame(istream& is) {
        sizeType nr;
        char comma;
        is >> nr;
        resize(nr);
        for(sizeType i=0; i<nr; i++) {
            Vec3Type& pos=get(i)._pos;
            is >> comma >> pos.x();
            ASSERT(comma == ',');
            is >> comma >> pos.y();
            ASSERT(comma == ',');
            is >> comma >> pos.z();
            ASSERT(comma == ',');
        }
        /*for(sizeType i=0;i<nr;i++)
        {
        	vec3Type& vel=get(i)._vel;
        	is >> comma >> vel.x();ASSERT(comma == ',');
        	is >> comma >> vel.y();ASSERT(comma == ',');
        	is >> comma >> vel.z();ASSERT(comma == ',');
        }*/
        return is.good();
    }
    bool writeFrame(ostream& os) const {
        os << size() << ",";
        for(sizeType i=0; i<(sizeType)_pSet.size(); i++) {
            const Vec3Type& pos=get(i)._pos;
            os << pos.x() << "," << pos.y() << "," << pos.z() << ",";
        }
        /*for(sizeType i=0;i<(sizeType)_pSet.size();i++)
        {
        	const vec3Type& vel=get(i)._vel;
        	os << vel.x() << "," << vel.y() << "," << vel.z();
        	if(i < (sizeType)_pSet.size()-1)
        		os << ",";
        }*/
        return os.good();
    }
	void writeVTK(const std::string& path) const {
        VTKWriter<scalarType> writer("Particles",boost::filesystem::path(path),true);
		writeVTK(writer);
	}
    void writeVTK(VTKWriter<scalarType>& writer) const {
        std::vector<Vec3Type,Eigen::aligned_allocator<Vec3Type> > pos;
        const sizeType n=size();
        for(sizeType i=0; i<n; i++)
            pos.push_back(_pSet[i]._pos);

        writer.setRelativeIndex();
		writer.appendPoints(pos.begin(),pos.end());
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,1,0);
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),1,0);
        writer.appendCells(itBeg,itEnd,VTKWriter<scalarType>::POINT,true);
    }
    void writeNormalVTK(const string& path,scalar len) const
    {
        std::vector<Vec3Type,Eigen::aligned_allocator<Vec3Type> > pos;
        const sizeType n=size();
        for(sizeType i=0; i<n; i++)
        {
            pos.push_back(_pSet[i]._pos);
            pos.push_back(_pSet[i]._pos+_pSet[i]._normal*len);
        }

        VTKWriter<scalarType> writer("Particle Normals",boost::filesystem::path(path),true);
        writer.appendPoints(pos.begin(),pos.end());
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,2,0);
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),2,0);
        writer.appendCells(itBeg,itEnd,
                           VTKWriter<scalarType>::LINE);
    }
    void writeVelVTK(const string& path,const scalarType len) const
    {
        std::vector<Vec3Type,Eigen::aligned_allocator<Vec3Type> > pos;
		vector<scalarType> css;
        const sizeType n=size();
        for(sizeType i=0; i<n; i++)
        {
            pos.push_back(_pSet[i]._pos);
            pos.push_back(_pSet[i]._pos+_pSet[i]._vel*len);
			css.push_back(0.0f);
			css.push_back(1.0f);
        }

        VTKWriter<scalarType> writer("Particle Normals",boost::filesystem::path(path),true);
        writer.appendPoints(pos.begin(),pos.end());
		writer.appendCustomPointData("Color",css.begin(),css.end());
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,2,0);
        typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),2,0);
        writer.appendCells(itBeg,itEnd,VTKWriter<scalarType>::LINE);
    }
    sizeType size() const {
        return _pSet.size();
    }
    FORCE_INLINE P_TYPE& operator[](const sizeType& i) {
        return get(i);
    }
    FORCE_INLINE const P_TYPE& operator[](const sizeType& i) const {
        return get(i);
    }
    FORCE_INLINE P_TYPE& get(const sizeType& i) {
        return _pSet[i];
    }
    FORCE_INLINE const P_TYPE& get(const sizeType& i) const {
        return _pSet[i];
    }
    //P_TYPE* getPtr(){return &(_pSet[0]);}
    //const P_TYPE* getPtr() const{return &(_pSet[0]);}
    void clear() {
        _pSet.clear();
    }
    void resize(const sizeType& n) {
        _pSet.resize(n);
    }
    void addParticle(const P_TYPE& p) {
        _pSet.push_back(p);
    }
    void swap(ParticleSetTpl& other) {
        _pSet.swap(other._pSet);
    }
    void append(const ParticleSetTpl& other) {
        _pSet.insert(_pSet.end(),other.begin(),other.end());
    }
    typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::iterator begin() {
        return _pSet.begin();
    }
    typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::iterator end() {
        return _pSet.end();
    }
    typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::const_iterator begin() const {
        return _pSet.begin();
    }
    typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::const_iterator end() const {
        return _pSet.end();
    }
    template<typename P_TYPE_FROM>
    ParticleSetTpl& copy(const ParticleSetTpl<P_TYPE_FROM>& other) {
        resize(other.size());
        for(sizeType i=0; i<other.size(); i++)
            get(i).copy<P_TYPE_FROM>(other.get(i));
        return *this;
    }
    template<typename COMP>
    void sort(const COMP& cmp) {
        std::sort(_pSet.begin(),_pSet.end(),cmp);
    }
protected:
    std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> > _pSet;
};
template <typename T>
struct ParticleBase : public Serializable {
    typedef T scalarType;
    ALIGN_16 typename ScalarUtil<T>::ScalarVec3 _pos;
    ParticleBase():_pos(0.0f,0.0f,0.0f),Serializable(-1) {}
    template<typename OTHER>
    ParticleBase& copy(const OTHER& p) {
        CopyPos<ParticleBase,OTHER,HasPos<OTHER>::value>::copy(*this,p);
        return *this;
    }
    virtual bool read(istream& is) {
        return readBinaryData(_pos,is).good();
    }
    virtual bool write(ostream& os) const {
        return writeBinaryData(_pos,os).good();
    }
    virtual bool operator==(const ParticleBase& other) const{
        return other._pos == _pos;
    }
};
template <typename T>
struct Particle : public ParticleBase<T> {
    ALIGN_16 typename ScalarUtil<T>::ScalarVec3 _vel;
    Particle():_vel(0.0f,0.0f,0.0f) {}
    template<typename OTHER>
    Particle& copy(const OTHER& p) {
        ParticleBase<T>::copy(p);
        CopyVel<Particle,OTHER,HasVel<OTHER>::value>::copy(*this,p);
        return *this;
    }
    virtual bool read(istream& is) {
        return ParticleBase<T>::read(is) && readBinaryData(_vel,is).good();
    }
    virtual bool write(ostream& os) const {
        return ParticleBase<T>::write(os) && writeBinaryData(_vel,os).good();
    }
    virtual bool operator==(const Particle& other) const{
        return ParticleBase<T>::operator==(other) && other._vel == _vel;
    }
};
template <typename T>
struct ParticleN : public Particle<T> {
    ALIGN_16 typename ScalarUtil<T>::ScalarVec3 _normal;
    ParticleN():_normal(0.0f,0.0f,0.0f) {}
    template<typename OTHER>
    ParticleN& copy(const OTHER& p) {
        Particle<T>::copy(p);
        CopyNormal<ParticleN,OTHER,HasNormal<OTHER>::value>::copy(*this,p);
        return *this;
    }
    virtual bool read(istream& is) {
        return Particle<T>::read(is) && readBinaryData(_normal,is).good();
    }
    virtual bool write(ostream& os) const {
        return Particle<T>::write(os) && writeBinaryData(_normal,os).good();
    }
    virtual bool operator==(const ParticleN& other) const{
        return Particle<T>::operator==(other) && other._normal == _normal;
    }
};

//============================================================================
// For VTK Writer
template <typename PS_TYPE,typename VALUE_TYPE>
struct PSetIter 
{
    PSetIter(const PS_TYPE& ps,sizeType id):_ps(ps),_id(id){}
    void operator++() {
        _id++;
    }
    bool operator!=(const PSetIter& other) const {
        return _id != other._id;
    }
    virtual VALUE_TYPE operator*() const=0;
	const PS_TYPE& _ps;
	sizeType _id;
};
template <typename PS_TYPE>
struct PSetIterPos : public PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>
{
	typedef typename PS_TYPE::Vec3Type value_type;
    PSetIterPos(const PS_TYPE& ps,sizeType id):PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>(ps,id){}
    virtual value_type operator*() const{return PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_id]._pos;}
};
template <typename PS_TYPE>
struct PSetIterPosAddVel : public PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>
{
	typedef typename PS_TYPE::Vec3Type value_type;
    PSetIterPosAddVel(const PS_TYPE& ps,sizeType id):PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>(ps,id){}
    virtual value_type operator*() const{return PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_id]._pos
                +PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::Vec3Type>::_id]._vel;}
};
template <typename PS_TYPE>
struct PSetIterNormalX : public PSetIter<PS_TYPE,typename PS_TYPE::scalarType>
{
	typedef typename PS_TYPE::scalarType value_type;
    PSetIterNormalX(const PS_TYPE& ps,sizeType id):PSetIter<PS_TYPE,typename PS_TYPE::scalarType>(ps,id){}
    virtual value_type operator*() const{return PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_id]._normal[0];}
};
template <typename PS_TYPE>
struct PSetIterNormalY : public PSetIter<PS_TYPE,typename PS_TYPE::scalarType>
{
	typedef typename PS_TYPE::scalarType value_type;
    PSetIterNormalY(const PS_TYPE& ps,sizeType id):PSetIter<PS_TYPE,typename PS_TYPE::scalarType>(ps,id){}
    virtual value_type operator*() const{return PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_id]._normal[1];}
};
template <typename PS_TYPE>
struct PSetIterNormalZ : public PSetIter<PS_TYPE,typename PS_TYPE::scalarType>
{
	typedef typename PS_TYPE::scalarType value_type;
    PSetIterNormalZ(const PS_TYPE& ps,sizeType id):PSetIter<PS_TYPE,typename PS_TYPE::scalarType>(ps,id){}
    virtual value_type operator*() const{return PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_ps[PSetIter<PS_TYPE,typename PS_TYPE::scalarType>::_id]._normal[2];}
};

typedef ParticleSetTpl<ParticleBase<scalarF> > ParticleSetBF;
typedef ParticleSetTpl<Particle<scalarF> > ParticleSetF;
typedef ParticleSetTpl<ParticleN<scalarF> > ParticleSetNF;

typedef ParticleSetTpl<ParticleBase<scalarD> > ParticleSetBD;
typedef ParticleSetTpl<Particle<scalarD> > ParticleSetD;
typedef ParticleSetTpl<ParticleN<scalarD> > ParticleSetND;

typedef ParticleSetTpl<ParticleBase<scalar> > ParticleSetB;
typedef ParticleSetTpl<Particle<scalar> > ParticleSet;
typedef ParticleSetTpl<ParticleN<scalar> > ParticleSetN;

template <typename T>
struct HasPos<ParticleBase<T> > : public boost::true_type {};

template <typename T>
struct HasPos<Particle<T> > : public boost::true_type {};
template <typename T>
struct HasVel<Particle<T> > : public boost::true_type {};

template <typename T>
struct HasPos<ParticleN<T> > : public boost::true_type {};
template <typename T>
struct HasVel<ParticleN<T> > : public boost::true_type {};
template <typename T>
struct HasNormal<ParticleN<T> > : public boost::true_type {};

PRJ_END

#endif
