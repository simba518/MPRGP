#ifndef VELOCITY_SET_H
#define VELOCITY_SET_H

#include "GridBasic.h"

PRJ_BEGIN

template <typename T,typename TI,typename TV=vector<T,Eigen::aligned_allocator<T> > >
class VelocitySet : public HasMagic
{
public:
    VelocitySet():HasMagic(0xFFCCFFCCFFCCFFCC) {}
    virtual ~VelocitySet() {}
    void read(istream& is) {
        if(!HasMagic::readMagic(is))
            return false;

        sizeType n;
        _frames.swap(std::vector< MACGrid<T,TI,TV> >());
        readBinaryData(n,is);

        _frames.resize(n);
        for(sizeType i=0; i<n; i++)
            _frames[i].read(is);

        return is.good();
    }
    bool write(ostream& os) const {
        if(!HasMagic::writeMagic(os))
            return false;

        sizeType n=_frames.size();
        writeBinaryData(n,os);

        for(sizeType i=0; i<n; i++)
            _frames[i].write(os);

        return os.good();
    }
    void makeUnit(MACGrid<T,TI,TV>& curr) const {
        T total=0.0f;

        if(curr.getDim() >= 1)
            total+=curr.getGu().dot(curr.getGu());
        if(curr.getDim() >= 2)
            total+=curr.getGv().dot(curr.getGv());
        if(curr.getDim() >= 3)
            total+=curr.getGw().dot(curr.getGw());
        total=sqrt(total);

        if(curr.getDim() >= 1)
            curr.getGu().mul(1.0f/std::max<T>(total,ScalarUtil<T>::scalar_eps));
        if(curr.getDim() >= 2)
            curr.getGv().mul(1.0f/std::max<T>(total,ScalarUtil<T>::scalar_eps));
        if(curr.getDim() >= 3)
            curr.getGw().mul(1.0f/std::max<T>(total,ScalarUtil<T>::scalar_eps));
    }
    void makeUnit() {
        sizeType n=_frames.size();
        for(sizeType i=0; i<n; i++) {
            MACGrid<T,TI,TV>& curr=_frames[i];
            makeUnit(curr);
        }
    }
    void applyMask(Grid<T,TI,TV>& curr,const Grid<T,TI,TV>& mask) const {
        const Vec3i nrPoint=curr.getNrPoint();
        OMP_PARALLEL_FOR_
        for(sizeType x=0; x<nrPoint.x(); x++) {
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++)
                    if(mask.get(Vec3i(x,y,z)) <= 0.0f)
                        curr.get(Vec3i(x,y,z))=0.0f;
        }
    }
    void applyMask(MACGrid<T,TI,TV>& mask) {
        sizeType n=_frames.size();
        for(sizeType i=0; i<n; i++) {
            MACGrid<T,TI,TV>& curr=_frames[i];
            if(curr.getDim() >= 1)
                applyMask(curr.getGu(),mask.getGu());
            if(curr.getDim() >= 2)
                applyMask(curr.getGv(),mask.getGv());
            if(curr.getDim() >= 3)
                applyMask(curr.getGw(),mask.getGw());
        }
    }
    void makeOrtho(MACGrid<T,TI,TV>& curr,const sizeType& i) const {
        for(sizeType j=0; j<i; j++) {
            T dot=0.0f;

            if(curr.getDim() >= 1)
                dot+=curr.getGu().dot(_frames[j].getGu());
            if(curr.getDim() >= 2)
                dot+=curr.getGv().dot(_frames[j].getGv());
            if(curr.getDim() >= 3)
                dot+=curr.getGw().dot(_frames[j].getGw());

            curr.addScaled(_frames[j],-dot);
        }
    }
    void makeOrtho() {
        sizeType n=_frames.size();
        for(sizeType i=0; i<n; i++) {
            MACGrid<T,TI,TV>& curr=_frames[i];
            makeOrtho(curr,i);
            makeUnit(curr);

            if(i >= 1)
                INFOV("Unity Check: %f, Ortho Check: %f",
                      curr.dot(curr),curr.dot(_frames[i-1]));
        }
    }
    void mul(const Eigen::Matrix<T,-1,1>& coef,MACGrid<T,TI,TV>& result) const {
        sizeType n=_frames.size();

        result.init(0.0f);
        for(sizeType i=0; i<n; i++)
            result.addScaled(_frames[i],coef(i,0));
    }
    std::vector< MACGrid<T,TI,TV> > _frames;
};

PRJ_END

#endif
