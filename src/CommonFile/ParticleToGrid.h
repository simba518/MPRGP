#ifndef PARTICLE_TO_GRID_H
#define PARTICLE_TO_GRID_H

#include "GridBasic.h"
#include "ParticleSet.h"

PRJ_BEGIN

template <typename T,typename TI,typename PS_TYPE>
class ParticleToGrid
{
public:
    typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
    virtual ~ParticleToGrid() {}
    void computePhi(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& narrowBand,const T& radius) const;
    void correctPhi(const Grid<T,TI>& nodalSolid,const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const;
    void extendPhi(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const;

    void correctPhi3D(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const;
    void correctPhi2D(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const;
    void extendPhi3D(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const;
    void extendPhi2D(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const;

    static FORCE_INLINE void updatePt(T& phi,const T& dist,const T& radius) {
        phi=min<T>(phi,dist-radius);
    }
};

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::computePhi(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& narrowBand,const T& radius) const
{
    if(phi.getDim() == 3) {
        phi.init(narrowBand);
        correctPhi3D(pSet,phi,radius);
    } else if(phi.getDim() == 2) {
        phi.init(narrowBand);
        correctPhi2D(pSet,phi,radius);
    } else ASSERT(false);
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::correctPhi(const Grid<T,TI>& nodalSolid,const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const
{
    if(phi.getDim() == 3)
        correctPhi3D(pSet,phi,radius);
    else if(phi.getDim() == 2)
        correctPhi2D(pSet,phi,radius);
    else ASSERT(false);
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::extendPhi(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const
{
    if(phi.getDim() == 3)
        extendPhi3D(nodalSolid,phi);
    else if(phi.getDim() == 2)
        extendPhi2D(nodalSolid,phi);
    else ASSERT(false);
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::correctPhi3D(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const
{
    const Vec3i maxIndex=phi.getNrCell()-Vec3i::Ones();
    const sizeType n=(sizeType)pSet.size();
    const sizeType nr=std::max<sizeType>((sizeType)(radius*2.0f/phi.getCellSize().maxCoeff()),1);
    for(sizeType i=0; i<n; i++) {
        const Vec3Type& pos=pSet.get(i)._pos;
        const Vec3i index=floor(phi.getIndexFracSafe(pos));
        const Vec3i minC=compMax((Vec3i)(index-Vec3i::Constant(nr)),Vec3i::Zero());
        const Vec3i maxC=compMin((Vec3i)(index+Vec3i::Constant(nr+1)),maxIndex);

        for(sizeType x=minC.x(); x<=maxC.x(); x++)
            for(sizeType y=minC.y(); y<=maxC.y(); y++)
                for(sizeType z=minC.z(); z<=maxC.z(); z++) {
                    Vec3Type dist=phi.getPt(Vec3i(x,y,z))-pos;
                    updatePt(phi.get(Vec3i(x,y,z)),dist.norm(),radius);
                }
    }
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::correctPhi2D(const PS_TYPE& pSet,Grid<T,TI>& phi,const T& radius) const
{
    const Vec3i maxIndex=phi.getNrCell()-Vec3i::Ones();
    const sizeType n=(sizeType)pSet.size();
    const sizeType nr=std::max<sizeType>((sizeType)(radius*2.0f/phi.getCellSize().maxCoeff()),1);
    for(sizeType i=0; i<n; i++) {
        const Vec3Type& pos=pSet.get(i)._pos;
        const Vec3i index=floor(phi.getIndexFracSafe(pos));
        const Vec3i minC=compMax((Vec3i)(index-Vec3i::Constant(nr)),Vec3i::Zero());
        const Vec3i maxC=compMin((Vec3i)(index+Vec3i::Constant(nr+1)),maxIndex);

        for(sizeType x=minC.x(); x<=maxC.x(); x++)
            for(sizeType y=minC.y(); y<=maxC.y(); y++) {
                Vec3Type dist=phi.getPt(Vec3i(x,y,0))-pos;
                dist.z()=0.0f;
                updatePt(phi.get(Vec3i(x,y,0)),dist.norm(),radius);
            }
    }
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::extendPhi3D(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const
{
    Grid<T,TI> tmp=phi;

    const T thres=phi.getCellSize().maxCoeff()*0.5f;
    const Vec3i nrCell=phi.getNrCell();
    for(sizeType x=0; x < nrCell.x(); ++x) {
        for(sizeType y=0; y < nrCell.y(); ++y)
            for(sizeType z=0; z < nrCell.z(); ++z) {
                if(phi.get(Vec3i(x,y,z)) < thres) {
                    T solidPhi=0.125f*(nodalSolid.get(Vec3i(x,y,z)+Vec3i(0,0,0))+nodalSolid.get(Vec3i(x,y,z)+Vec3i(0,0,1))+
                                       nodalSolid.get(Vec3i(x,y,z)+Vec3i(0,1,0))+nodalSolid.get(Vec3i(x,y,z)+Vec3i(0,1,1))+
                                       nodalSolid.get(Vec3i(x,y,z)+Vec3i(1,0,0))+nodalSolid.get(Vec3i(x,y,z)+Vec3i(1,0,1))+
                                       nodalSolid.get(Vec3i(x,y,z)+Vec3i(1,1,0))+nodalSolid.get(Vec3i(x,y,z)+Vec3i(1,1,1)));
                    if(solidPhi < 0.0f)
                        tmp.get(Vec3i(x,y,z))=-thres;
                }
            }
    }
    phi=tmp;
}

template <typename T,typename TI,typename PS_TYPE>
void ParticleToGrid<T,TI,PS_TYPE>::extendPhi2D(const Grid<T,TI>& nodalSolid,Grid<T,TI>& phi) const
{
    Grid<T,TI> tmp=phi;

    const T thres=phi.getCellSize().maxCoeff()*0.5f;
    const Vec3i nrCell=phi.getNrCell();
    for(sizeType x=0; x < nrCell.x(); ++x) {
        for(sizeType y=0; y < nrCell.y(); ++y) {
            if(phi.get(Vec3i(x,y,0)) < thres) {
                T solidPhi=0.25f*(nodalSolid.get(Vec3i(x,y,0)+Vec3i(0,0,0))+nodalSolid.get(Vec3i(x,y,0)+Vec3i(0,1,0))+
                                  nodalSolid.get(Vec3i(x,y,0)+Vec3i(1,0,0))+nodalSolid.get(Vec3i(x,y,0)+Vec3i(1,1,0)));
                if(solidPhi < 0.0f)
                    tmp.get(Vec3i(x,y,0))=-thres;
            }
        }
    }
    phi=tmp;
}

PRJ_END

#endif
