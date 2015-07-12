#ifndef POISSON_DISK_SAMPLING_H
#define POISSON_DISK_SAMPLING_H

#include "Config.h"
#include "MathBasic.h"
#include "GridBasic.h"
#include "ParticleSet.h"

PRJ_BEGIN

class PoissonDiskSampling
{
public:
    PoissonDiskSampling(bool initial) {
        _tSurface=30;
        _tVolume=initial ? 30 : 8;
        _e=1.085f;
        _tRelax=50;
        _kSurface=5;
        _kVolume=30;
    }
    virtual ~PoissonDiskSampling() {}
public:
    void buildGrid(const Grid<scalar,scalar>& constraint,
                   Grid<sizeType,scalar>& tag,const scalar& r) {
        if(constraint.getDim() == 3)
            buildGrid3D(constraint,tag,r);
        else
            buildGrid2D(constraint,tag,r);
    }
    void buildGrid2D(const Grid<scalar,scalar>& constraint,
                     Grid<sizeType,scalar>& tag,const scalar& r) {
        //construct assistant grid
        Vec3 szCell=Vec3::Constant(r/sqrt(2.0f));
        BBox<scalar> bb=constraint.getBB();
        Vec3 extent=bb.getExtent();
        Vec3i nrCell=ceil((Vec3)(extent.array()/szCell.array()).matrix());
        extent=Vec3(nrCell.x()*szCell.x(),nrCell.y()*szCell.y(),0.0f);

        //create assistant grid
        tag.reset(nrCell,BBox<scalar>(bb._minC,bb._minC+extent),-1,false);
    }
    void buildGrid3D(const Grid<scalar,scalar>& constraint,
                     Grid<sizeType,scalar>& tag,const scalar& r) {
        Vec3 szCell=Vec3::Constant(r/sqrt(3.0f));
        BBox<scalar> bb=constraint.getBB();
        Vec3 extent=bb.getExtent();
        Vec3i nrCell=ceil((Vec3)(extent.array()/szCell.array()).matrix());
        extent=Vec3(nrCell.x()*szCell.x(),nrCell.y()*szCell.y(),nrCell.z()*szCell.z());

        //create assistant grid
        tag.reset(nrCell,BBox<scalar>(bb._minC,bb._minC+extent),-1,false);
    }
    template<typename P_TYPE>
    void fillParticle(Grid<sizeType,scalar>& tag,const ParticleSetTpl<P_TYPE>& pSet) {
        const sizeType n=pSet.size();
        for(sizeType i=0; i<n; i++)
            tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos)))=i;
    }
public:
    template<typename P_TYPE>
    void surfaceSample(const Grid<scalar,scalar>& constraint,
                       Grid<sizeType,scalar>& tag,
                       const scalar& r,const scalar& phi0,const scalar& phi1,
                       ParticleSetTpl<P_TYPE>& pSet) {
        if(constraint.getDim() == 3)
            surfaceSample3D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet);
        else
            surfaceSample2D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet);
    }
    template<typename P_TYPE>
    void surfaceSample2D(const Grid<scalar,scalar>& constraint,
                         Grid<sizeType,scalar>& tag,
                         const scalar& r,const scalar& phi0,const scalar& phi1,
                         ParticleSetTpl<P_TYPE>& pSet) {
        //construct assistant grid
        Vec3i nrCell=tag.getNrCell();
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);
        safeBox._minC.z()=-EPS;
        safeBox._maxC.z()=EPS;

        //start algorithm
        P_TYPE p;
        P_TYPE pTmp;

        sizeType seed=rand();
        for(sizeType x=0; x<nrCell.x(); x++)
            for(sizeType y=0; y<nrCell.y(); y++) {
                const Vec3 minPt=tag.getPt(Vec3i(x  ,y  ,0));
                const Vec3 maxPt=tag.getPt(Vec3i(x+1,y+1,0));

                if(!changeSign2D(constraint.sampleSafe2D(minPt),
                                 constraint.sampleSafe2D(tag.getPt(Vec3i(x+1,y  ,0))),
                                 constraint.sampleSafe2D(tag.getPt(Vec3i(x  ,y+1,0))),
                                 constraint.sampleSafe2D(maxPt) ))
                    continue;

                bool found=false;
                for(sizeType att=0; att<_tSurface; att++) {
                    p._pos=randPt(seed,minPt,maxPt);
                    project(constraint,p._pos,phi0,phi1);
                    if(safeBox.contain(p._pos) && !checkNeigh2D(tag,pSet,p._pos,r)) {
                        ASSERT(tag.get(floor(tag.getIndexFrac(p._pos))) == -1)
                        tag.get(floor(tag.getIndexFrac(p._pos)))=pSet.size();
                        pSet.addParticle(p);
                        found=true;
                        break;
                    }
                }

                while(found) {
                    Vec3 normal=constraint.sampleSafe2DGrad(p._pos);
                    if(normal.norm() > EPS)
                        normal.normalize();
                    pTmp._pos=p._pos+randTangent2D(seed,normal)*_e*r;
                    project(constraint,pTmp._pos,phi0,phi1);
                    if(safeBox.contain(pTmp._pos) && !checkNeigh2D(tag,pSet,pTmp._pos,r)) {
                        ASSERT(tag.get(floor(tag.getIndexFrac(pTmp._pos))) == -1)
                        tag.get(floor(tag.getIndexFrac(pTmp._pos)))=pSet.size();
                        pSet.addParticle(pTmp);
                    } else found=false;
                }
            }
    }
    template<typename P_TYPE>
    void surfaceSample3D(const Grid<scalar,scalar>& constraint,
                         Grid<sizeType,scalar>& tag,
                         const scalar& r,const scalar& phi0,const scalar& phi1,
                         ParticleSetTpl<P_TYPE>& pSet) {
        //construct assistant grid
        Vec3i nrCell=tag.getNrCell();
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);

        //start algorithm
        P_TYPE p;
        P_TYPE pTmp;

        sizeType seed=rand();
        for(sizeType x=0; x<nrCell.x(); x++)
            for(sizeType y=0; y<nrCell.y(); y++)
                for(sizeType z=0; z<nrCell.z(); z++) {
                    const Vec3 minPt=tag.getPt(Vec3i(x  ,y  ,z));
                    const Vec3 maxPt=tag.getPt(Vec3i(x+1,y+1,z));

                    if(!changeSign3D(constraint.sampleSafe3D(minPt),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x+1,y  ,z  ))),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x  ,y+1,z  ))),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x+1,y+1,z  ))),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x  ,y  ,z+1))),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x+1,y  ,z+1))),
                                     constraint.sampleSafe3D(tag.getPt(Vec3i(x  ,y+1,z+1))),
                                     constraint.sampleSafe3D(maxPt )))
                        continue;

                    bool found=false;
                    for(sizeType att=0; att<_tSurface; att++) {
                        p._pos=randPt(seed,minPt,maxPt);
                        project(constraint,p._pos,phi0,phi1);
                        if(safeBox.contain(p._pos) && !checkNeigh3D(tag,pSet,p._pos,r)) {
                            ASSERT(tag.get(floor(tag.getIndexFrac(p._pos))) == -1)
                            tag.get(floor(tag.getIndexFrac(p._pos)))=pSet.size();
                            pSet.addParticle(p);
                            found=true;
                            break;
                        }
                    }

                    while(found) {
                        Vec3 normal=constraint.sampleSafe3DGrad(p._pos);
                        if(normal.norm() > EPS)
                            normal.normalize();
                        pTmp._pos=p._pos+randTangent3D(seed,normal)*_e*r;
                        project(constraint,pTmp._pos,phi0,phi1);
                        if(safeBox.contain(pTmp._pos) && !checkNeigh3D(tag,pSet,pTmp._pos,r)) {
                            ASSERT(tag.get(floor(tag.getIndexFrac(pTmp._pos))) == -1)
                            tag.get(floor(tag.getIndexFrac(pTmp._pos)))=pSet.size();
                            pSet.addParticle(pTmp);
                        } else found=false;
                    }
                }
    }
public:
    template<typename P_TYPE>
    void sampleVolume(const Grid<scalar,scalar>& constraint,
                      Grid<sizeType,scalar>& tag,
                      const scalar& r,const scalar& phi0,const scalar& phi1,
                      ParticleSetTpl<P_TYPE>& pSet,bool restart) {
        if(constraint.getDim() == 3)
            sampleVolume3D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet,restart);
        else
            sampleVolume2D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet,restart);
    }
    template<typename P_TYPE>
    void sampleVolume2D(const Grid<scalar,scalar>& constraint,
                        Grid<sizeType,scalar>& tag,
                        const scalar& r,const scalar& phi0,const scalar& phi1,
                        ParticleSetTpl<P_TYPE>& pSet,bool restart) {
        //construct assistant grid
        Vec3i nrCell=tag.getNrCell();
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);
        safeBox._minC.z()=-EPS;
        safeBox._maxC.z()=EPS;

        //start algorithm
        std::vector<sizeType> activeSet;
        sizeType seed=rand();

        //generate initial sample
        P_TYPE p;
        P_TYPE pTmp;
        if(restart) {
            while(true) {
                p._pos=randPt(seed,safeBox._minC,safeBox._maxC);
                p._pos.z()=0.0f;
                if(!needProject(constraint.sampleSafe2D(p._pos),phi0,phi1)) {
                    pSet.addParticle(p);
                    tag.get(floor(tag.getIndexFracSafe(p._pos)))=0;
                    activeSet.push_back(0);
                    break;
                }
            }
        } else {
            for(sizeType i=0; i<pSet.size(); i++)
                activeSet.push_back(i);
        }

        //update active set
        while(!activeSet.empty()) {
            const sizeType i=randId(seed,0,(sizeType)activeSet.size()-1);
            p=pSet.get(activeSet[i]);

            bool found=false;
            for(sizeType att=0; att<_tVolume; att++) {
                pTmp._pos=p._pos+randSphere2D(seed,r);
                pTmp._pos.z()=0.0f;
                if(safeBox.contain(pTmp._pos) && !needProject(constraint.sampleSafe2D(pTmp._pos),phi0,phi1) && !checkNeigh2D(tag,pSet,pTmp._pos,r)) {
                    tag.get((Vec3i)floor(tag.getIndexFracSafe(pTmp._pos)))=pSet.size();
                    activeSet.push_back(pSet.size());
                    pSet.addParticle(pTmp);
                    found=true;
                    break;
                }
            }

            if(!found) {
                activeSet[i]=activeSet.back();
                activeSet.pop_back();
            }
        }
    }
    template<typename P_TYPE>
    void sampleVolume3D(const Grid<scalar,scalar>& constraint,
                        Grid<sizeType,scalar>& tag,
                        const scalar& r,const scalar& phi0,const scalar& phi1,
                        ParticleSetTpl<P_TYPE>& pSet,bool restart) {
        //construct assistant grid
        Vec3i nrCell=tag.getNrCell();
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);

        //start algorithm
        std::vector<sizeType> activeSet;
        sizeType seed=rand();

        //generate initial sample
        P_TYPE p;
        P_TYPE pTmp;
        if(restart) {
            while(true) {
                p._pos=randPt(seed,safeBox._minC,safeBox._maxC);
                if(!needProject(constraint.sampleSafe3D(p._pos),phi0,phi1)) {
                    pSet.addParticle(p);
                    tag.get(floor(tag.getIndexFracSafe(p._pos)))=0;
                    activeSet.push_back(0);
                    break;
                }
            }
        } else {
            for(sizeType i=0; i<pSet.size(); i++)
                activeSet.push_back(i);
        }

        //update active set
        while(!activeSet.empty()) {
            const sizeType i=randId(seed,0,(sizeType)activeSet.size()-1);
            p=pSet.get(activeSet[i]);

            bool found=false;
            for(sizeType att=0; att<_tVolume; att++) {
                pTmp._pos=p._pos+randSphere3D(seed,r);
                if(safeBox.contain(pTmp._pos) && !needProject(constraint.sampleSafe3D(pTmp._pos),phi0,phi1) && !checkNeigh3D(tag,pSet,pTmp._pos,r)) {
                    tag.get((Vec3i)floor(tag.getIndexFracSafe(pTmp._pos)))=pSet.size();
                    activeSet.push_back(pSet.size());
                    pSet.addParticle(pTmp);
                    found=true;
                    break;
                }
            }

            if(!found) {
                activeSet[i]=activeSet.back();
                activeSet.pop_back();
            }
        }
    }
public:
    template<typename P_TYPE>
    void relaxSample(const Grid<scalar,scalar>& constraint,
                     Grid<sizeType,scalar>& tag,
                     const scalar& r,const scalar& phi0,const scalar& phi1,
                     ParticleSetTpl<P_TYPE>& pSet,
                     const sizeType& start,const sizeType& end) {
        if(constraint.getDim() == 3)
            relaxSample3D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet,start,end);
        else
            relaxSample2D<P_TYPE>(constraint,tag,r,phi0,phi1,pSet,start,end);
    }
    template<typename P_TYPE>
    void relaxSample2D(const Grid<scalar,scalar>& constraint,
                       Grid<sizeType,scalar>& tag,
                       const scalar& r,const scalar& phi0,const scalar& phi1,
                       ParticleSetTpl<P_TYPE>& pSet,
                       const sizeType& start,const sizeType& end) {
        ASSERT(start<=end && start>=0 && end<=pSet.size());
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);
        safeBox._minC.z()=-EPS;
        safeBox._maxC.z()=EPS;

        bool isSurface=(phi0 == phi1);
        sizeType k=isSurface ? _kSurface : _kVolume;
        sizeType seed=rand();

        for(sizeType kk=0; kk<k; kk++) {
            for(sizeType i=start; i<end; i++) {
                scalar d=minDist2D(tag,pSet,i,pSet.get(i)._pos,r*2.0f);
                scalar dd;
                scalar tau;
                ASSERT(d >= r)

                P_TYPE pNew=pSet.get(i);
                P_TYPE pCand;
                //adjust
                for(sizeType tt=0; tt<_tRelax; tt++) {
                    tau=(_tRelax-tt)/(scalar)_tRelax;
                    if(isSurface) {
                        Vec3 normal=constraint.sampleSafe2DGrad(pNew._pos);
                        if(normal.norm() > EPS)
                            normal.normalize();
                        pCand._pos=pNew._pos+randTangent2D(seed,normal)*(tau*r);
                    } else {
                        pCand._pos=pNew._pos+randDir2D(seed)*(tau*r);
                    }
                    project(constraint,pCand._pos,phi0,phi1);
                    dd=minDist2D(tag,pSet,i,pCand._pos,r*2.0f);
                    if(safeBox.contain(pCand._pos) && dd > d) {
                        d=dd;
                        pNew=pCand;
                    }
                }
                //reset
                tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos)))=-1;
                pSet.get(i)=pNew;
                ASSERT(tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos))) == -1)
                tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos)))=i;
            }
        }
    }
    template<typename P_TYPE>
    void relaxSample3D(const Grid<scalar,scalar>& constraint,
                       Grid<sizeType,scalar>& tag,
                       const scalar& r,const scalar& phi0,const scalar& phi1,
                       ParticleSetTpl<P_TYPE>& pSet,
                       const sizeType& start,const sizeType& end) {
        ASSERT(start<=end && start>=0 && end<=pSet.size());
        BBox<scalar> safeBox=tag.getBB();
        safeBox.enlarge(-EPS);

        bool isSurface=(phi0 == phi1);
        sizeType k=isSurface ? _kSurface : _kVolume;
        sizeType seed=rand();

        for(sizeType kk=0; kk<k; kk++) {
            for(sizeType i=start; i<end; i++) {
                scalar d=minDist3D(tag,pSet,i,pSet.get(i)._pos,r*2.0f);
                scalar dd;
                scalar tau;
                ASSERT(d >= r)

                P_TYPE pNew=pSet.get(i);
                P_TYPE pCand;
                //adjust
                for(sizeType tt=0; tt<_tRelax; tt++) {
                    tau=(_tRelax-tt)/(scalar)_tRelax;
                    if(isSurface) {
                        Vec3 normal=constraint.sampleSafe3DGrad(pNew._pos);
                        if(normal.norm() > EPS)
                            normal.normalize();
                        pCand._pos=pNew._pos+randTangent3D(seed,normal)*(tau*r);
                    } else {
                        pCand._pos=pNew._pos+randDir3D(seed)*(tau*r);
                    }
                    project(constraint,pCand._pos,phi0,phi1);
                    dd=minDist3D(tag,pSet,i,pCand._pos,r*2.0f);
                    if(safeBox.contain(pCand._pos) && dd > d) {
                        d=dd;
                        pNew=pCand;
                    }
                }
                //reset
                tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos)))=-1;
                pSet.get(i)=pNew;
                ASSERT(tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos))) == -1)
                tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos)))=i;
            }
        }
    }
public:
    template<typename P_TYPE>
    void checkValidity(const Grid<sizeType,scalar>& tag,
                       const ParticleSetTpl<P_TYPE>& pSet) {
        for(sizeType i=0; i<pSet.size(); i++) {
            ASSERT(tag.get(floor(tag.getIndexFracSafe(pSet.get(i)._pos))) == i)
        }
    }
protected:
    static FORCE_INLINE Vec3 randPt(sizeType& seed,const Vec3& minPt,const Vec3& maxPt) {
        Vec3 ret(RandHash::randhash(seed,minPt.x(),maxPt.x()),
                 RandHash::randhash(seed+1,minPt.y(),maxPt.y()),
                 RandHash::randhash(seed+2,minPt.z(),maxPt.z()));
        seed+=3;
        return ret;
    }
    static FORCE_INLINE bool changeSign2D(const scalar& a,const scalar& b,const scalar& c,const scalar& d) {
        return (a<0.0f || b<0.0f || c<0.0f || d<0.0f) &&
               (a>0.0f || b>0.0f || c>0.0f || d>0.0f);
    }
    static FORCE_INLINE bool changeSign3D(const scalar& a,const scalar& b,const scalar& c,const scalar& d,
                                          const scalar& e,const scalar& f,const scalar& g,const scalar& h) {
        return (a<0.0f || b<0.0f || c<0.0f || d<0.0f || e<0.0f || f<0.0f || g<0.0f || h<0.0f) &&
               (a>0.0f || b>0.0f || c>0.0f || d>0.0f || e>0.0f || f>0.0f || g>0.0f || h>0.0f);
    }
    template<typename P_TYPE>
    static FORCE_INLINE bool checkNeigh2D(const Grid<sizeType,scalar>& tag,const ParticleSetTpl<P_TYPE>& pSet,const Vec3& pos,const scalar& r) {
        Vec3i base=floor(tag.getIndexFracSafe(pos));
        for(sizeType x=base.x()-2; x<=base.x()+2; x++)
            for(sizeType y=base.y()-2; y<=base.y()+2; y++) {
                sizeType id=tag.getSafe(Vec3i(x,y,0));
                if(id>=0 && (pSet.get(id)._pos-pos).norm() < r)
                    return true;
            }
        return false;
    }
    template<typename P_TYPE>
    static FORCE_INLINE bool checkNeigh3D(const Grid<sizeType,scalar>& tag,const ParticleSetTpl<P_TYPE>& pSet,const Vec3& pos,const scalar& r) {
        Vec3i base=floor(tag.getIndexFracSafe(pos));
        for(sizeType x=base.x()-2; x<=base.x()+2; x++)
            for(sizeType y=base.y()-2; y<=base.y()+2; y++)
                for(sizeType z=base.z()-2; z<=base.z()+2; z++) {
                    sizeType id=tag.getSafe(Vec3i(x,y,z));
                    if(id>=0 && (pSet.get(id)._pos-pos).norm() < r)
                        return true;
                }
        return false;
    }
    static FORCE_INLINE void project(const Grid<scalar,scalar>& phi,Vec3& pos,const scalar& phi0,const scalar& phi1) {
        Vec3 grad=phi.sampleSafeGrad(pos);
        if(grad.norm() > EPS)
            grad.normalize();
        scalar val=phi.sampleSafe(pos);
        if(val > phi1)
            pos-=grad*(val-phi1);
        else if(val < phi0)
            pos+=grad*(phi0-val);
    }
    static FORCE_INLINE bool needProject(const scalar& phi,const scalar& phi0,const scalar& phi1) {
        return phi < phi0 || phi > phi1;
    }
    static FORCE_INLINE Vec3 randDir2D(sizeType& seed) {
        scalar alpha=RandHash::randhash(seed++,0.0f,M_PI*2.0f);
        Vec3 randVec(cos(alpha),sin(alpha),0.0f);	//a random direction in unit sphere
        return randVec;
    }
    static FORCE_INLINE Vec3 randDir3D(sizeType& seed) {
        scalar alpha=RandHash::randhash(seed++,0.0f,M_PI*2.0f);
        scalar beta=RandHash::randhash(seed++,-M_PI*0.5f,M_PI*0.5f);
        Vec3 randVec(cos(alpha)*cos(beta),sin(alpha)*cos(beta),sin(beta));	//a random direction in unit sphere
        return randVec;
    }
    static FORCE_INLINE Vec3 randTangent2D(sizeType& seed,const Vec3& normal) {
        scalar coef=rand()/(scalar)RAND_MAX > 0.5f ? 1.0f : -1.0f;
        Vec3 ret(-normal.y(),normal.x(),0.0f);
        if(ret.norm() > EPS)
            ret.normalize();
        return ret*coef;
    }
    static FORCE_INLINE Vec3 randTangent3D(sizeType& seed,const Vec3& normal) {
        Vec3 randVec=randDir3D(seed);
        Vec3 ret(randVec-normal*randVec.dot(normal));			//project to tangent space
        if(ret.norm() > EPS)
            ret.normalize();
        return ret;
    }
    static FORCE_INLINE Vec3 randSphere2D(sizeType& seed,const scalar& r) {
        return randDir2D(seed)*(r+RandHash::randhash(seed++,0.0f,r));
    }
    static FORCE_INLINE Vec3 randSphere3D(sizeType& seed,const scalar& r) {
        return randDir3D(seed)*(r+RandHash::randhash(seed++,0.0f,r));
    }
    static FORCE_INLINE sizeType randId(sizeType& seed,const sizeType& a,const sizeType& b) {
        sizeType diff=RandHash::randhashI(seed++)%(b-a+1);
        if(diff<0)
            diff+=(b-a+1);
        return diff+a;
    }
    template<typename P_TYPE>
    static FORCE_INLINE scalar minDist2D(const Grid<sizeType,scalar>& tag,const ParticleSetTpl<P_TYPE>& pSet,const sizeType i,const Vec3& pos,const scalar& r2) {
        scalar dist=r2;
        Vec3i base=floor(tag.getIndexFracSafe(pos));
        for(sizeType x=base.x()-2; x<=base.x()+2; x++)
            for(sizeType y=base.y()-2; y<=base.y()+2; y++) {
                sizeType id=tag.getSafe(Vec3i(x,y,0));
                if(id >= 0 && id != i)
                    dist=min<scalar>(dist,(pSet.get(id)._pos-pos).norm());
            }
        return dist;
    }
    template<typename P_TYPE>
    static FORCE_INLINE scalar minDist3D(const Grid<sizeType,scalar>& tag,const ParticleSetTpl<P_TYPE>& pSet,const sizeType& i,const Vec3& pos,const scalar& r2) {
        scalar dist=r2;
        Vec3i base=floor(tag.getIndexFracSafe(pos));
        for(sizeType x=base.x()-2; x<=base.x()+2; x++)
            for(sizeType y=base.y()-2; y<=base.y()+2; y++)
                for(sizeType z=base.z()-2; z<=base.z()+2; z++) {
                    sizeType id=tag.getSafe(Vec3i(x,y,z));
                    if(id >= 0 && id != i)
                        dist=min<scalar>(dist,(pSet.get(id)._pos-pos).norm());
                }
        return dist;
    }
protected:
    sizeType _tSurface;
    sizeType _tVolume;
    scalar _e;
    sizeType _tRelax;
    sizeType _kSurface;
    sizeType _kVolume;
};

PRJ_END

#endif
