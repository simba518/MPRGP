#ifndef LAPLACE_POISSON_OP_H
#define LAPLACE_POISSON_OP_H

#include "GridOp.h"
#include "solvers/LinearSolver.h"

PRJ_BEGIN

template<typename T,typename TI,typename TV=vector<T,Eigen::aligned_allocator<T> > >
class laplacePoissonOp
{
public:
    typedef Grid<T,TI,TV> GridType;
    typedef MACGrid<T,TI,TV> MACGridType;
    typedef Grid<unsigned char,TI> TagGridType;
    typedef typename GridType::indexType IndexType;
    typedef Kernel<scalarD> KT;
    typedef typename KT::Vec Vec;
    void laplacePoisson(GridType& grd,const TagGridType* tag=NULL,const GridType* centerWeight=NULL,const MACGridType* faceWeight=NULL,const GridType* reg=NULL,const T& regCoef=0.0f) {
        const IndexType cellSz=grd.getCellSize();
        const Vec3i nrPoint=grd.getNrPoint();
        const sizeType sz=nrPoint.prod();
        Vec3i stride;
        T weight;
        if(grd.getDim() == 3) {
            stride=Vec3i(nrPoint.y()*nrPoint.z(),nrPoint.z(),1);
            weight=cellSz.x()*cellSz.y()*cellSz.z();
        } else {
            stride=Vec3i(nrPoint.y(),1,0);
            weight=cellSz.x()*cellSz.y();
        }

        _matrix.resize(sz);
        Vec rhs(sz);

        _matrix.zero();
        rhs.setZero();

        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const Vec3i id(x,y,z);
                    const sizeType off=id.dot(stride);

                    if(tag && tag->get(id) == 0) {
                        rhs[off]=grd.get(id);
                        _matrix.setElement(off,off,1.0f);
                        continue;
                    }

                    const T wCtr=(centerWeight ? centerWeight->get(id) : 1.0f)*weight;
                    const T wL=(faceWeight ? faceWeight->getGu().get(id) : 1.0f)*weight;
                    const T wR=(faceWeight ? faceWeight->getGu().get(id+Vec3i(1,0,0)) : 1.0f)*weight;
                    const T wB=(faceWeight ? faceWeight->getGv().get(id) : 1.0f)*weight;
                    const T wT=(faceWeight ? faceWeight->getGv().get(id+Vec3i(0,1,0)) : 1.0f)*weight;
                    const T wP=(grd.getDim() == 3 && faceWeight ? faceWeight->getGw().get(id) : 1.0f)*weight;
                    const T wN=(grd.getDim() == 3 && faceWeight ? faceWeight->getGw().get(id+Vec3i(0,0,1)) : 1.0f)*weight;

                    if(reg) {
                        rhs[off]+=(regCoef*reg->get(id))*wCtr;
                        _matrix.addToElement(off,off,regCoef*wCtr);
                    }

                    if(grd.isSafeIndex(id+Vec3i(1,0,0)) && grd.isSafeIndex(id-Vec3i(1,0,0))) {
                        _matrix.addToElement(off,off,wR);
                        if(!tag || tag->get(id+Vec3i(1,0,0)) == 1)
                            _matrix.addToElement(off,off+stride.x(),-wR);
                        else rhs[off]+=wR*grd.get(id+Vec3i(1,0,0));

                        _matrix.addToElement(off,off,wL);
                        if(!tag || tag->get(id-Vec3i(1,0,0)) == 1)
                            _matrix.addToElement(off,off-stride.x(),-wL);
                        else rhs[off]+=wL*grd.get(id-Vec3i(1,0,0));
                    }

                    if(grd.isSafeIndex(id+Vec3i(0,1,0)) && grd.isSafeIndex(id-Vec3i(0,1,0))) {
                        _matrix.addToElement(off,off,wT);
                        if(!tag || tag->get(id+Vec3i(0,1,0)) == 1)
                            _matrix.addToElement(off,off+stride.y(),-wT);
                        else rhs[off]+=wT*grd.get(id+Vec3i(0,1,0));

                        _matrix.addToElement(off,off,wB);
                        if(!tag || tag->get(id-Vec3i(0,1,0)) == 1)
                            _matrix.addToElement(off,off-stride.y(),-wB);
                        else rhs[off]+=wB*grd.get(id-Vec3i(0,1,0));
                    }

                    if(grd.isSafeIndex(id+Vec3i(0,0,1)) && grd.isSafeIndex(id-Vec3i(0,0,1)) && grd.getDim() == 3) {
                        _matrix.addToElement(off,off,wN);
                        if(!tag || tag->get(id+Vec3i(0,0,1)) == 1)
                            _matrix.addToElement(off,off+stride.z(),-wN);
                        else rhs[off]+=wN*grd.get(id+Vec3i(0,0,1));

                        _matrix.addToElement(off,off,wP);
                        if(!tag || tag->get(id-Vec3i(0,0,1)) == 1)
                            _matrix.addToElement(off,off-stride.z(),-wP);
                        else rhs[off]+=wP*grd.get(id-Vec3i(0,0,1));
                    }
                }

        ASSERT(_matrix.isDiagonalDominant());
        //ASSERT(_matrix.isSymmetric(1E-9f));

        Vec X=rhs;
        _solver.setSolverParameters(1E-5f,1000);
        _solver.setMatrix(_matrix);
        _solver.solve(rhs,X);

        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++)
                for(sizeType z=0; z<nrPoint.z(); z++) {
                    const Vec3i id(x,y,z);
                    const sizeType off=id.dot(stride);
                    if(!tag || tag->get(id) == 1)
                        grd.get(id)=X[off];
                }
    }
protected:
    BiCGStabSolver<scalarD,KT> _solver;
    SparseMatrix<scalarD,KT> _matrix;
};

PRJ_END

#endif
