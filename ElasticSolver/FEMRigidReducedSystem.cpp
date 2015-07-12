#include "FEMRigidReducedSystem.h"
#include "FEMMesh.h"
#include "FEMUtils.h"
#include "FEMRotationUtil.h"

USE_PRJ_NAMESPACE

//M calculator for rigid coupled model
//we ignore all the centrifugal and coriolis force
RigidReducedMCalculator::RigidReducedMCalculator(FEMReducedSystem& sys,boost::shared_ptr<ReducedMCalculator> innerMCalc)
    :ReducedMCalculator(sys),_innerMCalc(innerMCalc) {}
void RigidReducedMCalculator::calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U)
{
    const SparseReducedBasis& basis=*(_body._basis);
    sizeType rows=basis._U.rows();
    sizeType cols=basis._U.cols();
    //rotation part
    Mat3d DRDWX,DRDWY,DRDWZ;
    Mat3d R=expWGrad<Mat3d>(_body._W.cast<scalarD>(),&DRDWX,&DRDWY,&DRDWZ);
    //calculate U
    _UDense.setZero(rows,cols+6);
    {
        FEMSystemMatrix innerU;
        innerU.reset(rows,cols,true);
        Vec innerS=S.block(0,0,_innerMCalc->size(),1);
        _innerMCalc->calculateM(innerS,M,Vec::Zero(_innerMCalc->sizeL()),RHS,0,0,&innerU);

        Matd innerUDense;
        innerU.getDense(innerUDense);
        _UDense.block(0,0,rows,cols)=innerUDense.block(0,0,rows,cols);
    }
    //rigid motion part
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<_body.nrV(); i++) {
        const FEMVertex& v=_body.getV(i);
        Vec3d dpos=R.transpose()*(v._pos-_body._X0).cast<scalarD>();
        _UDense.block(i*3,0,3,cols)=R*_UDense.block(i*3,0,3,cols);
        _UDense.block<3,3>(i*3,cols).setIdentity();
        _UDense.block<3,1>(i*3,cols+3)=DRDWX*dpos;
        _UDense.block<3,1>(i*3,cols+4)=DRDWY*dpos;
        _UDense.block<3,1>(i*3,cols+5)=DRDWZ*dpos;
    }
    RHS.block(0,0,_UDense.cols(),1)+=_UDense.transpose()*basis._M*MRHSL;
    if(U)U->addBlock(r,c,_UDense);
    //update UTMU upper left
    Matd tmp=M;
    M.resize(cols+6,cols+6);
    M.block(0,0,cols,cols)=tmp.block(0,0,cols,cols);
    //update UTMU other part
    Matd MIR=basis._M*_UDense.block(0,cols,rows,6);
    M.block(0,cols,cols+6,6)=_UDense.transpose()*MIR;
    M.block(cols,0,6,cols)=M.block(0,cols,cols,6).transpose();
    //Matd Mtmp=_UDense.transpose()*(basis._M*_UDense).eval();
}
RigidReducedMCalculator::Vec RigidReducedMCalculator::calculateFE(const Vec& S,Matd& K)
{
    TRIPS HTrips;
    Vec fFull=Vec::Zero(_body.nrV()*3);
    _energys.nonMaterialEnergyEval(&fFull,&HTrips,1.0f,1.0f,false);
    Eigen::SparseMatrix<scalarD,0,sizeType> HFull(_body.nrV()*3,_body.nrV()*3);
    HFull.setFromTriplets(HTrips.begin(),HTrips.end());

    K=_UDense.transpose()*(HFull*_UDense).eval();
    return _UDense.transpose()*fFull;
}
void RigidReducedMCalculator::getPos(Vec& P) const
{
    _innerMCalc->getPos(P);
    Vec tmp=P;
    P.resize(size());

    sizeType szI=_innerMCalc->size();
    P.block(0,0,szI,1)=tmp.block(0,0,szI,1);
    P.block(szI,0,3,1)=_body._X0.cast<scalarD>();
    P.block(szI+3,0,3,1)=_body._W.cast<scalarD>();
}
void RigidReducedMCalculator::getPosL(Vec& L) const
{
    _body.getPos(L);
}
RigidReducedMCalculator::Vec RigidReducedMCalculator::LtoF(const Vec& L) const
{
    return L;
}
void RigidReducedMCalculator::setPos(const Vec& P)
{
    sizeType szI=_innerMCalc->size();
    _innerMCalc->setPos(P.block(0,0,szI,1));
    _body._X0=P.block(szI,0,3,1).cast<scalar>();
    _body._W=P.block(szI+3,0,3,1).cast<scalar>();
    Mat3 R=expWGrad<Mat3d>(_body._W.cast<scalarD>()).cast<scalar>();

    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<_body.nrV(); i++) {
        FEMVertex& v=_body.getV(i);
        v._pos=R*v._pos+_body._X0;
    }
}
sizeType RigidReducedMCalculator::size() const
{
    return _innerMCalc->size()+6;
}
sizeType RigidReducedMCalculator::sizeL() const
{
    return _body.nrV()*3;
}