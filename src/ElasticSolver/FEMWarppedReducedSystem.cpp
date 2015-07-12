#include "FEMWarppedReducedSystem.h"
#include "FEMRotationUtil.h"
#include "FEMUtils.h"
#include "FEMMesh.h"

USE_PRJ_NAMESPACE

WrappedReducedMCalculator::WrappedReducedMCalculator(FEMReducedSystem& sys)
    :ReducedMCalculator(sys) {}
void WrappedReducedMCalculator::precompute()
{
    const SparseReducedBasis& basis=*(_body._basis);
    sizeType nrV=_body.nrV();
    sizeType nrC=_body.nrC();
    //build G
    TRIPS trips;
    Eigen::SparseMatrix<scalarD,0,sizeType> G;
    G.resize(nrC*9,nrV*3);
    _body.buildG(trips);
    G.setFromTriplets(trips.begin(),trips.end());

    //build W
    trips.clear();
    Eigen::SparseMatrix<scalarD,0,sizeType> W;
    W.resize(nrC*3,nrC*9);
    for(sizeType i=0; i<nrC; i++) {
#define ID(I,J) I+J*3
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+0,i*9+ID(2,1),+0.5f));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+0,i*9+ID(1,2),-0.5f));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+1,i*9+ID(0,2),+0.5f));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+1,i*9+ID(2,0),-0.5f));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+2,i*9+ID(1,0),+0.5f));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+2,i*9+ID(0,1),-0.5f));
#undef ID
    }
    W.setFromTriplets(trips.begin(),trips.end());

    //build Avg
    trips.clear();
    Eigen::SparseMatrix<scalarD,0,sizeType> A;
    A.resize(nrV*3,nrC*3);
    Vec nrNeigh(nrV);
    nrNeigh.setZero();
    for(sizeType i=0; i<nrC; i++) {
        const FEMCell& c=_body.getC(i);
        for(sizeType v=0; v<4; v++)
            if(c._v[v])nrNeigh[c._v[v]->_index]+=1.0f;
    }
    for(sizeType i=0; i<nrC; i++) {
        const FEMCell& c=_body.getC(i);
        for(sizeType v=0; v<4; v++)
            if(c._v[v])addI3x3(trips,c._v[v]->_index*3,i*3,1.0f/nrNeigh[c._v[v]->_index]);
    }
    A.setFromTriplets(trips.begin(),trips.end());
    _Gamma=A*W*G*basis._U;
}
void WrappedReducedMCalculator::calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U)
{
    //build rotation
    TRIPS trips;
    Vec W=_Gamma*S.block(0,0,_Gamma.cols(),1);
    _R.resize(_body.nrV()*3,_body.nrV()*3);
    for(sizeType i=0; i<_body.nrV(); i++)
        addBlock(trips,i*3,i*3,expWGrad<Mat3d>(W.block<3,1>(i*3,0)));
    _R.setFromTriplets(trips.begin(),trips.end());

    //build M and U
    const SparseReducedBasis& basis=*(_body._basis);
    M=basis._UTMU;
    _UDense=_R*basis._U;
    RHS.block(0,0,S.size(),1)+=basis._UTMU*MRHSL;
    if(U)U->addBlock(r,c,_UDense);
}
WrappedReducedMCalculator::Vec WrappedReducedMCalculator::calculateFE(const Vec& S,Matd& K)
{
    TRIPS HTrips;
    Vec fFull=Vec::Zero(_body.nrV()*3),f;
    _energys.nonMaterialEnergyEval(&fFull,&HTrips,1.0f,1.0f);
    Eigen::SparseMatrix<scalarD,0,sizeType> HFull(_body.nrV()*3,_body.nrV()*3);
    HFull.setFromTriplets(HTrips.begin(),HTrips.end());

    const SparseReducedBasis& basis=*(_body._basis);
    f=_UDense.transpose()*fFull;
    K=_UDense.transpose()*(HFull*_UDense).eval();
    return f;
}
void WrappedReducedMCalculator::setPos(const Vec& P)
{
    const SparseReducedBasis& basis=*(_body._basis);
    _body._S=P.block(0,0,size(),1);

    Vec W=_Gamma*_body._S;
    Vec dpos=basis._U*_body._S;
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<_body.nrV(); i++) {
        Eigen::Block<Vec,3,1> pi(dpos.block<3,1>(i*3,0));
        pi=expWGradInte(W.block<3,1>(i*3,0))*pi;
    }
    _body.setDPos(dpos);
}