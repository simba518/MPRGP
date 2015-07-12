#include "FEMRSReducedSystem.h"
#include "FEMRotationUtil.h"
#include "FEMCubatureSolver.h"
#include "FEMMesh.h"
#include "FEMUtils.h"
#include "PolyCoef.h"

USE_PRJ_NAMESPACE

//still we ignore centrifugal and coriolis force
CoupledRSReducedMCalculator::CoupledRSReducedMCalculator(FEMReducedSystem& sys)
    :ReducedMCalculator(sys),_constraints(sys.constraints()) {}
void CoupledRSReducedMCalculator::precompute()
{
    ReducedMCalculator::precompute();
    Eigen::SparseMatrix<scalarD,0,sizeType> PRJ;
    buildMatrix(PRJ,true);
    _invGTVG.compute(PRJ);
    if(_body._tree.get<bool>("debugRSDLDS",false))
        debugDLDS();
}
void CoupledRSReducedMCalculator::init() {}
void CoupledRSReducedMCalculator::calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U)
{
    if(!_body._tree.get<bool>("coupled")) {
        ReducedMCalculator::calculateM(S,M,MRHSL,RHS,r,c,U);
        return;
    }
    calculateMX(S);
    SparseReducedBasis& basis=*(_body._basis);
    RHS.block(0,0,_DLDS.cols(),1)+=_DLDS.transpose()*(basis._M*MRHSL);
    //calculate M and U
    M=_DLDS.transpose()*(basis._M*_DLDS).eval();
    if(U)U->addBlock(r,c,_DLDS);
}
CoupledRSReducedMCalculator::Vec CoupledRSReducedMCalculator::calculateFE(const Vec& S,Matd& K)
{
    if(_body._tree.get<bool>("coupled")) {
        TRIPS HTrips;
        Vec fFull=Vec::Zero(_body.nrV()*3);
        _energys.nonMaterialEnergyEval(&fFull,&HTrips,1.0f,1.0f,false);
        Eigen::SparseMatrix<scalarD,0,sizeType> HFull(_body.nrV()*3,_body.nrV()*3);
        HFull.setFromTriplets(HTrips.begin(),HTrips.end());

        K=_DLDS.transpose()*(HFull*_DLDS).eval();
        return _DLDS.transpose()*fFull;
    } else return ReducedMCalculator::calculateFE(S,K);
}
void CoupledRSReducedMCalculator::getPosL(Vec& L) const
{
    if(_body._tree.get<bool>("coupled"))
        L=_L;
    else ReducedMCalculator::getPosL(L);
}
CoupledRSReducedMCalculator::Vec CoupledRSReducedMCalculator::LtoF(const Vec& L) const
{
    if(_body._tree.get<bool>("coupled"))
        return L;
    else return ReducedMCalculator::LtoF(L);
}
void CoupledRSReducedMCalculator::setPos(const Vec& P)
{
    _body._S=P.block(0,0,size(),1);
    calculateMX(P);
}
sizeType CoupledRSReducedMCalculator::sizeL() const
{
    if(_body._tree.get<bool>("coupled"))
        return _body.nrV()*3;
    else return ReducedMCalculator::sizeL();
}
//private
void CoupledRSReducedMCalculator::setRSConstraint(Vec& X) const
{
    Vec tmp=X;
    X.setZero(_invGTVG.rows());
    X.block(0,0,tmp.rows(),1)=tmp;
}
void CoupledRSReducedMCalculator::calculateMX(const Vec& P)
{
    if(_lastS.size() == P.size() && (_lastS-P).isZero())
        return;
    else _lastS=P;

    //calculate DEXPDS
    Eigen::Matrix<scalarD,9,9> gradRSF;
    Vec FX=_GU*P.block(0,0,_GU.cols(),1);
    Matd DEXPDS=_GU;

    sizeType nrV=(sizeType)_body.nrV();
    OMP_PARALLEL_FOR_I(OMP_PRI(gradRSF))
    for(sizeType i=0; i<_body.nrC(); i++) {
        Eigen::Map<Mat3d> map(FX.data()+i*9);
        map=calcFGrad(map,&gradRSF)*_body.getC(i)._mass;
        DEXPDS.block(i*9,0,9,_GU.cols())=((gradRSF*_body.getC(i)._mass)*DEXPDS.block(i*9,0,9,_GU.cols())).eval();
    }
    //calculate UDense
    Vec X,tmp;
    setRSConstraint(X);
    _DLDS.setZero(nrV*3,size());
    for(sizeType i=0; i<_GU.cols(); i++) {
        X.block(0,0,nrV*3,1)=_G.transpose()*DEXPDS.col(i);
        tmp=_invGTVG.solve(X);
        _DLDS.col(i)=tmp.block(0,0,nrV*3,1);
    }
    //regenerate
    X.block(0,0,nrV*3,1)=_G.transpose()*FX;
    tmp=_invGTVG.solve(X);
    _L=tmp.block(0,0,nrV*3,1);
    _body.setDPos(_L);
}
void CoupledRSReducedMCalculator::buildMatrix(Eigen::SparseMatrix<scalarD,0,sizeType>& PRJ,bool cons)
{
    sizeType nrV=(sizeType)_body.nrV();
    sizeType nrC=(sizeType)_body.nrC();

    //build G
    TRIPS trips;
    _G.resize(nrC*9,nrV*3);
    _body.buildG(trips);
    _G.setFromTriplets(trips.begin(),trips.end());

    //build V
    trips.clear();
    _V.resize(_G.rows(),_G.rows());
    for(sizeType i=0; i<nrC; i++)
        addI(trips,i*9,i*9,Vec::Constant(9,_body.getC(i)._mass));
    _V.setFromTriplets(trips.begin(),trips.end());

    //build solver
    trips.clear();
    Eigen::SparseMatrix<scalarD,0,sizeType> GTVG=_G.transpose()*_V*_G;
    addSparseBlock(trips,0,0,GTVG);
    if(_body.dim() == 2)
        for(sizeType i=0; i<nrV; i++)
            trips.push_back(Eigen::Triplet<scalarD,sizeType>(i*3+2,i*3+2,1.0f));

    //build constraint
    sizeType szMat=nrV*3;
    if(cons) {
        if(_constraints.empty()) {
            for(sizeType i=0; i<nrV; i++) {
                const FEMVertex& v=_body.getV(i);
                addI3x3(trips,szMat,v._index*3,v._mass);
                addI3x3(trips,v._index*3,szMat,v._mass);
            }
            szMat+=3;
        } else {
            for(boost::unordered_map<sizeType,Vec3>::const_iterator beg=_constraints.begin(),end=_constraints.end(); beg!=end; beg++) {
                addI3x3(trips,szMat,beg->first*3,_body.getV(beg->first)._mass);
                addI3x3(trips,beg->first*3,szMat,_body.getV(beg->first)._mass);
                szMat+=3;
            }
        }
    }
    PRJ.resize(szMat,szMat);
    PRJ.setFromTriplets(trips.begin(), trips.end());

    //build GU
    const SparseReducedBasis& basis=*(_body._basis);
    _GU=_G*basis._U;
}
void CoupledRSReducedMCalculator::debugDLDS()
{
    for(sizeType i=0; i<10; i++) {
        Vec S=Vec::Random(size());
        calculateMX(S);
        Vec L0=_L;
        Matd DLDS0=_DLDS;

#define DELTA 1E-7f
        for(sizeType j=0; j<size(); j++) {
            Vec STmp=S+Vec::Unit(size(),j)*DELTA;
            calculateMX(STmp);
            Vec DLDSN=(_L-L0)/DELTA;
            INFOV("Error: %f %f",DLDS0.col(j).norm(),(DLDSN-DLDS0.col(j)).norm())
        }
#undef DELTA
    }
}

//now we use cubature to accelerate the computation
class KineticCubatureProb : public FEMCubatureProb
{
public:
    KineticCubatureProb(const FEMBody& body,const Matd& FtoB,const Matd& GU)
        :_body(body),_FtoB(FtoB),_GU(GU) {}
    virtual Matd operator()(const sizeType cid,const Matd& S) const {
        //build the interpolation matrix
        FEMCell& cell=*(_body.getCPtr(cid));
        Eigen::Matrix<scalarD,9,1> F;
        Matd GUBLK=_GU.block(cid*9,0,9,S.rows());
        Matd FtoBBLK=_FtoB.block(0,cid*9,_FtoB.rows(),9);
        Matd ret(_FtoB.rows()*S.cols(),nr()[1]);
        for(sizeType pose=0; pose<S.cols(); pose++) {
            F=GUBLK*S.col(pose);
            Eigen::Map<Mat3d> mapF(F.data());
            mapF=calcFGrad(mapF,NULL);
            ret.block(pose*_FtoB.rows(),0,_FtoB.rows(),1)=FtoBBLK*F;
        }
        return ret;
    }
    virtual void debugCallback(const Vec& b,const Matd& S) const {
        if(!_body._tree.get<bool>("debugKineticCubature",false))
            return;
        for(sizeType i=0; i<10; i++) {
            sizeType c=rand()%S.cols();
            Vec FX=_GU*S.col(c);
            sizeType nrV=(sizeType)_body.nrV();
            OMP_PARALLEL_FOR_
            for(sizeType i=0; i<_body.nrC(); i++) {
                Eigen::Map<Mat3d> map(FX.data()+i*9);
                map=calcFGrad(map,NULL)*_body.getC(i)._mass;
            }
            Vec BLKb=b.block(c*nr()[0],0,nr()[0],1);
            INFOV("DebugKineticCubature: %f %f",BLKb.norm(),(_FtoB*FX-BLKb).norm());
        }
    }
    virtual Vec2i nr() const {
        return Vec2i(_FtoB.rows(),1);
    }
protected:
    const FEMBody& _body;
    const Matd &_FtoB,&_GU;
};
CRCoupledRSReducedMCalculatorBase::CRCoupledRSReducedMCalculatorBase(FEMReducedSystem& sys)
    :CoupledRSReducedMCalculator(sys)
{
    if(_body._tree.find("basisSizeLimit") == _body._tree.not_found())
        _body._tree.put<sizeType>("basisSizeLimit",90);
    if(_body._tree.find("overwriteRS") == _body._tree.not_found())
        _body._tree.put<bool>("overwriteRS",false);
}
void CRCoupledRSReducedMCalculatorBase::buildBasis()
{
    //reset BASIS to LMA basis and initialize UTKU and UTMU
    Matd& BASIS=_body._basis->_B;
    if(BASIS.rows() == 0 || BASIS.cols() == 0)
        BASIS=_body._basis->_U;
    sizeType nrLMA=_body._basis->nrLMA(BASIS);
    //this means user has already provided their desired reduced space, we don't replace these
    if(!_body._tree.get<bool>("overwriteRS") && BASIS.cols() > nrLMA)
        return;

    //how many basis do we use
    sizeType colLMA=BASIS.cols();
    sizeType nrDV=colLMA*(colLMA+1)/2;
    sizeType nr=std::min(_body._tree.get<sizeType>("nrDU",0),_body._tree.get<sizeType>("basisSizeLimit")-nrLMA);
    if(nr > 0)
        nrDV=std::min(nrDV,nr);
    Matd deriv(BASIS.rows(),nrDV+colLMA);

    //precompute invGTVG
    Eigen::SparseMatrix<scalarD,0,sizeType> PRJ;
    buildMatrix(PRJ,true);
    _invGTVG.compute(PRJ);

    //filter model basis
    Vec X,tmp;
    setRSConstraint(X);
    sizeType nrC=_body.nrC();
    sizeType k=0;
    for(sizeType i=0; i<colLMA; i++,k++) {
        Vec FX=_G*BASIS.col(i);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<_body.nrC(); i++) {
            Eigen::Map<Mat3d> map(FX.data()+i*9);
            map=calcFGradDU(map)*_body.getC(i)._mass;
        }
        X.block(0,0,_G.cols(),1)=_G.transpose()*FX;
        tmp=_invGTVG.solve(X);
        deriv.col(k)=tmp.block(0,0,BASIS.rows(),1);
        if(_body._tree.get<bool>("debugRSB",false))
            debugRSB1(deriv.col(k),i);
    }

    //filter model derivative
    for(sizeType i=0; i<colLMA; i++)
        for(sizeType j=i; j<colLMA; j++,k++) {
            Vec FX=_G*BASIS.col(i);
            Vec FXV=_G*BASIS.col(j);
            OMP_PARALLEL_FOR_
            for(sizeType i=0; i<_body.nrC(); i++) {
                Eigen::Map<Mat3d> map(FX.data()+i*9);
                Eigen::Map<Mat3d> mapV(FXV.data()+i*9);
                if(i == j)
                    map=calcFGradDDU(map)*_body.getC(i)._mass;
                else map=calcFGradDUDV(map,mapV)*_body.getC(i)._mass;
            }
            X.block(0,0,_G.cols(),1)=_G.transpose()*FX;
            tmp=_invGTVG.solve(X);
            deriv.col(k)=tmp.block(0,0,BASIS.rows(),1);

            if(_body._tree.get<bool>("debugRSB",false))
                debugRSB2(deriv.col(k),i,j);
        }

    //normalize
    ASSERT(k == deriv.cols())
    for(sizeType i=0; i<k; i++)
        deriv.col(i).normalize();
    BASIS=deriv;
    boost::dynamic_pointer_cast<FEMReducedSystem>(_body._system)->writeBasisVTK("./basisDerivExp/",1.0f,false);
}
Matd CRCoupledRSReducedMCalculatorBase::calculateFtoB()
{
    Eigen::SparseMatrix<scalarD,0,sizeType> PRJ;
    buildMatrix(PRJ,false);

    //calculate invPBT using [Stanton2013]
    Matd& B=_body._basis->_B;
    const Matd& U0=_body._basis->_U;
    Matd P=B.transpose()*(PRJ*B).eval();
    Matd invP=P.inverse();

    Matd FtoB;
    //using cubature to further accelerate
    {
        Vec weight(_body.nrC());
        for(sizeType i=0; i<weight.size(); i++)
            weight[i]=_body.getC(i)._mass;

        _GU=_G*U0;
        FtoB=invP*(_G*B).eval().transpose();
        Eigen::SparseMatrix<scalarD,0,sizeType> mass;
        _body.getMass(mass,true);

        KineticCubatureProb::transformMetric(mass,B,FtoB);
        KineticCubatureProb obj(_body,FtoB,_GU);

        _body._tree.put<std::string>("pathCub","meshCubatureKinetic");
        FEMCubatureSolver sol(_body,obj);
        sol.generate(weight);
        TBEG("Compute Cubature")
        sol.solve(_body._cubatureKinetic,".CUBK");
        TEND
    }

    //get cubature cells
    _body.clearPos();
    const Cubature& cub=*(_body._cubatureKinetic);
    CFEMBody cBody(_body,cub);

    //build matrix
    sizeType nrV=(sizeType)cBody.nrV();
    sizeType nrC=(sizeType)cBody.nrC();

    TRIPS trips;
    cBody.buildG(trips);
    _G.resize(nrC*9,nrV*3);
    _G.setFromTriplets(trips.begin(),trips.end());
    _GU=_G*(cBody._cInterp*U0).eval();

    trips.clear();
    _V.resize(nrC*9,nrC*9);
    for(sizeType i=0; i<nrC; i++)
        addI(trips,i*9,i*9,Vec::Constant(9,cub._weight[i]));
    _V.setFromTriplets(trips.begin(),trips.end());

    FtoB=invP*B.transpose();
    FtoB=FtoB*cBody._cInterp.transpose();
    FtoB=FtoB*_G.transpose()*_V;
    return FtoB;
}
void CRCoupledRSReducedMCalculatorBase::debugRSB1(const Vec& col,sizeType i)
{
    Vec in=Vec::Unit(size(),i);
#define DELTA 1E-7f
    CoupledRSReducedMCalculator::calculateMX(in*DELTA);
    Vec dx=_DLDS.col(i);
    for(sizeType i=0; i<dx.size(); i++)
        INFOV("DebugRSB1: %f %f",col[i],std::abs(dx[i]-col[i]))
        //system("pause");
#undef DELTA
    }
void CRCoupledRSReducedMCalculatorBase::debugRSB2(const Vec& col,sizeType i,sizeType j)
{
    Vec inI=Vec::Unit(size(),i);
    Vec inJ=Vec::Unit(size(),j);
#define DELTA 1E-7f
    CoupledRSReducedMCalculator::calculateMX(inI*DELTA);
    Vec L1=_DLDS.col(j);
    CoupledRSReducedMCalculator::calculateMX(Vec::Zero(size()));
    Vec L0=_DLDS.col(j);
    Vec dx=(L1-L0)/DELTA;
    for(sizeType i=0; i<dx.size(); i++)
        INFOV("DebugRSB2: %f %f",col[i],std::abs(dx[i]-col[i]))
        //system("pause");
#undef DELTA
    }

//basic version for debug
struct ACCoupledBasisExtractor : public BasisExtractor {
    ACCoupledBasisExtractor(const Matd& DLTbDS,const Matd& BInvLT,const Cold& b)
        :_DLTbDS(DLTbDS),_BInvLT(BInvLT),_b(b) {}
    virtual Matd operator()(sizeType i) const {
        return _BInvLT.block(i*3,0,3,_BInvLT.cols())*_DLTbDS;
    }
    virtual Vec3d operator[](sizeType i) const {
        return _BInvLT.block(i*3,0,3,_BInvLT.cols())*_b;
    }
    const Matd &_DLTbDS,&_BInvLT;
    const Cold &_b;
};
CRCoupledRSReducedMCalculatorRef::CRCoupledRSReducedMCalculatorRef(FEMReducedSystem& sys)
    :CRCoupledRSReducedMCalculatorBase(sys)
{
    if(_body._tree.find("overwriteCub") == _body._tree.not_found())
        _body._tree.put<bool>("overwriteCub",false);
}
void CRCoupledRSReducedMCalculatorRef::precompute()
{
    ReducedMCalculator::precompute();
    Matd FtoB=calculateFtoB();
    const Matd& B=_body._basis->_B;
    _BTMB.reset(new Matd);
    *_BTMB=B.transpose()*(_body._basis->_M*B).eval();
    _BCub=FtoB;

    //precompute external force
    if(_body._tree.get<bool>("coupled")) {
        _UTf.setZero(sizeL());
        _UTHU.setZero(sizeL(),sizeL());
        _energys.lockNonMaterial(false);
        _energys.nonMaterialEnergyEval(MatrixBasisExtractor(getBasis(),Vec::Zero(sizeL())),&_UTf,&_UTHU,1.0f,1.0f);
        _energys.lockNonMaterial(true);
    }
    if(_body._tree.get<bool>("debugRSDLDS",false))
        debugDLDS();
}
void CRCoupledRSReducedMCalculatorRef::calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U)
{
    if(!_body._tree.get<bool>("coupled")) {
        ReducedMCalculator::calculateM(S,M,MRHSL,RHS,r,c,U);
        return;
    }
    const SparseReducedBasis& basis=*(_body._basis);
    calculateMX(S);
    M=_matM;
    if(_BTMB)
        RHS.block(0,0,size(),1)+=_DLDS.transpose()*(*_BTMB*MRHSL).eval();
    else RHS.block(0,0,size(),1)+=_DLDS.transpose()*MRHSL;
    if(U)U->addBlock(r,c,getBasis()*_DLDS);
}
CRCoupledRSReducedMCalculatorRef::Vec CRCoupledRSReducedMCalculatorRef::calculateFE(const Vec& S,Matd& K)
{
    if(_body._tree.get<bool>("coupled")) {
        Vec f=Vec::Zero(size());
        K.setZero(size(),size());

        f.block(0,0,_DLDS.cols(),1)=_DLDS.transpose()*_UTf;
        K.block(0,0,_DLDS.cols(),_DLDS.cols())=_DLDS.transpose()*(_UTHU*_DLDS).eval();
        _energys.nonMaterialEnergyEval(ACCoupledBasisExtractor(_DLDS,getBasis(),_L),&f,&K,1.0f,1.0f);
        return f;
    } else return ReducedMCalculator::calculateFE(S,K);
}
CRCoupledRSReducedMCalculatorRef::Vec CRCoupledRSReducedMCalculatorRef::LtoF(const Vec& L) const
{
    if(_body._tree.get<bool>("coupled"))
        return getBasis()*L;
    else return ReducedMCalculator::LtoF(L);
}
sizeType CRCoupledRSReducedMCalculatorRef::sizeL() const
{
    if(_body._tree.get<bool>("coupled"))
        return getBasis().cols();
    else return ReducedMCalculator::sizeL();
}
void CRCoupledRSReducedMCalculatorRef::calculateMX(const Vec& P)
{
    if(_lastS.size() == P.size() && (_lastS-P).isZero())
        return;
    else _lastS=P;

    sizeType nrB=getBasis().cols();
    _DLDS.setZero(nrB,size());
    _L=Vec::Zero(nrB);
    Vec FX=_GU*P.block(0,0,size(),1);
    Eigen::Matrix<scalarD,9,9> gradRSF;
    for(sizeType i=0; i<_GU.rows()/9; i++) {
        Eigen::Map<Mat3d> map(FX.data()+i*9);
        map=calcFGrad(map,&gradRSF);
        Eigen::Block<Matd> BLK=_BCub.block(0,i*9,nrB,9);
        _L+=BLK*Eigen::Map<Eigen::Matrix<scalarD,9,1> >(FX.data()+i*9);
        _DLDS+=BLK*(gradRSF*_GU.block(i*9,0,9,size())).eval();
    }
    if(_BTMB)
        _matM=_DLDS.transpose()*(*_BTMB*_DLDS).eval();
    else _matM=_DLDS.transpose()*_DLDS;
    if(_body._tree.get<bool>("reconstruct"))
        _body.setDPos(getBasis()*_L);
}
const Matd& CRCoupledRSReducedMCalculatorRef::getBasis() const
{
    return _body._basis->_B;
}

//multiplication with BTMB factored out
CRCoupledRSReducedMCalculator::CRCoupledRSReducedMCalculator(FEMReducedSystem& sys)
    :CRCoupledRSReducedMCalculatorRef(sys) {}
void CRCoupledRSReducedMCalculator::precompute()
{
    ReducedMCalculator::precompute();
    Matd FtoB=calculateFtoB();
    const Matd& B=_body._basis->_B;

    _BTMB.reset((Matd*)NULL);
    Matd BTMB=B.transpose()*(_body._basis->_M*B).eval();
    Matd LTBTMB=BTMB.llt().matrixL().transpose();
    _BCub=LTBTMB*FtoB;
    _BInvLT=B*LTBTMB.inverse();

    //precompute external force
    if(_body._tree.get<bool>("coupled")) {
        _UTf.setZero(sizeL());
        _UTHU.setZero(sizeL(),sizeL());
        _energys.lockNonMaterial(false);
        _energys.nonMaterialEnergyEval(MatrixBasisExtractor(getBasis(),Vec::Zero(sizeL())),&_UTf,&_UTHU,1.0f,1.0f);
        _energys.lockNonMaterial(true);
    }
    if(_body._tree.get<bool>("debugRSDLDS",false))
        debugDLDS();
}
const Matd& CRCoupledRSReducedMCalculator::getBasis() const
{
    return _BInvLT;
}

//we further accelerate it by taylor series
struct PadeStruct {
    void buildCoef(sizeType p,sizeType q) {
        _pq=Vec2i(p,q);
        _PCoef.resize(p+1);
        _QCoef.resize(q+1);
        std::cout << "Pade P: ";
        for(sizeType i=0; i<=p; i++) {
            _PCoef[i]=Ppq(p,q,i);
            std::cout << _PCoef[i] << "*W^" << i;
            if(i < p)std::cout << "+";
        }
        std::cout << std::endl;

        std::cout << "Pade Q: ";
        for(sizeType i=0; i<=q; i++) {
            _QCoef[i]=Qpq(p,q,i);
            std::cout << _QCoef[i] << "*W^" << i;
            if(i < q)std::cout << "+";
        }
        std::cout << std::endl;
    }
    static scalarD fac(sizeType j) {
        scalarD init=1.0f;
        for(sizeType i=2; i<=j; i++)
            init*=(scalarD)i;
        return init;
    }
    static scalarD Ppq(sizeType p,sizeType q,sizeType j) {
        scalarD n=fac(p+q-j)*fac(p);
        scalarD d=fac(p+q)*fac(j)*fac(p-j);
        return n/d;
    }
    static scalarD Qpq(sizeType p,sizeType q,sizeType j) {
        scalarD n=fac(p+q-j)*fac(q);
        scalarD d=fac(p+q)*fac(j)*fac(q-j);
        return (n/d)*(j%2==1?-1.0f:1.0f);
    }
    //Pade coef
    Vec2i _pq;
    Cold _PCoef,_QCoef;
    //other data
    Matd _GB,_GU;
    Matd _coefP,_coefQ;
    //temporary data
    Matd _VQGB;
};
struct PadePRSPolynomial : public PolyCoefEfficient {
    PadePRSPolynomial(Matd& coef,boost::shared_ptr<PadeStruct> info)
        :PolyCoefEfficient(coef,info->_pq[0]+1),_info(info),
         _nrB(info->_GB.cols()),_nrS(info->_GU.cols()) {}
    Cold findF(const Cold& S) const {
        const Matd& GU=_info->_GU;
        Cold SBlk=S.block(0,0,_nrS,1);
        Cold FX=GU*SBlk;
        for(sizeType i=0; i<GU.rows()/9; i++) {
            Eigen::Map<Mat3d> map(FX.data()+i*9);
            map=calcFGradPadeP(map);
        }
        return _info->_coefP*FX;
    }
    Mat3d calcFGradPadeP(const Mat3d& F) const {
        Mat3d W=(F-F.transpose())/2.0f;
        Mat3d S=(F+F.transpose())/2.0f;
        Mat3d EXPP=Mat3d::Zero();
        Mat3d EXPQ=Mat3d::Zero();
        Mat3d tmp=Mat3d::Identity();
        sizeType maxPQ=_info->_pq.maxCoeff();
        for(sizeType ord=0; ord<=maxPQ; ord++) {
            if(ord <= _info->_pq[0])
                EXPP+=tmp*_info->_PCoef[ord];
            if(ord <= _info->_pq[1])
                EXPQ+=tmp*_info->_QCoef[ord];
            tmp*=W;
        }
        return EXPP*(S+Mat3d::Identity())-EXPQ;
    }
    sizeType nrF() const {
        return _nrB;
    }
    sizeType nrX() const {
        return _nrS;
    }
    boost::shared_ptr<PadeStruct> _info;
    sizeType _nrB,_nrS;
};
struct PadeQRSPolynomial : public PolyCoefEfficient {
    PadeQRSPolynomial(Matd& coef,boost::shared_ptr<PadeStruct> info)
        :PolyCoefEfficient(coef,info->_pq[1]),_info(info),
         _nrB(info->_GB.cols()),_nrS(info->_GU.cols()) {
        info->_VQGB=info->_GB;
    }
    Cold findF(const Cold& S) const {
        const Matd& GU=_info->_GU;
        const Matd& coefQ=_info->_coefQ;
        Matd& VQGB=_info->_VQGB;
        Cold SBlk=S.block(0,0,_nrS,1);
        Cold FX=GU*SBlk;
        for(sizeType i=0; i<GU.rows()/9; i++) {
            Eigen::Map<Mat3d> map(FX.data()+i*9);
            map=calcFGradPadeQ(map);
            for(char c=0; c<9; c+=3)
                VQGB.block(i*9+c,0,3,_nrB)=map*coefQ.block(i*9+c,0,3,_nrB);
        }
        Matd out=_info->_GB.transpose()*VQGB;
        return Eigen::Map<Cold>(out.data(),out.size());
    }
    Mat3d calcFGradPadeQ(const Mat3d& F) const {
        Mat3d W=(F-F.transpose())/2.0f;
        Mat3d EXPQ=Mat3d::Zero();
        Mat3d tmp=Mat3d::Identity();
        for(sizeType ord=0; ord<=_info->_pq[1]; ord++) {
            EXPQ+=tmp*_info->_QCoef[ord];
            tmp*=W;
        }
        return EXPQ;
    }
    sizeType nrF() const {
        return _nrB*_nrB;
    }
    sizeType nrX() const {
        return _nrS;
    }
    boost::shared_ptr<PadeStruct> _info;
    sizeType _nrB,_nrS;
};
TRCoupledRSReducedMCalculator::TRCoupledRSReducedMCalculator(FEMReducedSystem& sys)
    :CRCoupledRSReducedMCalculator(sys)
{
    if(_body._tree.find("overwriteApproxRS") == _body._tree.not_found())
        _body._tree.put<bool>("overwriteApproxRS",false);
    if(_body._tree.find("PadeP") == _body._tree.not_found())
        _body._tree.put<sizeType>("PadeP",3);
    if(_body._tree.find("PadeQ") == _body._tree.not_found())
        _body._tree.put<sizeType>("PadeQ",0);
}
void TRCoupledRSReducedMCalculator::precompute()
{
    boost::shared_ptr<PadeStruct> info(new PadeStruct);
    info->buildCoef(_body._tree.get<sizeType>("PadeP"),_body._tree.get<sizeType>("PadeQ"));
    {
        sizeType nrV=(sizeType)_body.nrV();
        sizeType nrC=(sizeType)_body.nrC();

        TRIPS trips;
        _body.buildG(trips);
        _G.resize(nrC*9,nrV*3);
        _G.setFromTriplets(trips.begin(),trips.end());
        _GU=_G*_body._basis->_U;

        trips.clear();
        _V.resize(nrC*9,nrC*9);
        for(sizeType i=0; i<nrC; i++)
            addI(trips,i*9,i*9,Vec::Constant(9,_body.getC(i)._mass));
        _V.setFromTriplets(trips.begin(),trips.end());
    }
    info->_GB=_G*_body._basis->_B;
    info->_GU=_GU;
    info->_coefP=info->_GB.transpose()*_V;

    _BTMB.reset((Matd*)NULL);
    Matd BTMB=_body._basis->_B.transpose()*(_body._basis->_M*_body._basis->_B).eval();
    Matd LTBTMB=BTMB.llt().matrixL().transpose();
    info->_coefQ=_V*info->_GB*LTBTMB.inverse().eval();
    _coefP.reset(new PadePRSPolynomial(_body._basis->_RSP,info));
    _coefQ.reset(new PadeQRSPolynomial(_body._basis->_RSQ,info));
    if(_body._tree.get<sizeType>("PadeQ") > 0) {
        _coefP->findCoef(_body._tree.get<bool>("overwriteApproxRS"));
        _coefQ->findCoef(_body._tree.get<bool>("overwriteApproxRS"));
        boost::dynamic_pointer_cast<PadePRSPolynomial>(_coefP)->_info.reset((PadeStruct*)NULL);
        boost::dynamic_pointer_cast<PadeQRSPolynomial>(_coefQ)->_info.reset((PadeStruct*)NULL);
    } else {
        Cold Q=_coefQ->findF(Cold::Zero(info->_GU.cols()));
        Matd matQ=Eigen::Map<Matd>(Q.data(),info->_GB.cols(),info->_GB.cols());
        info->_coefP=matQ.inverse()*info->_coefP;
        _coefP->findCoef(_body._tree.get<bool>("overwriteApproxRS"));
        boost::dynamic_pointer_cast<PadePRSPolynomial>(_coefP)->_info.reset((PadeStruct*)NULL);
        _coefQ.reset((PadeQRSPolynomial*)NULL);
    }
    _BInvLT=_body._basis->_B*LTBTMB.inverse().eval();

    //precompute external force
    if(_body._tree.get<bool>("coupled")) {
        _UTf.setZero(sizeL());
        _UTHU.setZero(sizeL(),sizeL());
        _energys.lockNonMaterial(false);
        _energys.nonMaterialEnergyEval(MatrixBasisExtractor(getBasis(),Vec::Zero(sizeL())),&_UTf,&_UTHU,1.0f,1.0f);
        _energys.lockNonMaterial(true);
    }
    if(_body._tree.get<bool>("debugRSDLDS",false))
        debugDLDS();
}
void TRCoupledRSReducedMCalculator::calculateMX(const Vec& P)
{
    if(_lastS.size() == P.size() && (_lastS-P).isZero())
        return;
    else _lastS=P;

    sizeType nrB=_coefP->nrF();
    sizeType nrS=_coefP->nrX();
    if(!_coefQ) {
        _L=_coefP->findFK(P,_DLDS);
    } else {
        //|b||s|^4
        Matd DPDS;
        Vec p=_coefP->findFK(P,DPDS);
        //|b|^2|s|^2
        Matd DQDS;
        Cold Q=_coefQ->findFK(P,DQDS);
        //|b|^3
        Matd matQ=Eigen::Map<Matd>(Q.data(),nrB,nrB);
        Eigen::FullPivLU<Matd> luQ=matQ.fullPivLu();
        _L=luQ.solve(p);
        //|b|^2|s|
        for(sizeType i=0; i<nrS; i++)
            DPDS.col(i)-=Eigen::Map<Matd>(DQDS.col(i).data(),nrB,nrB)*_L;
        //|b|^2|s|
        _DLDS=luQ.solve(DPDS);
    }
    _matM=_DLDS.transpose()*_DLDS;
    if(_body._tree.get<bool>("reconstruct"))
        _body.setDPos(getBasis()*_L);
}