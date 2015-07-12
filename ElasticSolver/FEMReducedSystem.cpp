#include "FEMReducedSystem.h"
#include "FEMSparseReducedBasis.h"
#include "FEMCubatureSolver.h"
#include "FEMUtils.h"
#include "PolyCoef.h"
#include "FEMMesh.h"
#include "solvers/PMinresQLP.h"
#include <boost/unordered_set.hpp>
#include <boost/filesystem/operations.hpp>

//load FEM models
#include "FEMRigidReducedSystem.h"
#include "FEMWarppedReducedSystem.h"
#include "FEMRSReducedSystem.h"

USE_PRJ_NAMESPACE

//basic M calculator
ReducedMCalculator::ReducedMCalculator(FEMReducedSystem& sys)
    :_body(sys.body()),_energys(sys.energys()) {}
void ReducedMCalculator::buildBasis() {}
void ReducedMCalculator::precompute()
{
    const Matd& U=_body._basis->_U;
    sizeType nr=U.cols();
    _UTf.setZero(nr);
    _UTHU.setZero(nr,nr);

    _energys.lockNonMaterial(false);
    _energys.nonMaterialEnergyEval(MatrixBasisExtractor(U,Vec::Zero(nr)),&_UTf,&_UTHU,1.0f,1.0f);
    _energys.lockNonMaterial(true);
}
void ReducedMCalculator::init() {}
void ReducedMCalculator::calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U)
{
    const SparseReducedBasis& basis=*(_body._basis);
    M=basis._UTMU;
    RHS.block(0,0,basis._U.cols(),1)+=basis._UTMU*MRHSL;
    if(U)U->addDense(r,c,basis._U);
}
ReducedMCalculator::Vec ReducedMCalculator::calculateFE(const Vec& S,Matd& K)
{
    Vec f=Vec::Zero(S.size());
    K.setZero(S.size(),S.size());
    f.block(0,0,_UTf.size(),1)=_UTf;
    K.block(0,0,_UTHU.rows(),_UTHU.cols())=_UTHU;
    _energys.nonMaterialEnergyEval(MatrixBasisExtractor(_body._basis->_U,S),&f,&K,1.0f,1.0f);
    return f;
}
void ReducedMCalculator::getPos(Vec& P) const
{
    P=_body._S;
}
void ReducedMCalculator::getPosL(Vec& P) const
{
    return getPos(P);
}
ReducedMCalculator::Vec ReducedMCalculator::LtoF(const Vec& L) const
{
    SparseReducedBasis& basis=*(_body._basis);
    return basis._U*L.block(0,0,basis._U.cols(),1);
}
void ReducedMCalculator::setPos(const Vec& P)
{
    SparseReducedBasis& basis=*(_body._basis);
    _body._S=P.block(0,0,size(),1);
    if(_body._tree.get<bool>("reconstruct"))
        _body.setDPos(LtoF(_body._S));
}
sizeType ReducedMCalculator::size() const
{
    return _body._basis->_U.cols();
}
sizeType ReducedMCalculator::sizeL() const
{
    return size();
}

PRJ_BEGIN
//basic K calculator
class ReducedKCalculator
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    ReducedKCalculator(FEMReducedSystem& sys)
        :_body(sys.body()),_energys(sys.energys()) {}
    virtual void buildBasis() {}
    virtual void precompute() {}
    virtual void init() {}
    virtual Vec calculateFI(const Vec& S,Matd& K) const {
        const SparseReducedBasis& basis=*(_body._basis);
        K=basis._UTKU;
        return -K*S.block(0,0,K.cols(),1);
    }
    void debugFI() const {
        sizeType nrX=_body._S.size();
        for(sizeType i=0; i<10; i++) {
            Vec testS(nrX),testS2;
            testS.setRandom();

            Matd KSub,KSub2;
            Vec fSub=calculateFI(testS,KSub),fSub2;

            for(sizeType j=0; j<nrX; j++) {
#define DELTA 1E-7f
                testS2=testS+Vec::Unit(nrX,j)*DELTA;
                fSub2=calculateFI(testS2,KSub2);
                Vec DN=(fSub2-fSub)/DELTA;
                INFOV("Jacobian: %f %f",KSub.col(j).norm(),(DN+KSub.col(j)).norm());
                ASSERT((DN+KSub.col(j)).norm() < 1E-5f*KSub.col(j).norm())
#undef DELTA
            }
        }
    }
protected:
    FEMBody& _body;
    const EnergyPool& _energys;
};
//STVK K calculator
class ReducedKCalculatorSTVK : public ReducedKCalculator, protected PolyCoefEfficient
{
public:
    ReducedKCalculatorSTVK(FEMReducedSystem& sys)
        :ReducedKCalculator(sys),PolyCoefEfficient(_body._basis->_STVKU,3),_constraints(sys.constraints()) {
        if(_body._tree.find("overwriteSTVK") == _body._tree.not_found())
            _body._tree.put<bool>("overwriteSTVK",false);
        if(_body._tree.find("basisSizeLimit") == _body._tree.not_found())
            _body._tree.put<sizeType>("basisSizeLimit",90);
    }
    virtual void buildBasis() {
        ASSERT_MSGV(_energys.getMaterialType() == MaterialEnergy::STVK,
                    "Model derivative can only be calculated for STVK model!")
        //reset BASIS to LMA basis and initialize UTKU and UTMU
        Matd& BASIS=_body._basis->_U;
        sizeType nrLMA=_body._basis->nrLMA(BASIS);
        //this means user has already provided their desired reduced space, we don't replace these
        if(!_body._tree.get<bool>("overwriteSTVK") && BASIS.cols() > nrLMA)
            return;
        BASIS=BASIS.block(0,0,BASIS.rows(),nrLMA).eval();
        Vec UTKUD=(BASIS.transpose()*(_body._basis->_K*BASIS).eval()).diagonal();
        Vec UTMUD=(BASIS.transpose()*(_body._basis->_M*BASIS).eval()).diagonal();

        //setup solver and parameter
        boost::shared_ptr<Solver<scalarD,Kernel<scalarD> > > sol;
        if(_body._tree.get<bool>("useDirectFactorize",true)) {
            sol.reset(new PMINRESSolverQLP<scalarD,Kernel<scalarD>,DirectPreconSolver<scalarD> >);
            DirectPreconSolver<scalarD>* pre=(DirectPreconSolver<scalarD>*)(sol->getPre());
            Eigen::SparseMatrix<scalarD,0,sizeType> KPre=_body._basis->_K;
            EigenSolver::scaleDiagonal(KPre,1.01f,1E-2f);
            pre->setMatrix(KPre,true);
        } else sol.reset(new PMINRESSolverQLP<scalarD,Kernel<scalarD> >);

        //calculate reduced basis
        //for each pair of <i,j> we solve:
        //(K-\lambda_iM)U_{ij}=(MU_iU_i^T-I)\FPP{K}{z_j}U_i
        sizeType colLMA=BASIS.cols();
        sizeType nrDV;
        if(_constraints.empty())
            nrDV=colLMA*colLMA;
        else nrDV=colLMA*(colLMA+1)/2;
        sizeType nr=std::min(_body._tree.get<sizeType>("nrDU",0),_body._tree.get<sizeType>("basisSizeLimit")-nrLMA);
        if(nr > 0)
            nrDV=std::min(nrDV,nr);
        Vec RHS(BASIS.rows()),tmp(BASIS.rows());
        Matd deriv(BASIS.rows(),nrDV+colLMA);
        deriv.block(0,0,BASIS.rows(),BASIS.cols())=BASIS;
        for(sizeType i=0,k=colLMA; k<deriv.cols(); i++) {
            FixedSparseMatrix<scalarD,Kernel<scalarD> > LHS;
            scalarD lambda=UTKUD[i]/UTMUD[i];
            if(!_constraints.empty())
                lambda=0.0f;
            LHS.fromEigen((_body._basis->_K-lambda*_body._basis->_M).eval());

            sol->setSolverParameters(1E-7f,10000);
            sol->setMatrix(LHS,false);
            for(sizeType j=_constraints.empty()?0:i; j<colLMA; j++) {
                //to calculate this
                Vec DKDUjUi=calculateDKDUjUi(BASIS.col(j),BASIS.col(i));
                RHS=-DKDUjUi;
                if(_constraints.empty())
                    RHS+=_body._basis->_M*BASIS.col(i)*DKDUjUi.dot(BASIS.col(i));
                ASSERT(sol->solve(RHS,tmp) == PMINRESSolverQLP<scalarD>::SUCCESSFUL);
                INFOV("Model Derivative Basis: %f %f",RHS.norm(),tmp.norm())
                deriv.col(k++)=tmp.normalized();
                if(k == deriv.cols())
                    break;
            }
        }
        BASIS=deriv;
        boost::dynamic_pointer_cast<FEMReducedSystem>(_body._system)->writeBasisVTK("./basisDeriv/",1.0f,true);
    }
    virtual void precompute() {
        findCoef(_body._tree.get<bool>("overwriteSTVK"));
        //debugCoef();
    }
    virtual Vec calculateFI(const Vec& S,Matd& K) const {
        Vec f=findFK(S,K);
        K*=-1.0f;
        return f;
    }
protected:
    Vec calculateDKDUjUi(const Vec& Uj,const Vec& Ui) const {
#define DIFF_DELTA 0.5f
        Vec DKDjUi;
        Eigen::SparseMatrix<scalarD,0,sizeType> K(Uj.size(),Uj.size());
        {
            ReducedKCalculator::TRIPS trips;
            _body.setDPos(Uj*DIFF_DELTA);
            _energys.materialEnergyEval(NULL,&trips,0.0f,1.0f);
            _body._system->removeConstraint(&trips,NULL,NULL,0);
            K.setFromTriplets(trips.begin(),trips.end());
        }
        DKDjUi=K*Ui;
        {
            ReducedKCalculator::TRIPS trips;
            _body.setDPos(-Uj*DIFF_DELTA);
            _energys.materialEnergyEval(NULL,&trips,0.0f,1.0f);
            _body._system->removeConstraint(&trips,NULL,NULL,0);
            K.setFromTriplets(trips.begin(),trips.end());
        }
        DKDjUi-=K*Ui;
        DKDjUi/=(DIFF_DELTA*2.0f);
        return DKDjUi;
#undef DIFF_DELTA
    }
    virtual Cold findF(const Cold& S) const {
        const Matd& U=_body._basis->_U;
        Vec BLKS=S.block(0,0,nrX(),1);
        Vec dpos=U*BLKS;
        _body.setDPos(dpos);

        Vec f=Vec::Zero(_body.nrV()*3);
        _energys.materialEnergyEval(&f,NULL,1.0f,1.0f);
        return U.transpose()*f;
    }
    virtual sizeType nrF() const {
        return _body._basis->_U.cols();
    }
    virtual sizeType nrX() const {
        return _body._basis->_U.cols();
    }
    //data
    Matd _K;
    Vec _S0,_F0;
    const boost::unordered_map<sizeType,Vec3>& _constraints;
};
//Cubature K calculator
class ForceCubatureProb : public FEMCubatureProb
{
public:
    ForceCubatureProb(const FEMBody& body,const EnergyPool& pool)
        :_body(body),_U0(body._basis->_U) {
        for(EnergyPool::ConstIterator
                beg=pool.begin(typeid(MaterialEnergy)),
                end=pool.end(typeid(MaterialEnergy)); beg!=end; beg++) {
            boost::shared_ptr<MaterialEnergy> me=
                boost::dynamic_pointer_cast<MaterialEnergy>(*beg);
            _meMap[me->getCell()->_index]=me;
        }
    }
    Matd operator()(const sizeType cid,const Matd& S) const {
        Matd U=Matd::Zero(12,_U0.cols());
        FEMCell& cell=*(_body.getCPtr(cid));
        MatrixBasisExtractor UBlk(_U0,Vec::Zero(_U0.cols()));
        for(sizeType i=0; i<_body.dim()+1; i++)
            U.block(i*3,0,3,_U0.cols())=UBlk(cell._v[i]->_index);

        //build the pose matrix
        Matd ret(_U0.cols()*S.cols(),1);
        const MaterialEnergy& me=*(_meMap.find(cid)->second);
        Eigen::Matrix<scalarD,12,1> dpos,force;
        for(sizeType pose=0; pose<S.cols(); pose++) {
            dpos=U.block(0,0,12,S.rows())*S.col(pose);
            for(char v=0; v<_body.dim()+1; v++)
                cell._v[v]->_pos=cell._v[v]->_pos0+dpos.block<3,1>(v*3,0).cast<scalar>();

            force.setZero();
            me.eval(force,NULL);
            ret.block(pose*_U0.cols(),0,_U0.cols(),1)=U.transpose()*force;
        }
        return ret;
    }
    Vec2i nr() const {
        return Vec2i(_U0.cols(),1);
    }
    //data
    boost::unordered_map<sizeType,boost::shared_ptr<MaterialEnergy> > _meMap;
    const FEMBody& _body;
    const Matd& _U0;
};
class ReducedKCalculatorCubature : public ReducedKCalculator
{
public:
    ReducedKCalculatorCubature(FEMReducedSystem& sys)
        :ReducedKCalculator(sys) {
        if(_body._tree.find("overwriteCub") == _body._tree.not_found())
            _body._tree.put<bool>("overwriteCub",false);
    }
    virtual void precompute() {
        //compute cubature
        {
            //solve cubature problem
            ForceCubatureProb obj(_body,_energys);

            MaterialEnergy::TYPE type=_energys.getMaterialType();
            const_cast<EnergyPool&>(_energys).setMaterialType(MaterialEnergy::STVK);

            _body._tree.put<std::string>("pathCub","meshCubaturePotential");
            FEMCubatureSolver sol(_body,obj);
            sol.generate(Vec::Ones(_body.nrC()));
            TBEG("Compute Cubature")
            sol.solve(_body._cubaturePotential,".CUBP");
            TEND

            const_cast<EnergyPool&>(_energys).setMaterialType(type);
            //debugWriteProb(boost::filesystem::ofstream("./meshCubatureForce/prob.dat",ios::binary),_body,sol,obj);
        }

        //get cubature cells
        const Cubature& cub=*(_body._cubaturePotential);
        _body.clearPos();
        _cBody.reset(new CFEMBody(_body,cub));

        //get cubature energy
        _ess.resize(_cBody->nrC());
        for(sizeType i=0; i<_cBody->nrC(); i++) {
            _cBody->getC(i)._mass=1.0f;
            ostringstream ossb;
            ossb << "cell" << cub._tet[i*(_body.dim()+1)]._id;
            _ess[i]=MaterialEnergy(_cBody->getCPtr(i),_body._tree.get_child(ossb.str()));
        }

        //get cubature basis
        const SparseReducedBasis& basis=*(_body._basis);
        _UDense=boost::dynamic_pointer_cast<CFEMBody>(_cBody)->_cInterp*_body._basis->_U;
    }
    virtual Vec calculateFI(const Vec& S,Matd& K) const {
        _cBody->setDPos(_UDense*S.block(0,0,_UDense.cols(),1));
        Vec ret=Vec::Zero(_UDense.cols());
        K.setZero(_UDense.cols(),_UDense.cols());

        Eigen::Matrix<scalarD,12,1> fI;
        Eigen::Matrix<scalarD,12,12> KI;
        Eigen::Matrix<scalarD,12,-1> UBlk(12,_UDense.cols());
        UBlk.setZero();

        sizeType off=(_body.dim()+1)*3;
        sizeType nrC=_cBody->nrC();
        for(sizeType i=0; i<nrC; i++) {
            UBlk.block(0,0,off,_UDense.cols())=_UDense.block(i*off,0,off,_UDense.cols());
            fI.setZero();
            KI.setZero();
            _ess[i].eval(fI,&KI);
            ret+=UBlk.transpose()*fI;
            K+=UBlk.transpose()*(KI*UBlk).eval();
        }
        return ret;
    }
protected:
    boost::shared_ptr<FEMBody> _cBody;
    vector<MaterialEnergy> _ess;
    Matd _UDense;
};
//Cubature K calculator specially designed for Fung-Soft-Tissue energy
class ReducedKCalculatorFung : public ReducedKCalculatorCubature
{
public:
    struct FungTerm {
        Matd _K1,_K2;
        scalarD _c;
    };
    ReducedKCalculatorFung(FEMReducedSystem& system):ReducedKCalculatorCubature(system) {}
    virtual void precompute() {
        ReducedKCalculatorCubature::precompute();

        sizeType nrC=_cBody->nrC();
        sizeType off=(_body.dim()+1)*3;
        Eigen::Matrix<scalarD,12,1> fI;
        Eigen::Matrix<scalarD,12,12> KI;

        _essFung.resize(nrC);
        Eigen::Matrix<scalarD,12,-1> UBlk(12,_UDense.cols());
        UBlk.setZero();
        for(sizeType i=0; i<nrC; i++) {
            UBlk.block(0,0,off,_UDense.cols())=_UDense.block(i*off,0,off,_UDense.cols());
            _ess[i].setType(MaterialEnergy::LINEAR);
            _ess[i].eval(fI,&KI);
            _essFung[i]._K1=UBlk.transpose()*(KI*UBlk).eval();

            scalar mu=_ess[i].getMu();
            scalar lambda=_ess[i].getLambda();
            _ess[i].setMu(_ess[i].getMuFung());
            _ess[i].setLambda(_ess[i].getLambdaFung());
            _ess[i].eval(fI,&KI);
            _essFung[i]._K2=UBlk.transpose()*(KI*UBlk).eval();
            _ess[i].setMu(mu);
            _ess[i].setLambda(lambda);

            _essFung[i]._c=_ess[i].getCFung();
            _ess[i].setType(MaterialEnergy::FUNG);
        }
    }
    virtual Vec calculateFI(const Vec& S,Matd& K) const {
        Vec SBlk=S.block(0,0,_body._basis->_U.cols(),1);
        Vec f=Vec::Zero(SBlk.size()),K1S(SBlk.size()),K2S(SBlk.size());
        K.setZero(SBlk.size(),SBlk.size());
        sizeType nrC=_cBody->nrC();
        for(sizeType i=0; i<nrC; i++) {
            const FungTerm& term=_essFung[i];
            K1S=term._K1*SBlk;
            K2S=term._K2*SBlk;
            scalarD STK2S=std::exp(SBlk.dot(K2S))*term._c;
            f-=(K1S+K2S*STK2S)*_ess[i].getVol();
            K+=(term._K1+STK2S*(term._K2+2.0f*K2S*K2S.transpose()))*_ess[i].getVol();
        }
        return f;
    }
    vector<FungTerm> _essFung;
};
PRJ_END

//FEMReducedSystem
FEMReducedSystem::FEMReducedSystem(FEMBody& body)
    :FEMSystem(body)
{
    if(body._tree.find("reconstruct") == body._tree.not_found())
        body._tree.put<bool>("reconstruct",true);
}
boost::shared_ptr<FEMSystem> FEMReducedSystem::copy(FEMBody& other) const
{
    ASSERT_MSGV(_body.nrV() == other.nrV() && _body.nrC() == other.nrC(),"Bodies don't have same vertex/cell number!");
    boost::shared_ptr<FEMReducedSystem> ret(new FEMReducedSystem(other));
    //ret->_rho=_rho;
    //ret->_MType=_MType;
    //ret->_nrCub=_nrCub;
    //ret->_nrPose=_nrPose;
    ret->energys().copy(_energys);
    ret->energys().lock(false);
    ret->setConstraints(_constraints);
    return boost::dynamic_pointer_cast<FEMSystem>(ret);
}
void FEMReducedSystem::setReducedModel(KINETIC_MODEL_TYPE type,sizeType nrCub,sizeType nrPose,scalar tolSC)
{
    _body._tree.put<sizeType>("type",type);//_MType=type;
    _body._tree.put<sizeType>("nrCub",nrCub);//_nrCub=nrCub;
    _body._tree.put<sizeType>("nrPose",nrPose);//_nrPose=nrPose;
    _body._tree.put<scalar>("tolSC",tolSC);//_tolSC=tolSC;
    _dirty=true;
}
//for integrator
void FEMReducedSystem::buildSystem(scalar MCoef,scalar KCoef,scalar CCoef,FEMSystemMatrix& LHS,FEMSystemMatrix& U,
                                   const Vec& MRHSL,const Vec& CRHS,scalar FCoef,Vec& RHS,sizeType offVar) const
{
    scalar alpha=_body._tree.get<scalar>("alpha");
    scalar beta=_body._tree.get<scalar>("beta");

    RHS=Vec::Zero(size());
    Vec S;
    getPos(S);
    //mass term and matrix U
    Matd M;
    _MCalc->calculateM(S,M,MRHSL,RHS,_body._offset,offVar,_selfColl?&U:NULL);
    LHS.addBlock(offVar,offVar,M*(MCoef+CCoef*alpha));
    RHS.block(0,0,M.rows(),1)+=(M*CRHS.block(0,0,M.cols(),1))*alpha;
    //force and stiffness term
    //this must be put after mass term because
    //sometimes force term depends on it's result
    Matd K;
    Vec F=_KCalc->calculateFI(S,K);
    LHS.addBlock(offVar,offVar,K*(KCoef+CCoef*beta));
    RHS.block(0,0,F.rows(),1)+=F*FCoef;
    RHS.block(0,0,K.rows(),1)+=(K*CRHS.block(0,0,K.cols(),1))*beta;
    //account for non-material force terms
    F=_MCalc->calculateFE(S,K);
    LHS.addBlock(offVar,offVar,K*KCoef);
    RHS.block(0,0,F.rows(),1)+=F*FCoef;
}
//IO
void FEMReducedSystem::getPos(Vec& P) const
{
    if(_MCalc)_MCalc->getPos(P);
}
void FEMReducedSystem::getPosL(Vec& L) const
{
    if(_MCalc)_MCalc->getPosL(L);
}
FEMReducedSystem::Vec FEMReducedSystem::LtoF(const Vec& L) const
{
    if(!_MCalc)return Vec();
    else return _MCalc->LtoF(L);
}
void FEMReducedSystem::setPos(const Vec& P)
{
    if(_dirty)onDirty();
    _MCalc->setPos(P);
}
sizeType FEMReducedSystem::size() const
{
    if(!_MCalc)return 0;
    else return _MCalc->size();
}
sizeType FEMReducedSystem::sizeL() const
{
    if(!_MCalc)return 0;
    else return _MCalc->sizeL();
}
void FEMReducedSystem::debugU()
{
    FEMSystem::debugU();
    _KCalc->debugFI();
}
//(local) basis
void FEMReducedSystem::writeBasisVTK(const std::string& path,scalar off,bool U)
{
    BBox<scalar> bb;
    for(sizeType i=0; i<_body.nrV(); i++)
        bb.setUnion(_body.getV(i)._pos0);
    off*=bb.getExtent().norm();

    boost::filesystem::create_directory(path);
    const SparseReducedBasis& basis=*(_body._basis);
    sizeType nrB=U?basis._U.cols():basis._B.cols();
    for(sizeType i=0; i<nrB; i++) {
        Vec DPos=(U?basis._U.col(i):basis._B.col(i))*off;
        _body.setDPos(DPos);
        sizeType dim=_body.dim();
        while(true) {
            bool flipped=false;
            for(sizeType i=0; i<_body.nrC(); i++) {
                Mat3 F,FN,FR;
                _body.getC(i).buildF(F,FN,FR,NULL,false);
                if(F.block(0,0,dim,dim).determinant() < 0.0f) {
                    flipped=true;
                    break;
                }
            }
            if(!flipped)break;
            DPos*=0.5f;
            _body.setDPos(DPos);
        }

        ostringstream ossVTK;
        ossVTK << path << "/basis" << i << ".vtk";
        VTKWriter<scalar> os("Basis",ossVTK.str(),true);
        _body.writeVTK(os);

        ObjMesh obj;
        _body.writeObj(obj);
        obj.getIG().clear();

        ostringstream ossObj;
        ossObj << path << "/basis" << i << ".obj";
        obj.write(ossObj.str());

        ostringstream ossPov;
        ossPov << path << "/basis" << i << ".pov";
		boost::filesystem::ofstream f(ossPov.str());
        obj.writePov(f);
    }
}
void FEMReducedSystem::buildU()
{
    buildU(_body._tree.get<sizeType>("nrU"),_body._tree.get<bool>("overwriteU",false));
}
void FEMReducedSystem::buildU(sizeType nr,bool overwrite)
{
    if(!_body._basis)_body._basis.reset(new SparseReducedBasis);
    SparseReducedBasis& basis=*(_body._basis);

    //reset to rest state
    Vec pos;
    _body.clearPos();
    _body.getPos(pos);

    //calculate K
    TRIPS trips;
    _energys.materialEnergyEval(NULL,&trips,1.0f,1.0f);
    boost::unordered_set<sizeType> fixSet;
    bool excludeRigid=removeConstraint(&trips,NULL,&fixSet,0);
    basis._K.resize(_body.nrV()*3,_body.nrV()*3);
    basis._K.setFromTriplets(trips.begin(),trips.end());

    //calculate M
    trips.clear();
    Vec M;
    _body.getMass(M);
    M.array()*=_rho.array();
    addI(trips,0,0,M);
    basis._M.resize(_body.nrV()*3,_body.nrV()*3);
    basis._M.setFromTriplets(trips.begin(),trips.end());

    //solve
    if(overwrite || basis._U.cols() == 0 || basis._U.rows() != _body.nrV()*3) {
        Vec lambda(nr);
        EigenSolver::solveEigen(_body,fixSet,lambda,excludeRigid);
        removeRigidMode(basis,lambda);

        nr=basis._U.cols();
        for(sizeType i=0; i<basis._U.cols(); i++) {
            if(basis._U.col(i).sum() < 0.0f)
                basis._U.col(i)*=-1.0f;	//so that we ensure each time basis is the same
            basis._U.col(i).normalize();
        }
        writeBasisVTK("./basisLMA/",1.0f,true);
    } else {
        nr=basis._U.cols();
    }
    _dirty=true;
}
void FEMReducedSystem::onDirty()
{
    _dirty=false;
    _body.getMass(_M);
    _M.array()*=_rho.array();

    //choose model
    sizeType type=_body._tree.get<sizeType>("type");
    switch((KINETIC_MODEL_TYPE)type) {
    case LMA:
        _MCalc.reset(new ReducedMCalculator(*this));
        break;
    case WARPPED_LMA:
        _MCalc.reset(new WrappedReducedMCalculator(*this));
        break;
    case POST_PROCESS_RS:
        _MCalc.reset(new CoupledRSReducedMCalculator(*this));
        _body._tree.put<bool>("coupled",false);
        break;
    case COUPLED_RS:
        _MCalc.reset(new CoupledRSReducedMCalculator(*this));
        _body._tree.put<bool>("coupled",true);
        break;
    case POST_PROCESS_REDUCED_RS:
        if(_body._tree.get<bool>("useFastRSImpl",true))
            _MCalc.reset(new CRCoupledRSReducedMCalculator(*this));
        else _MCalc.reset(new CRCoupledRSReducedMCalculatorRef(*this));
        _body._tree.put<bool>("coupled",false);
        break;
    case CUBATURE_COUPLED_REDUCED_RS:
        if(_body._tree.get<bool>("useFastRSImpl",true))
            _MCalc.reset(new CRCoupledRSReducedMCalculator(*this));
        else _MCalc.reset(new CRCoupledRSReducedMCalculatorRef(*this));
        _body._tree.put<bool>("coupled",true);
        break;
    case TAYLOR_COUPLED_REDUCED_RS:
        _MCalc.reset(new TRCoupledRSReducedMCalculator(*this));
        _body._tree.put<bool>("coupled",true);
        break;
    }

    //choose potential model
    if(_energys.getMaterialType() == MaterialEnergy::LINEAR)
        _KCalc.reset(new ReducedKCalculator(*this));
    else if(_body._tree.get<sizeType>("nrCub",0) > 0) {
        if(_energys.getMaterialType() == MaterialEnergy::FUNG)
            _KCalc.reset(new ReducedKCalculatorFung(*this));
        else _KCalc.reset(new ReducedKCalculatorCubature(*this));
    } else if(_energys.getMaterialType() == MaterialEnergy::STVK)
        _KCalc.reset(new ReducedKCalculatorSTVK(*this));
    else ASSERT_MSG(false,"Unsupported energy type!")

        //build basis
        _MCalc->buildBasis();
    _KCalc->buildBasis();
    if(_body._tree.get<bool>("debugBasisMode",false))
        return;
    SparseReducedBasis& basis=*(_body._basis);
    basis._UTMU=basis._U.transpose()*(basis._M*basis._U).eval();
    basis._UTKU=basis._U.transpose()*(basis._K*basis._U).eval();

    _MCalc->precompute();
    _KCalc->precompute();

    //free moving body wrapper
    if(_body._tree.get<bool>("rigidDOF",false)) {
        boost::shared_ptr<ReducedMCalculator> innerMCalc=_MCalc;
        _MCalc.reset(new RigidReducedMCalculator(*this,innerMCalc));
    }

    //write preprocessed result
    if(_body._tree.find("writeMesh") != _body._tree.not_found()) {
        boost::filesystem::path path=_body._tree.get<std::string>("writeMesh");
        boost::filesystem::create_directory(path.parent_path());
        _body.writeABQ(path.string());
    }

    //reset pos
    setPos(Vec::Zero(size()));
    setVelL(Vec::Zero(sizeL()));
    setAccelL(Vec::Zero(sizeL()));
}
//private
void FEMReducedSystem::removeRigidMode(SparseReducedBasis& basis,Vec& lambda) const
{
    Matd& U=basis._U;
    sizeType nrSmall=0;
    for(; nrSmall<6; nrSmall++) {
        ASSERT_MSG(nrSmall < lambda.size(),"Basis number too small!")
        if(lambda[nrSmall] > 1E-3f)
            break;
    }
    Matd tmpU=U.block(0,nrSmall,U.rows(),U.cols()-nrSmall);
    Matd tmpLambda=lambda.block(nrSmall,0,lambda.size()-nrSmall,1);
    U=tmpU;
    lambda=tmpLambda;
}
