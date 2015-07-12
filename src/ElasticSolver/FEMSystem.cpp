#include "FEMSystem.h"
#include "FEMCollision.h"
#include "FEMCollider.h"
#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMRotationUtil.h"
#include "FEMUtils.h"
#include "ParticleSet.h"

USE_PRJ_NAMESPACE

//FEMFRAME
bool FEMFrameData::read(std::istream& is)
{
    readBinaryData(_X,is);
    readBinaryData(_VL,is);
    readBinaryData(_AL,is);
    return is.good();
}
bool FEMFrameData::write(std::ostream& os) const
{
    writeBinaryData(_X,os);
    writeBinaryData(_VL,os);
    writeBinaryData(_AL,os);
    return os.good();
}

//FEMSYSTEMMATRIX
FEMSystemMatrix::FEMSystemMatrix() {}
bool FEMSystemMatrix::isDense() const
{
    return (_dense!=NULL);
}
FEMSystemMatrix::Vec FEMSystemMatrix::solve(const Vec& RHS) const
{
    if(_dense)return _dense->llt().solve(RHS);
    else {
        Eigen::SparseMatrix<scalarD,0,sizeType> sparse;
        getSparse(sparse);
        return Eigen::SimplicialCholesky<Eigen::SparseMatrix<scalarD,0,sizeType> >(sparse).solve(RHS);
    }
}
void FEMSystemMatrix::getDense(Matd& dense) const
{
    if(_dense)dense=*_dense;
    else if(_sparse) {
        Eigen::SparseMatrix<scalarD,0,sizeType> sparse;
        getSparse(sparse);
        dense=sparse.toDense();
    }
}
void FEMSystemMatrix::getSparse(Eigen::SparseMatrix<scalarD,0,sizeType>& sparse) const
{
    if(_dense) {
        TRIPS trips;
        ::addBlock(trips,0,0,*_dense);
        sparse.resize(_dense->rows(),_dense->cols());
        sparse.setFromTriplets(trips.begin(),trips.end());
    } else if(_sparse) {
        sparse.resize(_size[0],_size[1]);
        sparse.setFromTriplets(_sparse->begin(),_sparse->end());
    }
}
void FEMSystemMatrix::clear()
{
    if(_dense)_dense->setZero();
    if(_sparse)_sparse->clear();
}
void FEMSystemMatrix::reset(sizeType r,sizeType c,bool dense)
{
    _size=Vec2i(r,c);
    if(dense) {
        _dense.reset(new Matd(r,c));
        _dense->setZero();
    } else _sparse.reset(new TRIPS);
}
void FEMSystemMatrix::addBlock(sizeType r,sizeType c,const Matd& block)
{
    if(_dense)_dense->block(r,c,block.rows(),block.cols())+=block;
    else if(_sparse)::addBlock(*_sparse,r,c,block);
}
void FEMSystemMatrix::addBlock(sizeType r,sizeType c,const TRIPS& trips,scalarD coef)
{
    if(_dense)
        for(sizeType i=0; i<(sizeType)trips.size(); i++)
            (*_dense)(r+trips[i].row(),c+trips[i].col())+=trips[i].value()*coef;
    else if(_sparse)
        for(sizeType i=0; i<(sizeType)trips.size(); i++)
            _sparse->push_back(Eigen::Triplet<scalarD,sizeType>
                               (r+trips[i].row(),c+trips[i].col(),trips[i].value()*coef));
}
void FEMSystemMatrix::addI(sizeType r,sizeType c,const Vec& diag,scalarD coef)
{
    if(_dense)
        for(sizeType i=0; i<diag.size(); i++)
            (*_dense)(r+i,c+i)+=diag[i]*coef;
    else if(_sparse)
        ::addI(*_sparse,r,c,diag*coef);
}
void FEMSystemMatrix::addDense(sizeType r,sizeType c,const Matd& other)
{
    if(_dense)_dense->block(r,c,other.rows(),other.cols())+=other;
    else if(_sparse)::addBlock(*_sparse,r,c,other);
}
void FEMSystemMatrix::addSparse(sizeType r,sizeType c,const Eigen::SparseMatrix<scalarD,0,sizeType>& other)
{
    if(_dense)_dense->block(r,c,other.rows(),other.cols())+=other.toDense();
    else if(_sparse)::addSparseBlock(*_sparse,r,c,other);
}

//FEMSYSTEM
//constructor
FEMSystem::FEMSystem(FEMBody& body)
    :_body(body),_energys(body),_dirty(true)
{
    resetParam(0.0f,0.0f,1E3f);
}
FEMBody& FEMSystem::body()
{
    return _body;
}
void FEMSystem::setDamping(scalar alpha,scalar beta)
{
    _body._tree.put<scalar>("alpha",alpha);
    _body._tree.put<scalar>("beta",beta);
}
void FEMSystem::resetParam(scalar alpha,scalar beta,scalar rho)
{
    setDamping(alpha, beta);
    _rho.setConstant(_body.nrV()*3,rho);
}
boost::shared_ptr<FEMSystem> FEMSystem::copy(FEMBody& other) const
{
    ASSERT_MSGV(_body.nrV() == other.nrV() && _body.nrC() == other.nrC(),"Bodies don't have same vertex/cell number!");
    boost::shared_ptr<FEMSystem> ret(new FEMSystem(other));
    ret->_rho=_rho;
    ret->energys().copy(_energys);
    ret->energys().lock(false);
    ret->setConstraints(_constraints);
    return ret;
}
void FEMSystem::applyTrans(const Mat4& R)
{
    for(boost::unordered_map<sizeType,Vec3>::iterator
            beg=_constraints.begin(),end=_constraints.end(); beg!=end; beg++)
        beg->second=R.block<3,3>(0,0)*beg->second+R.block<3,1>(0,3);
}
void FEMSystem::beforeCollision()
{
    if(_dirty)onDirty();
    _selfColl=false;
}
void FEMSystem::onCollision(sizeType id) {}
void FEMSystem::onCollisionSelf(sizeType id)
{
    _selfColl=true;
    onCollision(id);
}
void FEMSystem::afterCollision() {}
//for integrator
void FEMSystem::getSystemEnergy(scalarD& EK,scalarD& EP) const
{
    Vec V=LtoF(_VL);
    EK=0.5f*(V.array()*V.array()*_M.array()).sum();
    EP=_energys.energyEval(NULL,NULL,1.0f,1.0f,EnergyPool::TypeChooser(),false);
}
void FEMSystem::buildSystem(scalar MCoef,scalar KCoef,scalar CCoef,FEMSystemMatrix& LHS,FEMSystemMatrix& U,
                            const Vec& MRHSL,const Vec& CRHS,scalar FCoef,Vec& RHS,sizeType offVar) const
{
    scalar alpha=_body._tree.get<scalar>("alpha");
    scalar beta=_body._tree.get<scalar>("beta");

    TRIPS trips;
    RHS=Vec::Zero(size());
    //force and stiffness term
    _energys.materialEnergyEval(&RHS,&trips,FCoef,1.0f);
    Eigen::SparseMatrix<scalarD,0,sizeType> K(_body.nrV()*3,_body.nrV()*3);
    K.setFromTriplets(trips.begin(),trips.end());
    RHS+=(K*CRHS)*beta;
    for(sizeType i=0; i<(sizeType)trips.size(); i++)
        const_cast<scalarD&>(trips[i].value())*=(KCoef+CCoef*beta);
    //account for non-material force terms
    _energys.nonMaterialEnergyEval(&RHS,&trips,FCoef,KCoef);
    //mass term and U
    RHS+=(_M.array()*(MRHSL+CRHS*alpha).array()).matrix();
    addI(trips,0,0,_M*(MCoef+CCoef*alpha));
    //account for constraint
    removeConstraint(&trips,&RHS,NULL,0);
    LHS.addBlock(offVar,offVar,trips,1.0f);
    if(_selfColl)
        U.addI(_body._offset,offVar,Vec::Ones(_M.size()),1.0f);
}
void FEMSystem::buildU(Eigen::SparseMatrix<scalarD,0,sizeType>& U)
{
    FEMSystemMatrix mLHS,mU;
    Vec MRHSL(sizeL()),CRHS(size()),RHS(size());
    mU.reset(_body.nrV()*3,size(),false);
    _selfColl=true;
    _body._offset=0;
    buildSystem(0.0f,0.0f,0.0f,mLHS,mU,MRHSL,CRHS,0.0f,RHS,0);
    mU.getSparse(U);
}
sizeType FEMSystem::size() const
{
    return _body.nrV()*3;
}
sizeType FEMSystem::sizeL() const
{
    return size();
}
bool FEMSystem::selfColl() const
{
    return _selfColl;
}
void FEMSystem::debugU()
{
    if(_dirty)onDirty();
    beforeCollision();
    onCollision(rand()%_body.nrSV());
    afterCollision();
    //save P
    Vec Ptmp;
    getPos(Ptmp);
    boost::unordered_map<sizeType,Vec3> constraints=_constraints;
    _constraints.clear();
    //test at random P
    Vec P=Ptmp;
    P.setRandom();
    Vec X0,X1;
    setPos(P);
    _body.getPos(X0);
    //calculate U
    Eigen::SparseMatrix<scalarD,0,sizeType> U;
    buildU(U);
    //test each component of P
    for(sizeType i=0; i<P.size(); i++) {
#define DELTA 1E-7f
        scalarD Pitmp=P[i];
        P[i]+=DELTA;
        setPos(P);
        _body.getPos(X1);
        Vec UiN=(X1-X0)/DELTA;
        for(sizeType j=0; j<UiN.rows(); j++) {
            INFOV("A: %f,N: %f,ERR: %f",U.coeff(j,i),UiN[j],std::abs(U.coeff(j,i)-UiN[j]))
            if(std::abs(UiN[j]) > 1E-5f && std::abs((U.coeff(j,i)-UiN[j])/UiN[j]) > 1E-3f)
                ASSERT_MSG(false,"U error!");
        }
        P[i]=Pitmp;
#undef DELTA
    }
    //load P
    setPos(Ptmp);
    _constraints=constraints;
}
//IO
void FEMSystem::getFrame(FEMFrameData& data)
{
    if(_dirty)onDirty();
    getPos(data._X);
    getVelL(data._VL);
    getAccelL(data._AL);
}
void FEMSystem::setFrame(const FEMFrameData& data)
{
    if(_dirty)onDirty();
    setPos(data._X);
    setVelL(data._VL);
    setAccelL(data._AL);
}
void FEMSystem::getPos(Vec& P) const
{
    _body.getPos(P);
}
void FEMSystem::getPosL(Vec& L) const
{
    getPos(L);
}
void FEMSystem::getVelL(Vec& VL) const
{
    VL=_VL;
}
void FEMSystem::getAccelL(Vec& AL) const
{
    AL=_AL;
}
FEMSystem::Vec FEMSystem::LtoF(const Vec& F) const
{
    return F;
}
void FEMSystem::setPos(const Vec& P)
{
    _body.setPos(P);
    applyConstraint();
}
void FEMSystem::setVelL(const Vec& VL)
{
    _VL=VL;
}
void FEMSystem::setAccelL(const Vec& AL)
{
    _AL=AL;
}
//energy utility
void FEMSystem::addEnergyMaterial(scalar young,scalar poisson,MaterialEnergy::TYPE type,bool invertible)
{
    sizeType nrC=(sizeType)_body.nrC();
    for(sizeType i=0; i<nrC; i++) {
        ostringstream oss;
        oss << "cell" << i;
        _body._tree.put<scalar>(oss.str()+".young",young);
        _body._tree.put<scalar>(oss.str()+".poisson",poisson);
        _body._tree.put<sizeType>(oss.str()+".type",type);
        _body._tree.put<bool>(oss.str()+".invertible",invertible);
        boost::shared_ptr<FEMEnergy> e(new MaterialEnergy(_body.getCPtr(i),_body._tree.get_child(oss.str())));
        //e->debugEnergy(_body);
        boost::dynamic_pointer_cast<MaterialEnergy>(e)->setInvertible(invertible);
        _energys.addEnergy(e);
    }
    _dirty=true;
}
void FEMSystem::setMaterialType(MaterialEnergy::TYPE type)
{
    _energys.setMaterialType(type);
    _dirty=true;
}
void FEMSystem::addEnergyMass()
{
    Vec3 gravity
    (_body._tree.get<scalar>("gravX",0.0f),
     _body._tree.get<scalar>("gravY",0.0f),
     _body._tree.get<scalar>("gravZ",0.0f));
    addEnergyMass(gravity);
}
void FEMSystem::addEnergyMass(const Vec3& g,const BBox<scalar>* bbm)
{
    sizeType nrV=_body.nrV();
    for(sizeType i=0; i<nrV; i++)
        if(!bbm || bbm->contain(_body.getV(i)._pos,_body.dim())) {
            boost::shared_ptr<FEMEnergy> e(new FixedForceEnergy(_body.getVPtr(i),g*(scalar)_rho[i]));
            _energys.addEnergy(e);
        }
}
void FEMSystem::addEnergy(boost::shared_ptr<FEMEnergy> e)
{
    _energys.addEnergy(e);
}
void FEMSystem::delEnergy(boost::shared_ptr<FEMEnergy> e)
{
    _energys.delEnergy(e);
}
void FEMSystem::addForcePoint(sizeType id, const Vec3& force)
{
    boost::shared_ptr<FEMEnergy> e(new FixedForceEnergy(_body.getVPtr(id),force));
    _energys.addEnergy(e);
}
void FEMSystem::delForcePoint(sizeType id)
{
    _energys.delEnergy(typeid(FixedForceEnergy).name(),id);
}
EnergyPool& FEMSystem::energys()
{
    return _energys;
}
void FEMSystem::clearEnergy()
{
    _energys.clear();
    _dirty=true;
}
void FEMSystem::clearEnergy(const std::string& name)
{
    _energys.clear(name);
}
void FEMSystem::readEnergy(const std::string& path)
{
    MaterialEnergy::TYPE type=(MaterialEnergy::TYPE)_body._tree.get<sizeType>("mtype");
    bool invertible=_body._tree.get<bool>("invertible",false);
    bool constraints=_body._tree.get<bool>("readConstraints",false);
    readEnergy(path,type,invertible,constraints);
}
void FEMSystem::readEnergy(const std::string& path,MaterialEnergy::TYPE type,bool invertible,bool constraints)
{
    sizeType nrGroup=0,nrTet=0;
    Vec weight=Vec::Zero(_body.nrV()*3);
    _rho.setZero();

    std::string line;
    boost::filesystem::ifstream is(path);
    boost::unordered_map<sizeType,std::pair<Vec3,sizeType> > groups;
    while(getline(is,line).good()) {
        int nr,id,cvid;
        std::string tmp;
        char comma;
#ifdef _MSC_VER
        nr=sscanf_s(line.c_str(),"*MATERIAL, NAME=GROUP_%d",&id);
#else
        nr=sscanf(line.c_str(),"*MATERIAL, NAME=GROUP_%d",&id);
#endif
        if(nr == 1) {
            nrGroup++;
            std::pair<Vec3,sizeType> info;
            getline(is,line);	//type
            getline(is,line);	//young poisson
            {
                istringstream iss(line);
                iss >> info.first[0] >> comma >> info.first[1];
            }
            getline(is,line);	//density
            {
                istringstream iss(line);
                iss >> tmp >> info.first[2];
            }
            line.clear();
            getline(is,line);	//nr tet
            if(line.empty()) {
                _rho.setConstant(info.first[2]);
                weight.setConstant(1.0f);
                addEnergyMaterial(info.first[0],info.first[1],type,invertible);
                nrTet=_body.nrC();
            } else {
                istringstream iss(line);
                iss >> tmp >> info.second;
            }
            groups[id]=info;
            continue;
        }

#ifdef _MSC_VER
        nr=sscanf_s(line.c_str(),"*ELSET, ELSET=GROUP_%d",&id);
#else
        nr=sscanf(line.c_str(),"*ELSET, ELSET=GROUP_%d",&id);
#endif
        if(nr == 1) {
            std::pair<Vec3,sizeType> info=groups[id];
            nrTet+=info.second;
            getline(is,line);	//tets
            istringstream iss(line);
            for(sizeType i=0; i<info.second; i++) {
                iss >> id >> comma;
                ostringstream oss;
                oss << "cell" << i;
                _body._tree.put<scalar>(oss.str()+".young",info.first[0]);
                _body._tree.put<scalar>(oss.str()+".poisson",info.first[1]);
                _body._tree.put<sizeType>(oss.str()+".type",type);
                _body._tree.put<bool>(oss.str()+".invertible",false);
                boost::shared_ptr<FEMEnergy> e(new MaterialEnergy(_body.getCPtr(id),_body._tree.get_child(oss.str())));
                boost::dynamic_pointer_cast<MaterialEnergy>(e)->setInvertible(invertible);
                _energys.addEnergy(e);	//energy

                const FEMCell& c=_body.getC(id);
                for(sizeType v=0; v<4; v++)
                    if(c._v[v]) {
                        weight.block<3,1>(c._v[v]->_index*3,0).array()+=1.0f;
                        _rho.block<3,1>(c._v[v]->_index*3,0).array()+=info.first[2];	//density
                    }
            }
            continue;
        }

#ifdef _MSC_VER
        nr=sscanf_s(line.c_str(),"*CSET, CONS=%d",&id);
#else
        nr=sscanf(line.c_str(),"*CSET, CONS=%d",&id);
#endif
        Vec3 pos;
        if(nr == 1)
            for(sizeType i=0; i<id; i++) {
                getline(is,line);
                istringstream(line) >>
                                    cvid >> comma >>
                                    pos[0] >> comma >>
                                    pos[1] >> comma >>
                                    pos[2];
                _constraints[cvid]=pos;
            }
    }
    _rho.array()/=weight.array();
    INFOV("%d groups, %d cells, %d constraints",nrGroup,nrTet,_constraints.size());
    ASSERT_MSG(nrTet == _body.nrC(),"Cell number mismatch!")
}
void FEMSystem::writeEnergy(const std::string& pathIn,const std::string& pathOut) const
{
    sizeType nr,id;
    ostringstream oss;
    std::string line;
    boost::filesystem::ifstream is(pathIn);
    while(getline(is,line)) {
#ifdef _MSC_VER
        nr=sscanf_s(line.c_str(),"*CSET, CONS=%d",&id);
#else
        nr=sscanf(line.c_str(),"*CSET, CONS=%d",&id);
#endif
        if(nr == 1)
            for(sizeType i=0; i<id; i++)
                getline(is,line);
        else oss << line << std::endl;
        if(!is.good())
            break;
    }
    oss << std::endl;

    ParticleSet css;
    oss << "*CSET, CONS=" << _constraints.size() << std::endl;
    for(boost::unordered_map<sizeType,Vec3>::const_iterator
            beg=_constraints.begin(),end=_constraints.end(); beg!=end; beg++) {
        oss << beg->first << ',' <<
            beg->second[0] << ',' <<
            beg->second[1] << ',' <<
            beg->second[2] << std::endl;

        Particle<scalar> p;
        p._pos=beg->second;
        css.addParticle(p);
    }
    css.writeVTK("./constraints.vtk");
    boost::filesystem::ofstream(pathOut) << oss.str();
}
//constraint utility
void FEMSystem::setConstraints(const boost::unordered_map<sizeType,Vec3>& cons)
{
    _constraints=cons;
    _dirty=true;
}
void FEMSystem::addConstraintPoint(FEMInterp p,scalar K)
{
    boost::shared_ptr<FEMCell> c=_body.getCPtr(p._id);
    boost::shared_ptr<FEMEnergy> e(new FixPointEnergy(_body.dim(),c,p,K,c->getVert(p)));
    addEnergy(e);
    _dirty=true;
}
void FEMSystem::addConstraintPoint(sizeType id, const Vec3& pos)
{
    _constraints[id]=pos;
    clearVel(id);
    clearAccel(id);
    _dirty=true;
}
void FEMSystem::addConstraintPoint(sizeType id)
{
    _constraints[id]=_body.getV(id)._pos;
    clearVel(id);
    clearAccel(id);
    _dirty=true;
}
void FEMSystem::delConstraintPoint(sizeType id)
{
    _constraints.erase(id);
    _dirty=true;
}
boost::unordered_map<sizeType,Vec3>& FEMSystem::constraints()
{
    return _constraints;
}
void FEMSystem::clearConstraint()
{
    _constraints.clear();
    _dirty=true;
}
void FEMSystem::onDirty()
{
    _dirty=false;
    _body.getMass(_M);
    _M.array()*=_rho.array();
    _body.clearPos();
    Vec X;
    _body.getPos(X);
    setPos(X);
    setVelL(Vec::Zero(sizeL()));
    setAccelL(Vec::Zero(sizeL()));
}
//private
bool FEMSystem::removeConstraint(TRIPS* LHS,Vec* RHS,boost::unordered_set<sizeType>* fixSetRet,sizeType offsetLHS) const
{
    boost::unordered_set<sizeType> fixSet;
    for(boost::unordered_map<sizeType,Vec3>::const_iterator
            beg=_constraints.begin(),end=_constraints.end(); beg!=end; beg++) {
        const FEMVertex& v=_body.getV(beg->first);
        fixSet.insert(v._index*3+0);
        fixSet.insert(v._index*3+1);
        fixSet.insert(v._index*3+2);
    }
    bool excludeRigid=fixSet.empty();
    for(sizeType i=0; i<_body.nrV(); i++)
        if(_body.dim() == 2)
            fixSet.insert(i*3+2);

    if(RHS)
        for(boost::unordered_set<sizeType>::const_iterator beg=fixSet.begin(),end=fixSet.end(); beg!=end; beg++)
            (*RHS)[*beg]=0.0f;

    if(LHS)
        for(sizeType i=offsetLHS; i<(sizeType)LHS->size();) {
            if((*LHS)[i].row() == (*LHS)[i].col()) {
                i++;
                continue;
            }
            if(fixSet.find((*LHS)[i].row()) != fixSet.end() ||
                    fixSet.find((*LHS)[i].col()) != fixSet.end()) {
                (*LHS)[i]=LHS->back();
                LHS->pop_back();
            } else i++;
        }
    if(fixSetRet)*fixSetRet=fixSet;
    return excludeRigid;
}
void FEMSystem::applyConstraint()
{
    for(boost::unordered_map<sizeType,Vec3>::const_iterator
            beg=_constraints.begin(),end=_constraints.end(); beg!=end; beg++) {
        _body.getV(beg->first)._pos=beg->second;
        clearVel(beg->first);
        clearAccel(beg->first);
    }
}
void FEMSystem::clearVel(sizeType id)
{
    if(sizeL() == _body.nrV()*3)
        _VL.block<3,1>(id*3,0).setZero();
}
void FEMSystem::clearAccel(sizeType id)
{
    if(sizeL() == _body.nrV()*3)
        _AL.block<3,1>(id*3,0).setZero();
}

//FEMSOLVER
//constructor
FEMSolver::FEMSolver(sizeType dim,sizeType cOption)
{
    boost::shared_ptr<FEMCollision> coll;
    switch(cOption) {
    case 0:
        coll.reset(new FEMCollision);
        break;
    case 1:
        coll.reset(new BVHFEMCollision);
        break;
    case 2:
        coll.reset(new SBVHFEMCollision);
        break;
    }
    _mesh.reset(new FEMMesh(dim,coll));

    _tree.put<bool>("denseSystem",true);
    resetImplicitEuler();
    setCollK(1E5f);
    _mesh->setCellSz(1.0f);
    setSelfColl(true);
}
void FEMSolver::resetParam(scalar beta,scalar beta2,scalar gamma,scalar eps,sizeType maxIter)
{
    _tree.put<scalar>("beta",beta);//_beta=beta;
    _tree.put<scalar>("beta2",beta2);//_beta2=beta2;
    _tree.put<scalar>("gamma",gamma);//_gamma=gamma;
    _tree.put<scalar>("eps",eps);//_eps=eps;
    _tree.put<sizeType>("maxIter",maxIter);//_maxIter=maxIter;
}
void FEMSolver::setCollK(scalar collK)
{
    _tree.put<scalar>("collK",collK);//_collK=collK;
}
FEMMesh& FEMSolver::getMesh()
{
    return *_mesh;
}
const FEMMesh& FEMSolver::getMesh() const
{
    return *_mesh;
}
void FEMSolver::setSelfColl(bool selfColl)
{
    _tree.put<bool>("selfColl",selfColl);//_selfColl=selfColl;
}
//solver utility
#define BLK(IN,I) (IN).block(_offVar[i],0,_offVar[i+1]-_offVar[i],1)
#define BLKL(IN,I) (IN).block(_offVarL[i],0,_offVarL[i+1]-_offVarL[i],1)
#define BLKF(IN,I) (IN).block(_offVarF[i],0,_offVarF[i+1]-_offVarF[i],1)
void FEMSolver::advance(scalar dt)
{
    //initialize
    _mesh->buildOffset();
    scalar collK=_tree.get<scalar>("collK");
    scalar beta2=_tree.get<scalar>("beta2");
    scalar beta=_tree.get<scalar>("beta");
    scalar gamma=_tree.get<scalar>("gamma");
    scalar eps=_tree.get<scalar>("eps");
    sizeType maxIter=_tree.get<sizeType>("maxIter");

    //handle collision detection
    bool selfColl=_tree.get<bool>("selfColl");
    for(sizeType i=0; i<_mesh->nrB(); i++)
        _mesh->getB(i)._system->beforeCollision();
    DefaultFEMCollider coll(*_mesh,collK,1.0f,1.0f);
    if(_geom)_mesh->getColl().collideGeom(*_geom,coll,true);
    if(selfColl)_mesh->getColl().collideMesh(coll,true);
    Vec fE=coll.getFE();
    Eigen::SparseMatrix<scalarD,0,sizeType> HE;
    coll.getHE(HE);
    for(sizeType i=0; i<_mesh->nrB(); i++)
        _mesh->getB(i)._system->afterCollision();

    //find variable offset
    assembleOffVar();

    //initialize velocity and position
    Vec x0(nrVar()),x1(nrVar());
    Vec X0(nrVarL()),X1(nrVarL()),PHI(nrVarL()),PSI(nrVarL());
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        const FEMSystem& sys=*(_mesh->getB(i)._system);
        Vec xB,XB,VB,AB;
        sys.getPos(xB);
        sys.getPosL(XB);
        sys.getVelL(VB);
        sys.getAccelL(AB);

        BLK(x0,i)=xB;
        BLKL(X0,i)=XB;
        BLKL(PHI,i)=VB+(1.0f-gamma)*dt*AB;
        BLKL(PSI,i)=XB+dt*VB+(1.0f-2.0f*beta2)*dt*dt*0.5f*AB;
    }

    if(_tree.get<bool>("debug",false))
        INFO("Newton Iteration")
        //main loop: we use Implicit Newmark Scheme
        Vec RHS(nrVar()),DELTA;
    _LHS.reset(nrVar(),nrVar(),_tree.get<bool>("denseSystem"));
    _U.reset(nrVarF(),nrVar(),_tree.get<bool>("denseSystem"));
    for(sizeType i=0; i<maxIter; i++) {
        _LHS.clear();
        _U.clear();
        for(sizeType i=0; i<_mesh->nrB(); i++) {
            const FEMSystem& sys=*(_mesh->getB(i)._system);
            Vec xB,XB;
            sys.getPos(xB);
            sys.getPosL(XB);
            BLK(x1,i)=xB;
            BLKL(X1,i)=XB;
        }

        //build LHS and RHS
        for(sizeType i=0; i<_mesh->nrB(); i++) {
            //tell system to build: M+(beta*dt*dt)*K+(gamma*dt)*C with off
            //tell system to build: U with offU
            //tell system to build: M(dSn+((1-gamma)*dt)*ddSn-dS_{n+1})+(gamma*dt)*f_I-C*dS_{n+1} with offU[1]
            const FEMSystem& sys=*(_mesh->getB(i)._system);
            Vec RHSB=Vec::Zero(sys.size());
            Vec MRHSLB=BLKL(PSI-X1,i);
            Vec CRHSB=BLK(x0-x1,i)*beta*dt;
            sys.buildSystem(1.0f,beta*dt*dt,beta*dt,_LHS,_U,
                            MRHSLB,CRHSB,beta*dt*dt,RHSB,_offVar[i]);
            BLK(RHS,i)=RHSB;
        }

        //add external term
        if(selfColl)
            if(_U.isDense()) {
                Matd UDense;
                _U.getDense(UDense);
                RHS+=UDense.transpose()*(fE-HE*LtoF(X1-X0))*(beta*dt*dt);
                _LHS.addDense(0,0,UDense.transpose()*(HE*UDense).eval()*(beta*dt*dt));
            } else {
                Eigen::SparseMatrix<scalarD,0,sizeType> USparse;
                _U.getSparse(USparse);
                RHS+=USparse.transpose()*(fE-HE*LtoF(X1-X0))*(beta*dt*dt);
                _LHS.addSparse(0,0,USparse.transpose()*(HE*USparse).eval()*(beta*dt*dt));
            }

        //solve
        DELTA=_LHS.solve(RHS);
        x1+=DELTA;
        for(sizeType i=0; i<_mesh->nrB(); i++) {
            FEMSystem& sys=*(_mesh->getB(i)._system);
            sys.setPos(BLK(x1,i));
        }

        //exit test
        if(DELTA.size() == 0)
            break;
        if(_tree.get<bool>("debug",false))
            INFOV("\tResidue: %f",DELTA.cwiseAbs().maxCoeff())
            if(DELTA.cwiseAbs().maxCoeff() < eps)
                break;
    }

    //update velocity and acceleration
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        FEMSystem& sys=*(_mesh->getB(i)._system);
        Vec XB,VB,AB;
        sys.getPosL(XB);
        AB=(XB-BLKL(PSI,i))/(beta*dt*dt);
        VB=BLKL(PHI,i)+AB*gamma*dt;

        sys.setVelL(VB);
        sys.setAccelL(AB);
    }
    _mesh->updateMesh();
}
FEMSolver::Vec FEMSolver::LtoF(const Vec& L) const
{
    Vec ret(nrVarF());
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        const FEMSystem& sys=*(_mesh->getB(i)._system);
        BLKF(ret,i)=sys.LtoF(BLKL(L,i));
    }
    return ret;
}
#undef BLK
#undef BLKL
#undef BLKF
void FEMSolver::resetBody(sizeType bid)
{
    FEMSystem& sys=getSystem<FEMSystem>(bid);
    sys.beforeCollision();
    sys.afterCollision();
    sys.setPos(Vec::Zero(sys.size()));
    sys.setVelL(Vec::Zero(sys.sizeL()));
    sys.setAccelL(Vec::Zero(sys.sizeL()));
}
void FEMSolver::getSystemEnergy(scalarD& EK,scalarD& EP) const
{
    EK=EP=0.0f;
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        scalarD EKB=0.0f,EPB=0.0f;
        _mesh->getB(i)._system->getSystemEnergy(EKB,EPB);
        EK+=EKB;
        EP+=EPB;
    }
}
//IO
void FEMSolver::writeFrame(std::ostream& os)
{
    sizeType nrB=_mesh->nrB();
    writeBinaryData(nrB,os);
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        FEMFrameData dat;
        _mesh->getB(i)._system->getFrame(dat);
        dat.write(os);
    }
}
void FEMSolver::readFrame(std::istream& is)
{
    sizeType nrB;
    readBinaryData(nrB,is);
    ASSERT_MSG(nrB == _mesh->nrB(),"Body number mismatch!")
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        FEMFrameData dat;
        dat.read(is);
        _mesh->getB(i)._system->setFrame(dat);
    }
    _mesh->updateMesh();
}
//energy utility
class FEMConstraintAdder : public FEMCollider
{
public:
    FEMConstraintAdder(const FEMMesh& mesh):_mesh(mesh) {}
    virtual void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n) {
        b->_system->addConstraintPoint(v->_index);
    }
    const FEMMesh& _mesh;
};
void FEMSolver::addConstraintPoint(FEMGeom& geom)
{
    FEMConstraintAdder temp(*_mesh);
    _mesh->updateMesh();
    _mesh->getColl().collideGeom(geom,temp,false);
}
void FEMSolver::addConstraintPoint(const BBox<scalar>& bb)
{
    FEMGeom geom(_mesh->dim());
    geom.addGeomBox(Mat4::Identity(),bb);
    geom.assemble();
    addConstraintPoint(geom);
}
void FEMSolver::addEnergy(boost::shared_ptr<FixPointEnergy> e)
{
    for(sizeType i=0; i<_mesh->nrB(); i++)
        if(_mesh->getB(i).getCPtr(e->getCell()->_index) == e->getCell()) {
            _mesh->getB(i)._system->addEnergy(e);
            return;
        }
}
void FEMSolver::delEnergy(boost::shared_ptr<FixPointEnergy> e)
{
    for(sizeType i=0; i<_mesh->nrB(); i++)
        if(_mesh->getB(i).getCPtr(e->getCell()->_index) == e->getCell()) {
            _mesh->getB(i)._system->delEnergy(e);
            return;
        }
}
void FEMSolver::assembleOffVar()
{
    _offVar.assign(_mesh->nrB()+1,0);
    _offVarL.assign(_mesh->nrB()+1,0);
    _offVarF.assign(_mesh->nrB()+1,0);
    for(sizeType i=0; i<_mesh->nrB(); i++) {
        const FEMBody& b=_mesh->getB(i);
        FEMSystem& sys=*(b._system);
        _offVar[i+1]=_offVar[i]+sys.size();
        _offVarL[i+1]=_offVarL[i]+sys.sizeL();
        _offVarF[i+1]=_offVarF[i]+b.nrV()*3;
    }
}
sizeType FEMSolver::nrVar() const
{
    return _offVar.back();
}
sizeType FEMSolver::nrVarL() const
{
    return _offVarL.back();
}
sizeType FEMSolver::nrVarF() const
{
    return _offVarF.back();
}
