#include "FEMEnergy.h"
#include "FEMMesh.h"
#include "FEMUtils.h"
#include "FEMSparseReducedBasis.h"

USE_PRJ_NAMESPACE

//BASIS EXTRACTOR
MatrixBasisExtractor::MatrixBasisExtractor(const Matd& U,const Cold& S):_U(U)
{
    _S=S.block(0,0,_U.cols(),1);
}
Matd MatrixBasisExtractor::operator()(sizeType i) const
{
    return _U.block(i*3,0,3,_U.cols());
}
Vec3d MatrixBasisExtractor::operator[](sizeType i) const
{
    Cold BLKs=_S.block(0,0,_U.cols(),1);
    return _U.block(i*3,0,3,_U.cols())*BLKs;
}
//FEMENERGY
void FEMEnergy::debugEnergy(const FEMBody& body)
{
    //force
    Vec f,f2;
    f.resize(body.nrV()*3);
    f2.resize(body.nrV()*3);
    f.setZero();
    //hessian
    TRIPS H;
    scalarD E=eval(&f,&H,-1.0f,1.0f);
    //build
    Eigen::SparseMatrix<scalarD,0,sizeType> M(body.nrV()*3,body.nrV()*3);
    M.setFromTriplets(H.begin(),H.end());
#define DELTA 1E-7f
    //test force
    for(sizeType i=0; i<(sizeType)body.nrV()*3; i++) {
        if(body.dim() == 2 && i%3 == 2)
            continue;
        if(std::abs(f[i]) > 0.0f) {
            Vec3& pt=body.getVPtr(i/3)->_pos;
            Vec3 ptTmp=pt;
            pt[i%3]+=DELTA;
            H.clear();
            scalarD EN=(eval(&f2,&H,-1.0f,1.0f)-E)/DELTA;
            INFOV("A: %f,N: %f,ERR: %f",f[i],EN,std::abs(f[i]-EN))
            if(std::abs((f[i]-EN)/EN) > 1E-3f)
                ASSERT_MSG(false,"Energy error!");
            pt=ptTmp;
        }
    }
    //test hessian
    for(int k=0; k<M.outerSize(); ++k)
        for(Eigen::SparseMatrix<scalarD,0,sizeType>::InnerIterator it(M,k); it; ++it) {
            if(body.dim() == 2 && (it.row()%3 == 2 || it.col()%3 == 2))
                continue;
            Vec3& pt=body.getVPtr(it.col()/3)->_pos;
            Vec3 ptTmp=pt;
            pt[it.col()%3]+=DELTA;
            H.clear();
            f2.setZero();
            eval(&f2,&H,-1.0f,1.0f);
            scalarD EN=(f2[it.row()]-f[it.row()])/DELTA;
            INFOV("A: %f,N: %f,ERR: %f",it.value(),EN,std::abs(it.value()-EN))
            if(std::abs((it.value()-EN)/EN) > 1E-3f)
                ASSERT_MSG(false,"Energy error!");
            pt=ptTmp;
        }
#undef DELTA
}
scalarD FEMEnergy::eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const
{
    ASSERT_MSG(false,"This FEMEnergy doesn't support reduced evaluation!")
    return 0.0f;
}
void FEMEnergy::flagVertex(vector<bool>& flag) const
{
    ASSERT_MSG(false,"This FEMEnergy doesn't support vertex flagging!")
}
boost::shared_ptr<FEMEnergy> FEMEnergy::copy(const FEMBody& other) const
{
    ASSERT_MSG(false,"This FEMEnergy doesn't support interbody copy!")
    return boost::shared_ptr<FEMEnergy>((FEMEnergy*)NULL);
}
void FEMEnergy::evalHelper(const BasisExtractor& U,Vec* f,Matd* H,Eigen::Matrix<scalarD,12,1>* fB,Eigen::Matrix<scalarD,12,12>* HB,const FEMCell& cell)
{
    int nr=cell._v[3] == NULL ? 3 : 4;
    vector<Matd> UBlk(nr);
    for(int i=0; i<nr; i++)
        UBlk[i]=U(cell(i));
    for(int i=0; i<nr; i++) {
        if(f)(*f)+=UBlk[i].transpose()*fB->block<3,1>(i*3,0);
        if(H)for(int j=0; j<nr; j++)
                (*H)+=UBlk[i].transpose()*(HB->block<3,3>(i*3,j*3)*UBlk[j]).eval();
    }
}
void FEMEnergy::evalHelper(Vec* f,TRIPS* H,Eigen::Matrix<scalarD,12,1>* fB,Eigen::Matrix<scalarD,12,12>* HB,const FEMCell& cell)
{
    int nr=cell._v[3] == NULL ? 3 : 4;
    for(int i=0; i<nr; i++) {
        if(f)f->block<3,1>(cell(i)*3,0)+=fB->block<3,1>(i*3,0);
        if(H)for(int j=0; j<nr; j++)
                addBlock(*H,cell(i)*3,cell(j)*3,HB->block<3,3>(i*3,j*3));
    }
}
//force energy
FixedForceEnergy::FixedForceEnergy(boost::shared_ptr<FEMVertex> v,Vec3 a):_v(v)
{
    _f=a*v->_mass;
}
scalarD FixedForceEnergy::eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const
{
    //energy=-v.dot(f)
    if(f)f->block<3,1>(_v->_index*3,0)+=_f.cast<scalarD>()*CF;
    return -_v->_pos.dot(_f);
}
scalarD FixedForceEnergy::eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const
{
    if(f)(*f)+=U(_v->_index).transpose()*_f.cast<scalarD>()*CF;
    return -U[_v->_index].cast<scalar>().dot(_f);
}
void FixedForceEnergy::flagVertex(vector<bool>& flag) const
{
    flag[_v->_index]=true;
}
boost::shared_ptr<FEMEnergy> FixedForceEnergy::copy(const FEMBody& other) const
{
    return boost::shared_ptr<FEMEnergy>(new FixedForceEnergy(other.getVPtr(_v->_index),_f/_v->_mass));
}
//fix vertex energy
FixVertexEnergy::FixVertexEnergy(sizeType dim,boost::shared_ptr<FEMVertex> v,scalar K,const Vec3& p0)
    :_dim(dim),_pos(p0),_K(K),_v(v) {}
scalarD FixVertexEnergy::eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const
{
    Vec3d delta=(_v->_pos-_pos).cast<scalarD>();
    if(f)f->block<3,1>(_v->_index*3,0)-=delta*CF*_K;
    if(H)addI3x3(*H,_v->_index*3,_v->_index*3,CH*_K);
    return 0.5f*delta.squaredNorm()*_K;
}
scalarD FixVertexEnergy::eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const
{
    Vec3d delta=U[_v->_index]-_pos.cast<scalarD>();
    Matd UBlk=U(_v->_index);
    if(f)(*f)-=UBlk.transpose()*delta*CF*_K;
    if(H)(*H)+=UBlk.transpose()*(Mat3d::Identity()*UBlk).eval()*(CH*_K);
    return 0.5f*delta.squaredNorm()*_K;
}
void FixVertexEnergy::flagVertex(vector<bool>& flag) const
{
    flag[_v->_index]=true;
}
boost::shared_ptr<FEMEnergy> FixVertexEnergy::copy(const FEMBody& other) const
{
    return boost::shared_ptr<FEMEnergy>(new FixVertexEnergy(_dim,other.getVPtr(_v->_index),_K,_pos));
}
Vec3 FixVertexEnergy::getPos() const
{
    return _pos;
}
void FixVertexEnergy::setPos(const Vec3& pos)
{
    _pos=pos;
}
boost::shared_ptr<FEMVertex> FixVertexEnergy::getVert() const
{
    return _v;
}
//fix point energy
FixPointEnergy::FixPointEnergy(sizeType dim,boost::shared_ptr<FEMCell> c,const FEMInterp& I,scalar K,const Vec3& p0)
    :_dim(dim),_pos(p0),_K(K),_I(new FEMInterp(I)),_c(c)
{
    _HB.setZero();
    for(int i=0; i<(_dim+1); i++)
        for(int j=0; j<(_dim+1); j++)
            _HB.block<3,3>(i*3,j*3)=Mat3d::Identity()*(I._coef[i]*I._coef[j]*K);
}
scalarD FixPointEnergy::eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const
{
    //energy=0.5*_K*|\sum _I._v[i] _I._b[i] - p|^2
    Vec3 delta=-_pos;
    for(int i=0; i<(_dim+1); i++)
        delta+=(*_c)[i]*_I->_coef[i];

    Eigen::Matrix<scalarD,12,12> HB=_HB*CH;
    Eigen::Matrix<scalarD,12,1> fB;
    for(int i=0; i<(_dim+1); i++)
        if(f)fB.block<3,1>(i*3,0)=-delta.cast<scalarD>()*(CF*_K*_I->_coef[i]);
    evalHelper(f,H,&fB,&HB,*_c);
    return 0.5f*delta.squaredNorm()*_K;
}
scalarD FixPointEnergy::eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const
{
    //energy=0.5*_K*|\sum _I._v[i] _I._b[i] - p|^2
    Vec3 delta=-_pos;
    for(int i=0; i<(_dim+1); i++)
        delta+=U[(*_c)(i)].cast<scalar>()*_I->_coef[i];

    Eigen::Matrix<scalarD,12,12> HB=_HB*CH;
    Eigen::Matrix<scalarD,12,1> fB;
    for(int i=0; i<(_dim+1); i++)
        if(f)fB.block<3,1>(i*3,0)=-delta.cast<scalarD>()*(CF*_K*_I->_coef[i]);
    evalHelper(U,f,H,&fB,&HB,*_c);
    return 0.5f*delta.squaredNorm()*_K;
}
void FixPointEnergy::flagVertex(vector<bool>& flag) const
{
    for(char v=0; v<_dim+1; v++)
        flag[(*_c)(v)]=true;
}
boost::shared_ptr<FEMEnergy> FixPointEnergy::copy(const FEMBody& other) const
{
    return boost::shared_ptr<FEMEnergy>(new FixPointEnergy(_dim,other.getCPtr(_c->_index),*_I,_K,_pos));
}
Vec3 FixPointEnergy::getPos() const
{
    return _pos;
}
void FixPointEnergy::setPos(const Vec3& pos)
{
    _pos=pos;
}
const FEMInterp& FixPointEnergy::getI() const
{
    return *_I;
}
boost::shared_ptr<FEMCell> FixPointEnergy::getCell() const
{
    return _c;
}
//STVK energy
MaterialEnergy::MaterialEnergy() {}
MaterialEnergy::MaterialEnergy(boost::shared_ptr<FEMCell> c,const boost::property_tree::ptree& tree):_c(c)
{
    setParam(tree);
    _dimMask=c->_v[3] == NULL ? 0.0f : 1.0f;
    if(c->_v[3] == NULL)
        calcGComp2D(_GComp,c->_d);
    else calcGComp3D(_GComp,c->_d);
}
scalarD MaterialEnergy::eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const
{
    Eigen::Matrix<scalarD,12,1> fB;
    Eigen::Matrix<scalarD,12,12> HB;
    scalarD E=eval(fB,H?&HB:NULL);
    fB*=CF;
    HB*=CH;
    evalHelper(f,H,&fB,&HB,*_c);
    return E;
}
scalarD MaterialEnergy::eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const
{
    for(char v=0; v<4; v++)
        if(_c->_v[v])_c->_v[v]->_pos=U[(*_c)(v)].cast<scalar>();
    Eigen::Matrix<scalarD,12,1> fB;
    Eigen::Matrix<scalarD,12,12> HB;
    scalarD E=eval(fB,H?&HB:NULL);
    fB*=CF;
    HB*=CH;
    evalHelper(U,f,H,&fB,&HB,*_c);
    return E;
}
scalarD MaterialEnergy::eval(Eigen::Matrix<scalarD,12,1>& fB,Eigen::Matrix<scalarD,12,12>* HB) const
{
    scalarD E;
    Eigen::Matrix<scalarD,9,9> FHB,DRDF;
    Eigen::Matrix<scalarD,9,1> FfB;
    Eigen::Matrix<scalarD,3,4> x;
    x.col(0)=(*_c)[0].cast<scalarD>();
    x.col(1)=(*_c)[1].cast<scalarD>();
    x.col(2)=(*_c)[2].cast<scalarD>();
    if(_c->_v[3])x.col(3)=(*_c)[3].cast<scalarD>();
    else x.col(3).setZero();

    Mat3 F,FN,FR;
    _c->buildF(F,FN,FR,_type == COROTATIONAL_EXACT ? &DRDF : NULL,_invertible);
    scalarD vol=_c->_mass;
    switch(_type) {
    case STVK:
        E=evalF(F,_lambda,_mu,vol,_dimMask,FfB,&FHB);
        break;
    case LINEAR:
        E=evalFLinear(F,_lambda,_mu,vol,_dimMask,FfB,&FHB);
        break;
    case COROTATIONAL:
        E=evalFCorot(F,FR,_lambda,_mu,vol,_dimMask,FfB,&FHB);
        break;
    case COROTATIONAL_EXACT:
        E=evalFCorotE(F,FR,DRDF,_lambda,_mu,vol,_dimMask,FfB,&FHB);
        break;
    case NONHK:
        E=evalFNonHK(F,_lambda,_mu,vol,1.0f-_dimMask,0.01f,FfB,&FHB);
        break;
    case FUNG:
        E=evalFFung(F,_lambda,_mu,_lambdaFung,_muFung,_cFung,vol,_dimMask,FfB,&FHB);
        break;
    default :
        E = 0.0;
        break;
    }
    fB=-FfB.transpose()*_GComp;
    if(HB)*HB=_GComp.transpose()*(FHB*_GComp).eval();
    return E;
}
void MaterialEnergy::flagVertex(vector<bool>& flag) const
{
    for(char v=0; v<4; v++)
        if(_c->_v[v])flag[_c->_v[v]->_index]=true;
}
boost::shared_ptr<FEMEnergy> MaterialEnergy::copy(const FEMBody& other) const
{
    ostringstream oss;
    oss << "cell" << _c->_index;
    boost::property_tree::ptree child=other._tree.get_child(oss.str());
    return boost::shared_ptr<FEMEnergy>(new MaterialEnergy(other.getCPtr(_c->_index),child));
}
void MaterialEnergy::setInvertible(bool invertible)
{
    _invertible=invertible;
}
MaterialEnergy::TYPE MaterialEnergy::getType() const
{
    return _type;
}
void MaterialEnergy::setType(TYPE type)
{
    _type=type;
}
scalar MaterialEnergy::getVol() const
{
    return _c->_mass;
}
scalar MaterialEnergy::getMu() const
{
    return _mu;
}
void MaterialEnergy::setMu(scalar mu)
{
    _mu=mu;
}
scalar MaterialEnergy::getLambda() const
{
    return _lambda;
}
void MaterialEnergy::setLambda(scalar lambda)
{
    _lambda=lambda;
}
scalar MaterialEnergy::getMuFung() const
{
    return _muFung;
}
void MaterialEnergy::setMuFung(scalar muFung)
{
    _muFung=muFung;
}
scalar MaterialEnergy::getLambdaFung() const
{
    return _lambdaFung;
}
void MaterialEnergy::setLambdaFung(scalar lambdaFung)
{
    _lambdaFung=lambdaFung;
}
scalar MaterialEnergy::getCFung() const
{
    return _cFung;
}
void MaterialEnergy::setCFung(scalar cFung)
{
    _cFung=cFung;
}
boost::shared_ptr<FEMCell> MaterialEnergy::getCell() const
{
    return _c;
}
void MaterialEnergy::setParam(const boost::property_tree::ptree& tree)
{
    if(tree.find("young") != tree.not_found() &&
            tree.find("poisson") != tree.not_found()) {
        scalar young=tree.get<scalar>("young");
        scalar poisson=tree.get<scalar>("poisson");
        _mu=young/(2.0f*(1.0f+poisson));
        _lambda=poisson*young/((1.0f+poisson)*(1.0f-2.0f*poisson));
    } else {
        _mu=tree.get<scalar>("mu");
        _lambda=tree.get<scalar>("lambda");
    }
    _muFung=tree.get<scalar>("muFung",_mu*1E-3f);
    _lambdaFung=tree.get<scalar>("lambdaFung",_lambda*1E-3f);
    _cFung=tree.get<scalar>("cFung",1E-4f);
    _type=(TYPE)tree.get<sizeType>("type");
    _invertible=tree.get<bool>("invertible");
}
scalarD MaterialEnergy::evalF(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                              Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    scalarD E;
    //input
    //scalarD V;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD mu;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;
    scalarD tt61;
    scalarD tt62;
    scalarD tt63;
    scalarD tt64;
    scalarD tt65;
    scalarD tt66;
    scalarD tt67;
    scalarD tt68;
    scalarD tt69;
    scalarD tt70;
    scalarD tt71;
    scalarD tt72;
    scalarD tt73;
    scalarD tt74;
    scalarD tt75;
    scalarD tt76;
    scalarD tt77;

    tt1=pow(f(0,0),2);
    tt2=pow(f(1,0),2);
    tt3=pow(f(2,0),2);
    tt4=tt3+tt2+tt1-1;
    tt5=f(2,0)*f(2,1)+f(1,0)*f(1,1)+f(0,0)*f(0,1);
    tt6=pow(f(0,1),2);
    tt7=pow(f(1,1),2);
    tt8=pow(f(2,1),2);
    tt9=tt8+tt7+tt6-1;
    tt10=f(2,0)*f(2,2)+f(1,0)*f(1,2)+f(0,0)*f(0,2);
    tt11=f(2,1)*f(2,2)+f(1,1)*f(1,2)+f(0,1)*f(0,2);
    tt12=pow(f(0,2),2);
    tt13=pow(f(1,2),2);
    tt14=pow(f(2,2),2);
    tt15=-dimMask+tt14+tt13+tt12;
    tt16=0.5*tt15+0.5*tt9+0.5*tt4;
    E=V*(pow(tt16,2)*lambda/2.0+mu*(0.25*pow(tt15,2)+0.5*pow(tt11,2)+0.5*pow(tt10,2)+0.25*pow(tt9,2)+0.5*pow(tt5,2)+0.25*pow(tt4,2)));
    grad[0]=V*(1.0*f(0,0)*tt16*lambda+(1.0*f(0,2)*tt10+1.0*f(0,1)*tt5+1.0*f(0,0)*tt4)*mu);
    grad[1]=V*(1.0*f(1,0)*tt16*lambda+(1.0*f(1,2)*tt10+1.0*f(1,1)*tt5+1.0*f(1,0)*tt4)*mu);
    grad[2]=V*(1.0*f(2,0)*tt16*lambda+(1.0*f(2,2)*tt10+1.0*f(2,1)*tt5+1.0*f(2,0)*tt4)*mu);
    grad[3]=V*(1.0*f(0,1)*tt16*lambda+(1.0*f(0,2)*tt11+1.0*f(0,1)*tt9+1.0*f(0,0)*tt5)*mu);
    grad[4]=V*(1.0*f(1,1)*tt16*lambda+(1.0*f(1,2)*tt11+1.0*f(1,1)*tt9+1.0*f(1,0)*tt5)*mu);
    grad[5]=V*(1.0*f(2,1)*tt16*lambda+(1.0*f(2,2)*tt11+1.0*f(2,1)*tt9+1.0*f(2,0)*tt5)*mu);
    grad[6]=V*(1.0*f(0,2)*tt16*lambda+mu*(1.0*f(0,2)*tt15+1.0*f(0,1)*tt11+1.0*f(0,0)*tt10));
    grad[7]=V*(1.0*f(1,2)*tt16*lambda+mu*(1.0*f(1,2)*tt15+1.0*f(1,1)*tt11+1.0*f(1,0)*tt10));
    grad[8]=V*(1.0*f(2,2)*tt16*lambda+mu*(1.0*f(2,2)*tt15+1.0*f(2,1)*tt11+1.0*f(2,0)*tt10));
    if(!hess)return E;

    tt17=1.0*tt6;
    tt18=1.0*tt12;
    tt19=1.0*tt4;
    tt20=1.0*tt16*lambda;
    tt21=1.0*f(0,1)*f(1,1);
    tt22=1.0*f(0,2)*f(1,2);
    tt23=V*(1.0*f(0,0)*f(1,0)*lambda+(tt22+tt21+2.0*f(0,0)*f(1,0))*mu);
    tt24=1.0*f(0,1)*f(2,1);
    tt25=1.0*f(0,2)*f(2,2);
    tt26=V*(1.0*f(0,0)*f(2,0)*lambda+(tt25+tt24+2.0*f(0,0)*f(2,0))*mu);
    tt27=1.0*tt5;
    tt28=V*(1.0*f(0,0)*f(0,1)*lambda+(tt27+1.0*f(0,0)*f(0,1))*mu);
    tt29=V*(1.0*f(0,0)*f(1,1)*lambda+1.0*f(0,1)*f(1,0)*mu);
    tt30=V*(1.0*f(0,0)*f(2,1)*lambda+1.0*f(0,1)*f(2,0)*mu);
    tt31=1.0*tt10;
    tt32=V*(1.0*f(0,0)*f(0,2)*lambda+(tt31+1.0*f(0,0)*f(0,2))*mu);
    tt33=V*(1.0*f(0,0)*f(1,2)*lambda+1.0*f(0,2)*f(1,0)*mu);
    tt34=V*(1.0*f(0,0)*f(2,2)*lambda+1.0*f(0,2)*f(2,0)*mu);
    tt35=1.0*tt7;
    tt36=1.0*tt13;
    tt37=1.0*f(1,1)*f(2,1);
    tt38=1.0*f(1,2)*f(2,2);
    tt39=V*(1.0*f(1,0)*f(2,0)*lambda+(tt38+tt37+2.0*f(1,0)*f(2,0))*mu);
    tt40=V*(1.0*f(0,1)*f(1,0)*lambda+1.0*f(0,0)*f(1,1)*mu);
    tt41=V*(1.0*f(1,0)*f(1,1)*lambda+(tt27+1.0*f(1,0)*f(1,1))*mu);
    tt42=V*(1.0*f(1,0)*f(2,1)*lambda+1.0*f(1,1)*f(2,0)*mu);
    tt43=V*(1.0*f(0,2)*f(1,0)*lambda+1.0*f(0,0)*f(1,2)*mu);
    tt44=V*(1.0*f(1,0)*f(1,2)*lambda+(tt31+1.0*f(1,0)*f(1,2))*mu);
    tt45=V*(1.0*f(1,0)*f(2,2)*lambda+1.0*f(1,2)*f(2,0)*mu);
    tt46=1.0*tt8;
    tt47=1.0*tt14;
    tt48=V*(1.0*f(0,1)*f(2,0)*lambda+1.0*f(0,0)*f(2,1)*mu);
    tt49=V*(1.0*f(1,1)*f(2,0)*lambda+1.0*f(1,0)*f(2,1)*mu);
    tt50=V*(1.0*f(2,0)*f(2,1)*lambda+(tt27+1.0*f(2,0)*f(2,1))*mu);
    tt51=V*(1.0*f(0,2)*f(2,0)*lambda+1.0*f(0,0)*f(2,2)*mu);
    tt52=V*(1.0*f(1,2)*f(2,0)*lambda+1.0*f(1,0)*f(2,2)*mu);
    tt53=V*(1.0*f(2,0)*f(2,2)*lambda+(tt31+1.0*f(2,0)*f(2,2))*mu);
    tt54=1.0*tt1;
    tt55=1.0*tt9;
    tt56=1.0*f(0,0)*f(1,0);
    tt57=V*(1.0*f(0,1)*f(1,1)*lambda+(tt22+2.0*f(0,1)*f(1,1)+tt56)*mu);
    tt58=1.0*f(0,0)*f(2,0);
    tt59=V*(1.0*f(0,1)*f(2,1)*lambda+(tt25+2.0*f(0,1)*f(2,1)+tt58)*mu);
    tt60=1.0*tt11;
    tt61=V*(1.0*f(0,1)*f(0,2)*lambda+(tt60+1.0*f(0,1)*f(0,2))*mu);
    tt62=V*(1.0*f(0,1)*f(1,2)*lambda+1.0*f(0,2)*f(1,1)*mu);
    tt63=V*(1.0*f(0,1)*f(2,2)*lambda+1.0*f(0,2)*f(2,1)*mu);
    tt64=1.0*tt2;
    tt65=1.0*f(1,0)*f(2,0);
    tt66=V*(1.0*f(1,1)*f(2,1)*lambda+(tt38+2.0*f(1,1)*f(2,1)+tt65)*mu);
    tt67=V*(1.0*f(0,2)*f(1,1)*lambda+1.0*f(0,1)*f(1,2)*mu);
    tt68=V*(1.0*f(1,1)*f(1,2)*lambda+(tt60+1.0*f(1,1)*f(1,2))*mu);
    tt69=V*(1.0*f(1,1)*f(2,2)*lambda+1.0*f(1,2)*f(2,1)*mu);
    tt70=1.0*tt3;
    tt71=V*(1.0*f(0,2)*f(2,1)*lambda+1.0*f(0,1)*f(2,2)*mu);
    tt72=V*(1.0*f(1,2)*f(2,1)*lambda+1.0*f(1,1)*f(2,2)*mu);
    tt73=V*(1.0*f(2,1)*f(2,2)*lambda+(tt60+1.0*f(2,1)*f(2,2))*mu);
    tt74=1.0*tt15;
    tt75=V*(1.0*f(0,2)*f(1,2)*lambda+(2.0*f(0,2)*f(1,2)+tt21+tt56)*mu);
    tt76=V*(1.0*f(0,2)*f(2,2)*lambda+(2.0*f(0,2)*f(2,2)+tt24+tt58)*mu);
    tt77=V*(1.0*f(1,2)*f(2,2)*lambda+(2.0*f(1,2)*f(2,2)+tt37+tt65)*mu);
    (*hess)(0,0)=V*(tt20+1.0*tt1*lambda+(tt19+tt18+tt17+2.0*tt1)*mu);
    (*hess)(0,1)=tt23;
    (*hess)(0,2)=tt26;
    (*hess)(0,3)=tt28;
    (*hess)(0,4)=tt29;
    (*hess)(0,5)=tt30;
    (*hess)(0,6)=tt32;
    (*hess)(0,7)=tt33;
    (*hess)(0,8)=tt34;
    (*hess)(1,0)=tt23;
    (*hess)(1,1)=V*(tt20+1.0*tt2*lambda+(tt19+tt36+tt35+2.0*tt2)*mu);
    (*hess)(1,2)=tt39;
    (*hess)(1,3)=tt40;
    (*hess)(1,4)=tt41;
    (*hess)(1,5)=tt42;
    (*hess)(1,6)=tt43;
    (*hess)(1,7)=tt44;
    (*hess)(1,8)=tt45;
    (*hess)(2,0)=tt26;
    (*hess)(2,1)=tt39;
    (*hess)(2,2)=V*(tt20+1.0*tt3*lambda+(tt47+tt46+tt19+2.0*tt3)*mu);
    (*hess)(2,3)=tt48;
    (*hess)(2,4)=tt49;
    (*hess)(2,5)=tt50;
    (*hess)(2,6)=tt51;
    (*hess)(2,7)=tt52;
    (*hess)(2,8)=tt53;
    (*hess)(3,0)=tt28;
    (*hess)(3,1)=tt40;
    (*hess)(3,2)=tt48;
    (*hess)(3,3)=V*(tt20+1.0*tt6*lambda+(tt55+tt18+2.0*tt6+tt54)*mu);
    (*hess)(3,4)=tt57;
    (*hess)(3,5)=tt59;
    (*hess)(3,6)=tt61;
    (*hess)(3,7)=tt62;
    (*hess)(3,8)=tt63;
    (*hess)(4,0)=tt29;
    (*hess)(4,1)=tt41;
    (*hess)(4,2)=tt49;
    (*hess)(4,3)=tt57;
    (*hess)(4,4)=V*(tt20+1.0*tt7*lambda+(tt55+tt36+2.0*tt7+tt64)*mu);
    (*hess)(4,5)=tt66;
    (*hess)(4,6)=tt67;
    (*hess)(4,7)=tt68;
    (*hess)(4,8)=tt69;
    (*hess)(5,0)=tt30;
    (*hess)(5,1)=tt42;
    (*hess)(5,2)=tt50;
    (*hess)(5,3)=tt59;
    (*hess)(5,4)=tt66;
    (*hess)(5,5)=V*(tt20+1.0*tt8*lambda+(tt47+tt55+2.0*tt8+tt70)*mu);
    (*hess)(5,6)=tt71;
    (*hess)(5,7)=tt72;
    (*hess)(5,8)=tt73;
    (*hess)(6,0)=tt32;
    (*hess)(6,1)=tt43;
    (*hess)(6,2)=tt51;
    (*hess)(6,3)=tt61;
    (*hess)(6,4)=tt67;
    (*hess)(6,5)=tt71;
    (*hess)(6,6)=V*(tt20+1.0*tt12*lambda+mu*(tt74+2.0*tt12+tt17+tt54));
    (*hess)(6,7)=tt75;
    (*hess)(6,8)=tt76;
    (*hess)(7,0)=tt33;
    (*hess)(7,1)=tt44;
    (*hess)(7,2)=tt52;
    (*hess)(7,3)=tt62;
    (*hess)(7,4)=tt68;
    (*hess)(7,5)=tt72;
    (*hess)(7,6)=tt75;
    (*hess)(7,7)=V*(tt20+1.0*tt13*lambda+mu*(tt74+2.0*tt13+tt35+tt64));
    (*hess)(7,8)=tt77;
    (*hess)(8,0)=tt34;
    (*hess)(8,1)=tt45;
    (*hess)(8,2)=tt53;
    (*hess)(8,3)=tt63;
    (*hess)(8,4)=tt69;
    (*hess)(8,5)=tt73;
    (*hess)(8,6)=tt76;
    (*hess)(8,7)=tt77;
    (*hess)(8,8)=V*(tt20+1.0*tt14*lambda+mu*(tt74+2.0*tt14+tt46+tt70));
    return E;
}
scalarD MaterialEnergy::evalFLinear(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                                    Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    scalarD E;
    //input
    //scalarD V;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD mu;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;

    tt1=1.0*f(0,0);
    tt2=tt1-1;
    tt3=f(1,0)+f(0,1);
    tt4=1.0*f(1,1);
    tt5=tt4-1;
    tt6=f(2,0)+f(0,2);
    tt7=f(2,1)+f(1,2);
    tt8=1.0*f(2,2);
    tt9=-dimMask;
    tt10=tt9+tt8;
    tt11=tt9+tt8+tt4+tt1-2;
    tt12=1.0*tt11*lambda;
    tt13=1.0*tt3*mu*V;
    tt14=1.0*tt6*mu*V;
    tt15=1.0*tt7*mu*V;
    E=V*(pow(tt11,2)*lambda/2.0+mu*(pow(tt10,2)+0.5*pow(tt7,2)+0.5*pow(tt6,2)+pow(tt5,2)+0.5*pow(tt3,2)+pow(tt2,2)));
    grad[0]=V*(tt12+2.0*tt2*mu);
    grad[1]=tt13;
    grad[2]=tt14;
    grad[3]=tt13;
    grad[4]=V*(tt12+2.0*tt5*mu);
    grad[5]=tt15;
    grad[6]=tt14;
    grad[7]=tt15;
    grad[8]=V*(tt12+2.0*mu*tt10);
    if(!hess)return E;

    tt16=V*(1.0*lambda+2.0*mu);
    tt17=1.0*V*lambda;
    tt18=1.0*mu*V;
    (*hess)(0,0)=tt16;
    (*hess)(0,1)=0;
    (*hess)(0,2)=0;
    (*hess)(0,3)=0;
    (*hess)(0,4)=tt17;
    (*hess)(0,5)=0;
    (*hess)(0,6)=0;
    (*hess)(0,7)=0;
    (*hess)(0,8)=tt17;
    (*hess)(1,0)=0;
    (*hess)(1,1)=tt18;
    (*hess)(1,2)=0;
    (*hess)(1,3)=tt18;
    (*hess)(1,4)=0;
    (*hess)(1,5)=0;
    (*hess)(1,6)=0;
    (*hess)(1,7)=0;
    (*hess)(1,8)=0;
    (*hess)(2,0)=0;
    (*hess)(2,1)=0;
    (*hess)(2,2)=tt18;
    (*hess)(2,3)=0;
    (*hess)(2,4)=0;
    (*hess)(2,5)=0;
    (*hess)(2,6)=tt18;
    (*hess)(2,7)=0;
    (*hess)(2,8)=0;
    (*hess)(3,0)=0;
    (*hess)(3,1)=tt18;
    (*hess)(3,2)=0;
    (*hess)(3,3)=tt18;
    (*hess)(3,4)=0;
    (*hess)(3,5)=0;
    (*hess)(3,6)=0;
    (*hess)(3,7)=0;
    (*hess)(3,8)=0;
    (*hess)(4,0)=tt17;
    (*hess)(4,1)=0;
    (*hess)(4,2)=0;
    (*hess)(4,3)=0;
    (*hess)(4,4)=tt16;
    (*hess)(4,5)=0;
    (*hess)(4,6)=0;
    (*hess)(4,7)=0;
    (*hess)(4,8)=tt17;
    (*hess)(5,0)=0;
    (*hess)(5,1)=0;
    (*hess)(5,2)=0;
    (*hess)(5,3)=0;
    (*hess)(5,4)=0;
    (*hess)(5,5)=tt18;
    (*hess)(5,6)=0;
    (*hess)(5,7)=tt18;
    (*hess)(5,8)=0;
    (*hess)(6,0)=0;
    (*hess)(6,1)=0;
    (*hess)(6,2)=tt18;
    (*hess)(6,3)=0;
    (*hess)(6,4)=0;
    (*hess)(6,5)=0;
    (*hess)(6,6)=tt18;
    (*hess)(6,7)=0;
    (*hess)(6,8)=0;
    (*hess)(7,0)=0;
    (*hess)(7,1)=0;
    (*hess)(7,2)=0;
    (*hess)(7,3)=0;
    (*hess)(7,4)=0;
    (*hess)(7,5)=tt18;
    (*hess)(7,6)=0;
    (*hess)(7,7)=tt18;
    (*hess)(7,8)=0;
    (*hess)(8,0)=tt17;
    (*hess)(8,1)=0;
    (*hess)(8,2)=0;
    (*hess)(8,3)=0;
    (*hess)(8,4)=tt17;
    (*hess)(8,5)=0;
    (*hess)(8,6)=0;
    (*hess)(8,7)=0;
    (*hess)(8,8)=tt16;
    return E;
}
scalarD MaterialEnergy::evalFCorot(const Mat3& f,const Mat3& r,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                                   Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    scalarD E;
    //input
    //scalarD V;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD mu;
    //Mat3d r;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;
    scalarD tt61;
    scalarD tt62;
    scalarD tt63;
    scalarD tt64;
    scalarD tt65;
    scalarD tt66;
    scalarD tt67;
    scalarD tt68;
    scalarD tt69;
    scalarD tt70;
    scalarD tt71;
    scalarD tt72;
    scalarD tt73;
    scalarD tt74;

    tt1=0.5*(2*f(2,0)*r(2,0)+2*f(1,0)*r(1,0)+2*f(0,0)*r(0,0));
    tt2=tt1-1;
    tt3=f(2,0)*r(2,1)+r(2,0)*f(2,1)+f(1,0)*r(1,1)+r(1,0)*f(1,1)+f(0,0)*r(0,1)+r(0,0)*f(0,1);
    tt4=0.5*(2*f(2,1)*r(2,1)+2*f(1,1)*r(1,1)+2*f(0,1)*r(0,1));
    tt5=tt4-1;
    tt6=f(2,0)*r(2,2)+r(2,0)*f(2,2)+f(1,0)*r(1,2)+r(1,0)*f(1,2)+f(0,0)*r(0,2)+r(0,0)*f(0,2);
    tt7=f(2,1)*r(2,2)+r(2,1)*f(2,2)+f(1,1)*r(1,2)+r(1,1)*f(1,2)+f(0,1)*r(0,2)+r(0,1)*f(0,2);
    tt8=0.5*(2*f(2,2)*r(2,2)+2*f(1,2)*r(1,2)+2*f(0,2)*r(0,2));
    tt9=-dimMask;
    tt10=tt9+tt8;
    tt11=tt9+tt8+tt4+tt1-2;
    E=V*(pow(tt11,2)*lambda/2.0+mu*(pow(tt10,2)+0.5*pow(tt7,2)+0.5*pow(tt6,2)+pow(tt5,2)+0.5*pow(tt3,2)+pow(tt2,2)));
    grad[0]=V*(1.0*r(0,0)*tt11*lambda+(1.0*r(0,2)*tt6+1.0*r(0,1)*tt3+2.0*r(0,0)*tt2)*mu);
    grad[1]=V*(1.0*r(1,0)*tt11*lambda+(1.0*r(1,2)*tt6+1.0*r(1,1)*tt3+2.0*r(1,0)*tt2)*mu);
    grad[2]=V*(1.0*r(2,0)*tt11*lambda+(1.0*r(2,2)*tt6+1.0*r(2,1)*tt3+2.0*r(2,0)*tt2)*mu);
    grad[3]=V*(1.0*r(0,1)*tt11*lambda+(1.0*r(0,2)*tt7+2.0*r(0,1)*tt5+1.0*r(0,0)*tt3)*mu);
    grad[4]=V*(1.0*r(1,1)*tt11*lambda+(1.0*r(1,2)*tt7+2.0*r(1,1)*tt5+1.0*r(1,0)*tt3)*mu);
    grad[5]=V*(1.0*r(2,1)*tt11*lambda+(1.0*r(2,2)*tt7+2.0*r(2,1)*tt5+1.0*r(2,0)*tt3)*mu);
    grad[6]=V*(1.0*r(0,2)*tt11*lambda+mu*(2.0*r(0,2)*tt10+1.0*r(0,1)*tt7+1.0*r(0,0)*tt6));
    grad[7]=V*(1.0*r(1,2)*tt11*lambda+mu*(2.0*r(1,2)*tt10+1.0*r(1,1)*tt7+1.0*r(1,0)*tt6));
    grad[8]=V*(1.0*r(2,2)*tt11*lambda+mu*(2.0*r(2,2)*tt10+1.0*r(2,1)*tt7+1.0*r(2,0)*tt6));
    if(!hess)return E;

    tt12=pow(r(0,0),2);
    tt13=pow(r(0,1),2);
    tt14=1.0*tt13;
    tt15=pow(r(0,2),2);
    tt16=1.0*tt15;
    tt17=1.0*r(0,1)*r(1,1);
    tt18=1.0*r(0,2)*r(1,2);
    tt19=V*(1.0*r(0,0)*r(1,0)*lambda+(tt18+tt17+2.0*r(0,0)*r(1,0))*mu);
    tt20=1.0*r(0,1)*r(2,1);
    tt21=1.0*r(0,2)*r(2,2);
    tt22=V*(1.0*r(0,0)*r(2,0)*lambda+(tt21+tt20+2.0*r(0,0)*r(2,0))*mu);
    tt23=V*(1.0*r(0,0)*r(0,1)*lambda+1.0*r(0,0)*r(0,1)*mu);
    tt24=V*(1.0*r(0,0)*r(1,1)*lambda+1.0*r(0,1)*r(1,0)*mu);
    tt25=V*(1.0*r(0,0)*r(2,1)*lambda+1.0*r(0,1)*r(2,0)*mu);
    tt26=V*(1.0*r(0,0)*r(0,2)*lambda+1.0*r(0,0)*r(0,2)*mu);
    tt27=V*(1.0*r(0,0)*r(1,2)*lambda+1.0*r(0,2)*r(1,0)*mu);
    tt28=V*(1.0*r(0,0)*r(2,2)*lambda+1.0*r(0,2)*r(2,0)*mu);
    tt29=pow(r(1,0),2);
    tt30=pow(r(1,1),2);
    tt31=1.0*tt30;
    tt32=pow(r(1,2),2);
    tt33=1.0*tt32;
    tt34=1.0*r(1,1)*r(2,1);
    tt35=1.0*r(1,2)*r(2,2);
    tt36=V*(1.0*r(1,0)*r(2,0)*lambda+(tt35+tt34+2.0*r(1,0)*r(2,0))*mu);
    tt37=V*(1.0*r(0,1)*r(1,0)*lambda+1.0*r(0,0)*r(1,1)*mu);
    tt38=V*(1.0*r(1,0)*r(1,1)*lambda+1.0*r(1,0)*r(1,1)*mu);
    tt39=V*(1.0*r(1,0)*r(2,1)*lambda+1.0*r(1,1)*r(2,0)*mu);
    tt40=V*(1.0*r(0,2)*r(1,0)*lambda+1.0*r(0,0)*r(1,2)*mu);
    tt41=V*(1.0*r(1,0)*r(1,2)*lambda+1.0*r(1,0)*r(1,2)*mu);
    tt42=V*(1.0*r(1,0)*r(2,2)*lambda+1.0*r(1,2)*r(2,0)*mu);
    tt43=pow(r(2,0),2);
    tt44=pow(r(2,1),2);
    tt45=1.0*tt44;
    tt46=pow(r(2,2),2);
    tt47=1.0*tt46;
    tt48=V*(1.0*r(0,1)*r(2,0)*lambda+1.0*r(0,0)*r(2,1)*mu);
    tt49=V*(1.0*r(1,1)*r(2,0)*lambda+1.0*r(1,0)*r(2,1)*mu);
    tt50=V*(1.0*r(2,0)*r(2,1)*lambda+1.0*r(2,0)*r(2,1)*mu);
    tt51=V*(1.0*r(0,2)*r(2,0)*lambda+1.0*r(0,0)*r(2,2)*mu);
    tt52=V*(1.0*r(1,2)*r(2,0)*lambda+1.0*r(1,0)*r(2,2)*mu);
    tt53=V*(1.0*r(2,0)*r(2,2)*lambda+1.0*r(2,0)*r(2,2)*mu);
    tt54=1.0*tt12;
    tt55=1.0*r(0,0)*r(1,0);
    tt56=V*(1.0*r(0,1)*r(1,1)*lambda+(tt18+2.0*r(0,1)*r(1,1)+tt55)*mu);
    tt57=1.0*r(0,0)*r(2,0);
    tt58=V*(1.0*r(0,1)*r(2,1)*lambda+(tt21+2.0*r(0,1)*r(2,1)+tt57)*mu);
    tt59=V*(1.0*r(0,1)*r(0,2)*lambda+1.0*r(0,1)*r(0,2)*mu);
    tt60=V*(1.0*r(0,1)*r(1,2)*lambda+1.0*r(0,2)*r(1,1)*mu);
    tt61=V*(1.0*r(0,1)*r(2,2)*lambda+1.0*r(0,2)*r(2,1)*mu);
    tt62=1.0*tt29;
    tt63=1.0*r(1,0)*r(2,0);
    tt64=V*(1.0*r(1,1)*r(2,1)*lambda+(tt35+2.0*r(1,1)*r(2,1)+tt63)*mu);
    tt65=V*(1.0*r(0,2)*r(1,1)*lambda+1.0*r(0,1)*r(1,2)*mu);
    tt66=V*(1.0*r(1,1)*r(1,2)*lambda+1.0*r(1,1)*r(1,2)*mu);
    tt67=V*(1.0*r(1,1)*r(2,2)*lambda+1.0*r(1,2)*r(2,1)*mu);
    tt68=1.0*tt43;
    tt69=V*(1.0*r(0,2)*r(2,1)*lambda+1.0*r(0,1)*r(2,2)*mu);
    tt70=V*(1.0*r(1,2)*r(2,1)*lambda+1.0*r(1,1)*r(2,2)*mu);
    tt71=V*(1.0*r(2,1)*r(2,2)*lambda+1.0*r(2,1)*r(2,2)*mu);
    tt72=V*(1.0*r(0,2)*r(1,2)*lambda+(2.0*r(0,2)*r(1,2)+tt17+tt55)*mu);
    tt73=V*(1.0*r(0,2)*r(2,2)*lambda+(2.0*r(0,2)*r(2,2)+tt20+tt57)*mu);
    tt74=V*(1.0*r(1,2)*r(2,2)*lambda+(2.0*r(1,2)*r(2,2)+tt34+tt63)*mu);
    (*hess)(0,0)=V*(1.0*tt12*lambda+(tt16+tt14+2.0*tt12)*mu);
    (*hess)(0,1)=tt19;
    (*hess)(0,2)=tt22;
    (*hess)(0,3)=tt23;
    (*hess)(0,4)=tt24;
    (*hess)(0,5)=tt25;
    (*hess)(0,6)=tt26;
    (*hess)(0,7)=tt27;
    (*hess)(0,8)=tt28;
    (*hess)(1,0)=tt19;
    (*hess)(1,1)=V*(1.0*tt29*lambda+(tt33+tt31+2.0*tt29)*mu);
    (*hess)(1,2)=tt36;
    (*hess)(1,3)=tt37;
    (*hess)(1,4)=tt38;
    (*hess)(1,5)=tt39;
    (*hess)(1,6)=tt40;
    (*hess)(1,7)=tt41;
    (*hess)(1,8)=tt42;
    (*hess)(2,0)=tt22;
    (*hess)(2,1)=tt36;
    (*hess)(2,2)=V*(1.0*tt43*lambda+(tt47+tt45+2.0*tt43)*mu);
    (*hess)(2,3)=tt48;
    (*hess)(2,4)=tt49;
    (*hess)(2,5)=tt50;
    (*hess)(2,6)=tt51;
    (*hess)(2,7)=tt52;
    (*hess)(2,8)=tt53;
    (*hess)(3,0)=tt23;
    (*hess)(3,1)=tt37;
    (*hess)(3,2)=tt48;
    (*hess)(3,3)=V*(1.0*tt13*lambda+(tt16+2.0*tt13+tt54)*mu);
    (*hess)(3,4)=tt56;
    (*hess)(3,5)=tt58;
    (*hess)(3,6)=tt59;
    (*hess)(3,7)=tt60;
    (*hess)(3,8)=tt61;
    (*hess)(4,0)=tt24;
    (*hess)(4,1)=tt38;
    (*hess)(4,2)=tt49;
    (*hess)(4,3)=tt56;
    (*hess)(4,4)=V*(1.0*tt30*lambda+(tt33+2.0*tt30+tt62)*mu);
    (*hess)(4,5)=tt64;
    (*hess)(4,6)=tt65;
    (*hess)(4,7)=tt66;
    (*hess)(4,8)=tt67;
    (*hess)(5,0)=tt25;
    (*hess)(5,1)=tt39;
    (*hess)(5,2)=tt50;
    (*hess)(5,3)=tt58;
    (*hess)(5,4)=tt64;
    (*hess)(5,5)=V*(1.0*tt44*lambda+(tt47+2.0*tt44+tt68)*mu);
    (*hess)(5,6)=tt69;
    (*hess)(5,7)=tt70;
    (*hess)(5,8)=tt71;
    (*hess)(6,0)=tt26;
    (*hess)(6,1)=tt40;
    (*hess)(6,2)=tt51;
    (*hess)(6,3)=tt59;
    (*hess)(6,4)=tt65;
    (*hess)(6,5)=tt69;
    (*hess)(6,6)=V*(1.0*tt15*lambda+(2.0*tt15+tt14+tt54)*mu);
    (*hess)(6,7)=tt72;
    (*hess)(6,8)=tt73;
    (*hess)(7,0)=tt27;
    (*hess)(7,1)=tt41;
    (*hess)(7,2)=tt52;
    (*hess)(7,3)=tt60;
    (*hess)(7,4)=tt66;
    (*hess)(7,5)=tt70;
    (*hess)(7,6)=tt72;
    (*hess)(7,7)=V*(1.0*tt32*lambda+(2.0*tt32+tt31+tt62)*mu);
    (*hess)(7,8)=tt74;
    (*hess)(8,0)=tt28;
    (*hess)(8,1)=tt42;
    (*hess)(8,2)=tt53;
    (*hess)(8,3)=tt61;
    (*hess)(8,4)=tt67;
    (*hess)(8,5)=tt71;
    (*hess)(8,6)=tt73;
    (*hess)(8,7)=tt74;
    (*hess)(8,8)=V*(1.0*tt46*lambda+(2.0*tt46+tt45+tt68)*mu);
    return E;
}
scalarD MaterialEnergy::evalFCorotE(const Mat3& f,const Mat3& r,const Eigen::Matrix<scalarD,9,9>& DRDF,
                                    scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                                    Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    Eigen::Matrix<scalarD,9,18> hessTmp;
    scalarD E=evalFCorotE(f,r,lambda,mu,V,dimMask,grad,hess?&hessTmp:NULL);
    if(hess) {
        *hess=hessTmp.block<9,9>(0,0);
        *hess+=hessTmp.block<9,9>(0,9)*DRDF;
    }
    return E;
}
scalarD MaterialEnergy::evalFCorotE(const Mat3& f,const Mat3& r,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                                    Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,18>* hess)
{
    //input
    scalarD E;
    //scalarD V;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD mu;
    //Mat3d r;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;
    scalarD tt61;
    scalarD tt62;
    scalarD tt63;
    scalarD tt64;
    scalarD tt65;
    scalarD tt66;
    scalarD tt67;
    scalarD tt68;
    scalarD tt69;
    scalarD tt70;
    scalarD tt71;
    scalarD tt72;
    scalarD tt73;
    scalarD tt74;
    scalarD tt75;
    scalarD tt76;
    scalarD tt77;
    scalarD tt78;
    scalarD tt79;
    scalarD tt80;
    scalarD tt81;
    scalarD tt82;
    scalarD tt83;
    scalarD tt84;
    scalarD tt85;
    scalarD tt86;
    scalarD tt87;
    scalarD tt88;
    scalarD tt89;
    scalarD tt90;
    scalarD tt91;
    scalarD tt92;
    scalarD tt93;
    scalarD tt94;
    scalarD tt95;
    scalarD tt96;
    scalarD tt97;
    scalarD tt98;
    scalarD tt99;
    scalarD tt100;
    scalarD tt101;
    scalarD tt102;
    scalarD tt103;
    scalarD tt104;
    scalarD tt105;
    scalarD tt106;
    scalarD tt107;
    scalarD tt108;

    tt1=0.5*(2*f(2,0)*r(2,0)+2*f(1,0)*r(1,0)+2*f(0,0)*r(0,0));
    tt2=tt1-1;
    tt3=f(2,0)*r(2,1)+r(2,0)*f(2,1)+f(1,0)*r(1,1)+r(1,0)*f(1,1)+f(0,0)*r(0,1)+r(0,0)*f(0,1);
    tt4=0.5*(2*f(2,1)*r(2,1)+2*f(1,1)*r(1,1)+2*f(0,1)*r(0,1));
    tt5=tt4-1;
    tt6=f(2,0)*r(2,2)+r(2,0)*f(2,2)+f(1,0)*r(1,2)+r(1,0)*f(1,2)+f(0,0)*r(0,2)+r(0,0)*f(0,2);
    tt7=f(2,1)*r(2,2)+r(2,1)*f(2,2)+f(1,1)*r(1,2)+r(1,1)*f(1,2)+f(0,1)*r(0,2)+r(0,1)*f(0,2);
    tt8=0.5*(2*f(2,2)*r(2,2)+2*f(1,2)*r(1,2)+2*f(0,2)*r(0,2));
    tt9=-dimMask;
    tt10=tt9+tt8;
    tt11=tt9+tt8+tt4+tt1-2;
    E=V*(pow(tt11,2)*lambda/2.0+mu*(pow(tt10,2)+0.5*pow(tt7,2)+0.5*pow(tt6,2)+pow(tt5,2)+0.5*pow(tt3,2)+pow(tt2,2)));
    grad[0]=V*(1.0*r(0,0)*tt11*lambda+(1.0*r(0,2)*tt6+1.0*r(0,1)*tt3+2.0*r(0,0)*tt2)*mu);
    grad[1]=V*(1.0*r(1,0)*tt11*lambda+(1.0*r(1,2)*tt6+1.0*r(1,1)*tt3+2.0*r(1,0)*tt2)*mu);
    grad[2]=V*(1.0*r(2,0)*tt11*lambda+(1.0*r(2,2)*tt6+1.0*r(2,1)*tt3+2.0*r(2,0)*tt2)*mu);
    grad[3]=V*(1.0*r(0,1)*tt11*lambda+(1.0*r(0,2)*tt7+2.0*r(0,1)*tt5+1.0*r(0,0)*tt3)*mu);
    grad[4]=V*(1.0*r(1,1)*tt11*lambda+(1.0*r(1,2)*tt7+2.0*r(1,1)*tt5+1.0*r(1,0)*tt3)*mu);
    grad[5]=V*(1.0*r(2,1)*tt11*lambda+(1.0*r(2,2)*tt7+2.0*r(2,1)*tt5+1.0*r(2,0)*tt3)*mu);
    grad[6]=V*(1.0*r(0,2)*tt11*lambda+mu*(2.0*r(0,2)*tt10+1.0*r(0,1)*tt7+1.0*r(0,0)*tt6));
    grad[7]=V*(1.0*r(1,2)*tt11*lambda+mu*(2.0*r(1,2)*tt10+1.0*r(1,1)*tt7+1.0*r(1,0)*tt6));
    grad[8]=V*(1.0*r(2,2)*tt11*lambda+mu*(2.0*r(2,2)*tt10+1.0*r(2,1)*tt7+1.0*r(2,0)*tt6));
    if(!hess)return E;

    tt12=pow(r(0,0),2);
    tt13=pow(r(0,1),2);
    tt14=1.0*tt13;
    tt15=pow(r(0,2),2);
    tt16=1.0*tt15;
    tt17=1.0*r(0,1)*r(1,1);
    tt18=1.0*r(0,2)*r(1,2);
    tt19=V*(1.0*r(0,0)*r(1,0)*lambda+(tt18+tt17+2.0*r(0,0)*r(1,0))*mu);
    tt20=1.0*r(0,1)*r(2,1);
    tt21=1.0*r(0,2)*r(2,2);
    tt22=V*(1.0*r(0,0)*r(2,0)*lambda+(tt21+tt20+2.0*r(0,0)*r(2,0))*mu);
    tt23=V*(1.0*r(0,0)*r(0,1)*lambda+1.0*r(0,0)*r(0,1)*mu);
    tt24=V*(1.0*r(0,0)*r(1,1)*lambda+1.0*r(0,1)*r(1,0)*mu);
    tt25=V*(1.0*r(0,0)*r(2,1)*lambda+1.0*r(0,1)*r(2,0)*mu);
    tt26=V*(1.0*r(0,0)*r(0,2)*lambda+1.0*r(0,0)*r(0,2)*mu);
    tt27=V*(1.0*r(0,0)*r(1,2)*lambda+1.0*r(0,2)*r(1,0)*mu);
    tt28=V*(1.0*r(0,0)*r(2,2)*lambda+1.0*r(0,2)*r(2,0)*mu);
    tt29=1.0*f(0,1)*r(0,1);
    tt30=1.0*f(0,2)*r(0,2);
    tt31=2.0*tt2;
    tt32=1.0*tt11*lambda;
    tt33=1.0*r(0,1)*f(1,1);
    tt34=1.0*r(0,2)*f(1,2);
    tt35=1.0*r(0,1)*f(2,1);
    tt36=1.0*r(0,2)*f(2,2);
    tt37=1.0*tt3;
    tt38=1.0*tt6;
    tt39=pow(r(1,0),2);
    tt40=pow(r(1,1),2);
    tt41=1.0*tt40;
    tt42=pow(r(1,2),2);
    tt43=1.0*tt42;
    tt44=1.0*r(1,1)*r(2,1);
    tt45=1.0*r(1,2)*r(2,2);
    tt46=V*(1.0*r(1,0)*r(2,0)*lambda+(tt45+tt44+2.0*r(1,0)*r(2,0))*mu);
    tt47=V*(1.0*r(0,1)*r(1,0)*lambda+1.0*r(0,0)*r(1,1)*mu);
    tt48=V*(1.0*r(1,0)*r(1,1)*lambda+1.0*r(1,0)*r(1,1)*mu);
    tt49=V*(1.0*r(1,0)*r(2,1)*lambda+1.0*r(1,1)*r(2,0)*mu);
    tt50=V*(1.0*r(0,2)*r(1,0)*lambda+1.0*r(0,0)*r(1,2)*mu);
    tt51=V*(1.0*r(1,0)*r(1,2)*lambda+1.0*r(1,0)*r(1,2)*mu);
    tt52=V*(1.0*r(1,0)*r(2,2)*lambda+1.0*r(1,2)*r(2,0)*mu);
    tt53=1.0*f(0,1)*r(1,1);
    tt54=1.0*f(0,2)*r(1,2);
    tt55=1.0*f(1,1)*r(1,1);
    tt56=1.0*f(1,2)*r(1,2);
    tt57=1.0*r(1,1)*f(2,1);
    tt58=1.0*r(1,2)*f(2,2);
    tt59=pow(r(2,0),2);
    tt60=pow(r(2,1),2);
    tt61=1.0*tt60;
    tt62=pow(r(2,2),2);
    tt63=1.0*tt62;
    tt64=V*(1.0*r(0,1)*r(2,0)*lambda+1.0*r(0,0)*r(2,1)*mu);
    tt65=V*(1.0*r(1,1)*r(2,0)*lambda+1.0*r(1,0)*r(2,1)*mu);
    tt66=V*(1.0*r(2,0)*r(2,1)*lambda+1.0*r(2,0)*r(2,1)*mu);
    tt67=V*(1.0*r(0,2)*r(2,0)*lambda+1.0*r(0,0)*r(2,2)*mu);
    tt68=V*(1.0*r(1,2)*r(2,0)*lambda+1.0*r(1,0)*r(2,2)*mu);
    tt69=V*(1.0*r(2,0)*r(2,2)*lambda+1.0*r(2,0)*r(2,2)*mu);
    tt70=1.0*f(0,1)*r(2,1);
    tt71=1.0*f(0,2)*r(2,2);
    tt72=1.0*f(1,1)*r(2,1);
    tt73=1.0*f(1,2)*r(2,2);
    tt74=1.0*f(2,1)*r(2,1);
    tt75=1.0*f(2,2)*r(2,2);
    tt76=1.0*tt12;
    tt77=1.0*r(0,0)*r(1,0);
    tt78=V*(1.0*r(0,1)*r(1,1)*lambda+(tt18+2.0*r(0,1)*r(1,1)+tt77)*mu);
    tt79=1.0*r(0,0)*r(2,0);
    tt80=V*(1.0*r(0,1)*r(2,1)*lambda+(tt21+2.0*r(0,1)*r(2,1)+tt79)*mu);
    tt81=V*(1.0*r(0,1)*r(0,2)*lambda+1.0*r(0,1)*r(0,2)*mu);
    tt82=V*(1.0*r(0,1)*r(1,2)*lambda+1.0*r(0,2)*r(1,1)*mu);
    tt83=V*(1.0*r(0,1)*r(2,2)*lambda+1.0*r(0,2)*r(2,1)*mu);
    tt84=1.0*f(0,0)*r(0,0);
    tt85=2.0*tt5;
    tt86=1.0*r(0,0)*f(1,0);
    tt87=1.0*r(0,0)*f(2,0);
    tt88=1.0*tt7;
    tt89=1.0*tt39;
    tt90=1.0*r(1,0)*r(2,0);
    tt91=V*(1.0*r(1,1)*r(2,1)*lambda+(tt45+2.0*r(1,1)*r(2,1)+tt90)*mu);
    tt92=V*(1.0*r(0,2)*r(1,1)*lambda+1.0*r(0,1)*r(1,2)*mu);
    tt93=V*(1.0*r(1,1)*r(1,2)*lambda+1.0*r(1,1)*r(1,2)*mu);
    tt94=V*(1.0*r(1,1)*r(2,2)*lambda+1.0*r(1,2)*r(2,1)*mu);
    tt95=1.0*f(0,0)*r(1,0);
    tt96=1.0*f(1,0)*r(1,0);
    tt97=1.0*r(1,0)*f(2,0);
    tt98=1.0*tt59;
    tt99=V*(1.0*r(0,2)*r(2,1)*lambda+1.0*r(0,1)*r(2,2)*mu);
    tt100=V*(1.0*r(1,2)*r(2,1)*lambda+1.0*r(1,1)*r(2,2)*mu);
    tt101=V*(1.0*r(2,1)*r(2,2)*lambda+1.0*r(2,1)*r(2,2)*mu);
    tt102=1.0*f(0,0)*r(2,0);
    tt103=1.0*f(1,0)*r(2,0);
    tt104=1.0*f(2,0)*r(2,0);
    tt105=V*(1.0*r(0,2)*r(1,2)*lambda+(2.0*r(0,2)*r(1,2)+tt17+tt77)*mu);
    tt106=V*(1.0*r(0,2)*r(2,2)*lambda+(2.0*r(0,2)*r(2,2)+tt20+tt79)*mu);
    tt107=2.0*tt10;
    tt108=V*(1.0*r(1,2)*r(2,2)*lambda+(2.0*r(1,2)*r(2,2)+tt44+tt90)*mu);
    (*hess)(0,0)=V*(1.0*tt12*lambda+(tt16+tt14+2.0*tt12)*mu);
    (*hess)(0,1)=tt19;
    (*hess)(0,2)=tt22;
    (*hess)(0,3)=tt23;
    (*hess)(0,4)=tt24;
    (*hess)(0,5)=tt25;
    (*hess)(0,6)=tt26;
    (*hess)(0,7)=tt27;
    (*hess)(0,8)=tt28;
    (*hess)(0,9)=V*(tt32+1.0*f(0,0)*r(0,0)*lambda+(tt31+tt30+tt29+2.0*f(0,0)*r(0,0))*mu);
    (*hess)(0,10)=V*(1.0*r(0,0)*f(1,0)*lambda+(tt34+tt33+2.0*r(0,0)*f(1,0))*mu);
    (*hess)(0,11)=V*(1.0*r(0,0)*f(2,0)*lambda+(tt36+tt35+2.0*r(0,0)*f(2,0))*mu);
    (*hess)(0,12)=V*(1.0*r(0,0)*f(0,1)*lambda+(tt37+1.0*f(0,0)*r(0,1))*mu);
    (*hess)(0,13)=V*(1.0*r(0,0)*f(1,1)*lambda+1.0*r(0,1)*f(1,0)*mu);
    (*hess)(0,14)=V*(1.0*r(0,0)*f(2,1)*lambda+1.0*r(0,1)*f(2,0)*mu);
    (*hess)(0,15)=V*(1.0*r(0,0)*f(0,2)*lambda+(tt38+1.0*f(0,0)*r(0,2))*mu);
    (*hess)(0,16)=V*(1.0*r(0,0)*f(1,2)*lambda+1.0*r(0,2)*f(1,0)*mu);
    (*hess)(0,17)=V*(1.0*r(0,0)*f(2,2)*lambda+1.0*r(0,2)*f(2,0)*mu);
    (*hess)(1,0)=tt19;
    (*hess)(1,1)=V*(1.0*tt39*lambda+(tt43+tt41+2.0*tt39)*mu);
    (*hess)(1,2)=tt46;
    (*hess)(1,3)=tt47;
    (*hess)(1,4)=tt48;
    (*hess)(1,5)=tt49;
    (*hess)(1,6)=tt50;
    (*hess)(1,7)=tt51;
    (*hess)(1,8)=tt52;
    (*hess)(1,9)=V*(1.0*f(0,0)*r(1,0)*lambda+(tt54+tt53+2.0*f(0,0)*r(1,0))*mu);
    (*hess)(1,10)=V*(tt32+1.0*f(1,0)*r(1,0)*lambda+(tt31+tt56+tt55+2.0*f(1,0)*r(1,0))*mu);
    (*hess)(1,11)=V*(1.0*r(1,0)*f(2,0)*lambda+(tt58+tt57+2.0*r(1,0)*f(2,0))*mu);
    (*hess)(1,12)=V*(1.0*f(0,1)*r(1,0)*lambda+1.0*f(0,0)*r(1,1)*mu);
    (*hess)(1,13)=V*(1.0*r(1,0)*f(1,1)*lambda+(tt37+1.0*f(1,0)*r(1,1))*mu);
    (*hess)(1,14)=V*(1.0*r(1,0)*f(2,1)*lambda+1.0*r(1,1)*f(2,0)*mu);
    (*hess)(1,15)=V*(1.0*f(0,2)*r(1,0)*lambda+1.0*f(0,0)*r(1,2)*mu);
    (*hess)(1,16)=V*(1.0*r(1,0)*f(1,2)*lambda+(tt38+1.0*f(1,0)*r(1,2))*mu);
    (*hess)(1,17)=V*(1.0*r(1,0)*f(2,2)*lambda+1.0*r(1,2)*f(2,0)*mu);
    (*hess)(2,0)=tt22;
    (*hess)(2,1)=tt46;
    (*hess)(2,2)=V*(1.0*tt59*lambda+(tt63+tt61+2.0*tt59)*mu);
    (*hess)(2,3)=tt64;
    (*hess)(2,4)=tt65;
    (*hess)(2,5)=tt66;
    (*hess)(2,6)=tt67;
    (*hess)(2,7)=tt68;
    (*hess)(2,8)=tt69;
    (*hess)(2,9)=V*(1.0*f(0,0)*r(2,0)*lambda+(tt71+tt70+2.0*f(0,0)*r(2,0))*mu);
    (*hess)(2,10)=V*(1.0*f(1,0)*r(2,0)*lambda+(tt73+tt72+2.0*f(1,0)*r(2,0))*mu);
    (*hess)(2,11)=V*(tt32+1.0*f(2,0)*r(2,0)*lambda+(tt75+tt74+tt31+2.0*f(2,0)*r(2,0))*mu);
    (*hess)(2,12)=V*(1.0*f(0,1)*r(2,0)*lambda+1.0*f(0,0)*r(2,1)*mu);
    (*hess)(2,13)=V*(1.0*f(1,1)*r(2,0)*lambda+1.0*f(1,0)*r(2,1)*mu);
    (*hess)(2,14)=V*(1.0*r(2,0)*f(2,1)*lambda+(tt37+1.0*f(2,0)*r(2,1))*mu);
    (*hess)(2,15)=V*(1.0*f(0,2)*r(2,0)*lambda+1.0*f(0,0)*r(2,2)*mu);
    (*hess)(2,16)=V*(1.0*f(1,2)*r(2,0)*lambda+1.0*f(1,0)*r(2,2)*mu);
    (*hess)(2,17)=V*(1.0*r(2,0)*f(2,2)*lambda+(tt38+1.0*f(2,0)*r(2,2))*mu);
    (*hess)(3,0)=tt23;
    (*hess)(3,1)=tt47;
    (*hess)(3,2)=tt64;
    (*hess)(3,3)=V*(1.0*tt13*lambda+(tt16+2.0*tt13+tt76)*mu);
    (*hess)(3,4)=tt78;
    (*hess)(3,5)=tt80;
    (*hess)(3,6)=tt81;
    (*hess)(3,7)=tt82;
    (*hess)(3,8)=tt83;
    (*hess)(3,9)=V*(1.0*f(0,0)*r(0,1)*lambda+(tt37+1.0*r(0,0)*f(0,1))*mu);
    (*hess)(3,10)=V*(1.0*r(0,1)*f(1,0)*lambda+1.0*r(0,0)*f(1,1)*mu);
    (*hess)(3,11)=V*(1.0*r(0,1)*f(2,0)*lambda+1.0*r(0,0)*f(2,1)*mu);
    (*hess)(3,12)=V*(tt32+1.0*f(0,1)*r(0,1)*lambda+(tt85+tt30+2.0*f(0,1)*r(0,1)+tt84)*mu);
    (*hess)(3,13)=V*(1.0*r(0,1)*f(1,1)*lambda+(tt34+2.0*r(0,1)*f(1,1)+tt86)*mu);
    (*hess)(3,14)=V*(1.0*r(0,1)*f(2,1)*lambda+(tt36+2.0*r(0,1)*f(2,1)+tt87)*mu);
    (*hess)(3,15)=V*(1.0*r(0,1)*f(0,2)*lambda+(tt88+1.0*f(0,1)*r(0,2))*mu);
    (*hess)(3,16)=V*(1.0*r(0,1)*f(1,2)*lambda+1.0*r(0,2)*f(1,1)*mu);
    (*hess)(3,17)=V*(1.0*r(0,1)*f(2,2)*lambda+1.0*r(0,2)*f(2,1)*mu);
    (*hess)(4,0)=tt24;
    (*hess)(4,1)=tt48;
    (*hess)(4,2)=tt65;
    (*hess)(4,3)=tt78;
    (*hess)(4,4)=V*(1.0*tt40*lambda+(tt43+2.0*tt40+tt89)*mu);
    (*hess)(4,5)=tt91;
    (*hess)(4,6)=tt92;
    (*hess)(4,7)=tt93;
    (*hess)(4,8)=tt94;
    (*hess)(4,9)=V*(1.0*f(0,0)*r(1,1)*lambda+1.0*f(0,1)*r(1,0)*mu);
    (*hess)(4,10)=V*(1.0*f(1,0)*r(1,1)*lambda+(tt37+1.0*r(1,0)*f(1,1))*mu);
    (*hess)(4,11)=V*(1.0*r(1,1)*f(2,0)*lambda+1.0*r(1,0)*f(2,1)*mu);
    (*hess)(4,12)=V*(1.0*f(0,1)*r(1,1)*lambda+(tt54+2.0*f(0,1)*r(1,1)+tt95)*mu);
    (*hess)(4,13)=V*(tt32+1.0*f(1,1)*r(1,1)*lambda+(tt85+tt56+2.0*f(1,1)*r(1,1)+tt96)*mu);
    (*hess)(4,14)=V*(1.0*r(1,1)*f(2,1)*lambda+(tt58+2.0*r(1,1)*f(2,1)+tt97)*mu);
    (*hess)(4,15)=V*(1.0*f(0,2)*r(1,1)*lambda+1.0*f(0,1)*r(1,2)*mu);
    (*hess)(4,16)=V*(1.0*r(1,1)*f(1,2)*lambda+(tt88+1.0*f(1,1)*r(1,2))*mu);
    (*hess)(4,17)=V*(1.0*r(1,1)*f(2,2)*lambda+1.0*r(1,2)*f(2,1)*mu);
    (*hess)(5,0)=tt25;
    (*hess)(5,1)=tt49;
    (*hess)(5,2)=tt66;
    (*hess)(5,3)=tt80;
    (*hess)(5,4)=tt91;
    (*hess)(5,5)=V*(1.0*tt60*lambda+(tt63+2.0*tt60+tt98)*mu);
    (*hess)(5,6)=tt99;
    (*hess)(5,7)=tt100;
    (*hess)(5,8)=tt101;
    (*hess)(5,9)=V*(1.0*f(0,0)*r(2,1)*lambda+1.0*f(0,1)*r(2,0)*mu);
    (*hess)(5,10)=V*(1.0*f(1,0)*r(2,1)*lambda+1.0*f(1,1)*r(2,0)*mu);
    (*hess)(5,11)=V*(1.0*f(2,0)*r(2,1)*lambda+(tt37+1.0*r(2,0)*f(2,1))*mu);
    (*hess)(5,12)=V*(1.0*f(0,1)*r(2,1)*lambda+(tt71+2.0*f(0,1)*r(2,1)+tt102)*mu);
    (*hess)(5,13)=V*(1.0*f(1,1)*r(2,1)*lambda+(tt73+2.0*f(1,1)*r(2,1)+tt103)*mu);
    (*hess)(5,14)=V*(tt32+1.0*f(2,1)*r(2,1)*lambda+(tt75+tt85+2.0*f(2,1)*r(2,1)+tt104)*mu);
    (*hess)(5,15)=V*(1.0*f(0,2)*r(2,1)*lambda+1.0*f(0,1)*r(2,2)*mu);
    (*hess)(5,16)=V*(1.0*f(1,2)*r(2,1)*lambda+1.0*f(1,1)*r(2,2)*mu);
    (*hess)(5,17)=V*(1.0*r(2,1)*f(2,2)*lambda+(tt88+1.0*f(2,1)*r(2,2))*mu);
    (*hess)(6,0)=tt26;
    (*hess)(6,1)=tt50;
    (*hess)(6,2)=tt67;
    (*hess)(6,3)=tt81;
    (*hess)(6,4)=tt92;
    (*hess)(6,5)=tt99;
    (*hess)(6,6)=V*(1.0*tt15*lambda+(2.0*tt15+tt14+tt76)*mu);
    (*hess)(6,7)=tt105;
    (*hess)(6,8)=tt106;
    (*hess)(6,9)=V*(1.0*f(0,0)*r(0,2)*lambda+(tt38+1.0*r(0,0)*f(0,2))*mu);
    (*hess)(6,10)=V*(1.0*r(0,2)*f(1,0)*lambda+1.0*r(0,0)*f(1,2)*mu);
    (*hess)(6,11)=V*(1.0*r(0,2)*f(2,0)*lambda+1.0*r(0,0)*f(2,2)*mu);
    (*hess)(6,12)=V*(1.0*f(0,1)*r(0,2)*lambda+(tt88+1.0*r(0,1)*f(0,2))*mu);
    (*hess)(6,13)=V*(1.0*r(0,2)*f(1,1)*lambda+1.0*r(0,1)*f(1,2)*mu);
    (*hess)(6,14)=V*(1.0*r(0,2)*f(2,1)*lambda+1.0*r(0,1)*f(2,2)*mu);
    (*hess)(6,15)=V*(tt32+1.0*f(0,2)*r(0,2)*lambda+mu*(tt107+2.0*f(0,2)*r(0,2)+tt29+tt84));
    (*hess)(6,16)=V*(1.0*r(0,2)*f(1,2)*lambda+(2.0*r(0,2)*f(1,2)+tt33+tt86)*mu);
    (*hess)(6,17)=V*(1.0*r(0,2)*f(2,2)*lambda+(2.0*r(0,2)*f(2,2)+tt35+tt87)*mu);
    (*hess)(7,0)=tt27;
    (*hess)(7,1)=tt51;
    (*hess)(7,2)=tt68;
    (*hess)(7,3)=tt82;
    (*hess)(7,4)=tt93;
    (*hess)(7,5)=tt100;
    (*hess)(7,6)=tt105;
    (*hess)(7,7)=V*(1.0*tt42*lambda+(2.0*tt42+tt41+tt89)*mu);
    (*hess)(7,8)=tt108;
    (*hess)(7,9)=V*(1.0*f(0,0)*r(1,2)*lambda+1.0*f(0,2)*r(1,0)*mu);
    (*hess)(7,10)=V*(1.0*f(1,0)*r(1,2)*lambda+(tt38+1.0*r(1,0)*f(1,2))*mu);
    (*hess)(7,11)=V*(1.0*r(1,2)*f(2,0)*lambda+1.0*r(1,0)*f(2,2)*mu);
    (*hess)(7,12)=V*(1.0*f(0,1)*r(1,2)*lambda+1.0*f(0,2)*r(1,1)*mu);
    (*hess)(7,13)=V*(1.0*f(1,1)*r(1,2)*lambda+(tt88+1.0*r(1,1)*f(1,2))*mu);
    (*hess)(7,14)=V*(1.0*r(1,2)*f(2,1)*lambda+1.0*r(1,1)*f(2,2)*mu);
    (*hess)(7,15)=V*(1.0*f(0,2)*r(1,2)*lambda+(2.0*f(0,2)*r(1,2)+tt53+tt95)*mu);
    (*hess)(7,16)=V*(tt32+1.0*f(1,2)*r(1,2)*lambda+mu*(tt107+2.0*f(1,2)*r(1,2)+tt55+tt96));
    (*hess)(7,17)=V*(1.0*r(1,2)*f(2,2)*lambda+(2.0*r(1,2)*f(2,2)+tt57+tt97)*mu);
    (*hess)(8,0)=tt28;
    (*hess)(8,1)=tt52;
    (*hess)(8,2)=tt69;
    (*hess)(8,3)=tt83;
    (*hess)(8,4)=tt94;
    (*hess)(8,5)=tt101;
    (*hess)(8,6)=tt106;
    (*hess)(8,7)=tt108;
    (*hess)(8,8)=V*(1.0*tt62*lambda+(2.0*tt62+tt61+tt98)*mu);
    (*hess)(8,9)=V*(1.0*f(0,0)*r(2,2)*lambda+1.0*f(0,2)*r(2,0)*mu);
    (*hess)(8,10)=V*(1.0*f(1,0)*r(2,2)*lambda+1.0*f(1,2)*r(2,0)*mu);
    (*hess)(8,11)=V*(1.0*f(2,0)*r(2,2)*lambda+(tt38+1.0*r(2,0)*f(2,2))*mu);
    (*hess)(8,12)=V*(1.0*f(0,1)*r(2,2)*lambda+1.0*f(0,2)*r(2,1)*mu);
    (*hess)(8,13)=V*(1.0*f(1,1)*r(2,2)*lambda+1.0*f(1,2)*r(2,1)*mu);
    (*hess)(8,14)=V*(1.0*f(2,1)*r(2,2)*lambda+(tt88+1.0*r(2,1)*f(2,2))*mu);
    (*hess)(8,15)=V*(1.0*f(0,2)*r(2,2)*lambda+(2.0*f(0,2)*r(2,2)+tt70+tt102)*mu);
    (*hess)(8,16)=V*(1.0*f(1,2)*r(2,2)*lambda+(2.0*f(1,2)*r(2,2)+tt72+tt103)*mu);
    (*hess)(8,17)=V*(tt32+1.0*f(2,2)*r(2,2)*lambda+mu*(tt107+2.0*f(2,2)*r(2,2)+tt74+tt104));
    return E;
}
scalarD MaterialEnergy::evalFNonHK(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,scalarD minVol,
                                   Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    scalarD E;
    //input
    //scalarD V;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD mu;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;

    tt1=dimMask+f(2,2);
    tt2=f(1,0)*f(2,1)-f(1,1)*f(2,0);
    tt3=f(1,1)*tt1-f(1,2)*f(2,1);
    tt4=f(0,0)*tt3-f(0,1)*(f(1,0)*tt1-f(1,2)*f(2,0))+f(0,2)*tt2;
    tt5=log(std::max(tt4,minVol));
    tt6=1/tt4;
    tt7=f(0,2)*f(2,1)-f(0,1)*tt1;
    tt8=f(0,1)*f(1,2)-f(0,2)*f(1,1);
    tt9=f(1,2)*f(2,0)-f(1,0)*tt1;
    tt10=f(0,0)*tt1-f(0,2)*f(2,0);
    tt11=f(0,2)*f(1,0)-f(0,0)*f(1,2);
    tt12=f(0,1)*f(2,0)-f(0,0)*f(2,1);
    tt13=f(0,0)*f(1,1)-f(0,1)*f(1,0);
    E=V*(0.5*pow(tt5,2)*lambda+0.5*mu*(-2*tt5+pow(tt1,2)+pow(f(2,1),2)+pow(f(2,0),2)+pow(f(1,2),2)+pow(f(1,1),2)+pow(f(1,0),2)+pow(f(0,2),2)+pow(f(0,1),2)+pow(f(0,0),2)-3));
    grad[0]=V*(1.0*tt3*tt6*tt5*lambda+0.5*mu*(2*f(0,0)-2*tt3*tt6));
    grad[1]=V*(1.0*tt7*tt6*tt5*lambda+0.5*mu*(2*f(1,0)-2*tt7*tt6));
    grad[2]=V*(1.0*tt8*tt6*tt5*lambda+0.5*mu*(2*f(2,0)-2*tt8*tt6));
    grad[3]=V*(1.0*tt9*tt6*tt5*lambda+0.5*mu*(2*f(0,1)-2*tt9*tt6));
    grad[4]=V*(1.0*tt10*tt6*tt5*lambda+0.5*mu*(2*f(1,1)-2*tt10*tt6));
    grad[5]=V*(1.0*tt11*tt6*tt5*lambda+0.5*mu*(2*f(2,1)-2*tt11*tt6));
    grad[6]=V*(1.0*tt2*tt6*tt5*lambda+0.5*mu*(2*f(0,2)-2*tt2*tt6));
    grad[7]=V*(1.0*tt12*tt6*tt5*lambda+0.5*mu*(2*f(1,2)-2*tt12*tt6));
    grad[8]=V*(1.0*tt13*tt6*tt5*lambda+0.5*mu*(2*tt1-2*tt13*tt6));
    if(!hess)return E;

    tt14=pow(tt3,2);
    tt15=1/pow(tt4,2);
    tt16=V*(-1.0*tt7*tt3*tt15*tt5*lambda+1.0*tt7*tt3*tt15*lambda+1.0*mu*tt7*tt3*tt15);
    tt17=V*(-1.0*tt8*tt3*tt15*tt5*lambda+1.0*tt8*tt3*tt15*lambda+1.0*tt8*mu*tt3*tt15);
    tt18=V*(-1.0*tt9*tt3*tt15*tt5*lambda+1.0*tt9*tt3*tt15*lambda+1.0*mu*tt9*tt3*tt15);
    tt19=V*(1.0*tt1*tt6*tt5*lambda-1.0*tt10*tt3*tt15*tt5*lambda+1.0*tt10*tt3*tt15*lambda+0.5*mu*(2*tt10*tt3*tt15-2*tt1*tt6));
    tt20=V*(-1.0*f(1,2)*tt6*tt5*lambda-1.0*tt11*tt3*tt15*tt5*lambda+1.0*tt11*tt3*tt15*lambda+0.5*mu*(2*f(1,2)*tt6+2*tt11*tt3*tt15));
    tt21=V*(-1.0*tt2*tt3*tt15*tt5*lambda+1.0*tt2*tt3*tt15*lambda+1.0*tt2*mu*tt3*tt15);
    tt22=V*(-1.0*f(2,1)*tt6*tt5*lambda-1.0*tt12*tt3*tt15*tt5*lambda+1.0*tt12*tt3*tt15*lambda+0.5*mu*(2*f(2,1)*tt6+2*tt12*tt3*tt15));
    tt23=V*(1.0*f(1,1)*tt6*tt5*lambda-1.0*tt13*tt3*tt15*tt5*lambda+1.0*tt13*tt3*tt15*lambda+0.5*mu*(2*tt13*tt3*tt15-2*f(1,1)*tt6));
    tt24=pow(tt7,2);
    tt25=V*(-1.0*tt8*tt7*tt15*tt5*lambda+1.0*tt8*tt7*tt15*lambda+1.0*tt8*mu*tt7*tt15);
    tt26=-dimMask-f(2,2);
    tt27=V*(1.0*tt26*tt6*tt5*lambda-1.0*tt7*tt9*tt15*tt5*lambda+1.0*tt7*tt9*tt15*lambda+0.5*mu*(2*tt7*tt9*tt15-2*tt26*tt6));
    tt28=V*(-1.0*tt10*tt7*tt15*tt5*lambda+1.0*tt10*tt7*tt15*lambda+1.0*mu*tt10*tt7*tt15);
    tt29=V*(1.0*f(0,2)*tt6*tt5*lambda-1.0*tt11*tt7*tt15*tt5*lambda+1.0*tt11*tt7*tt15*lambda+0.5*mu*(2*tt11*tt7*tt15-2*f(0,2)*tt6));
    tt30=V*(1.0*f(2,1)*tt6*tt5*lambda-1.0*tt2*tt7*tt15*tt5*lambda+1.0*tt2*tt7*tt15*lambda+0.5*mu*(2*tt2*tt7*tt15-2*f(2,1)*tt6));
    tt31=V*(-1.0*tt12*tt7*tt15*tt5*lambda+1.0*tt12*tt7*tt15*lambda+1.0*tt12*mu*tt7*tt15);
    tt32=V*(-1.0*f(0,1)*tt6*tt5*lambda-1.0*tt13*tt7*tt15*tt5*lambda+1.0*tt13*tt7*tt15*lambda+0.5*mu*(2*f(0,1)*tt6+2*tt13*tt7*tt15));
    tt33=pow(tt8,2);
    tt34=V*(1.0*f(1,2)*tt6*tt5*lambda-1.0*tt8*tt9*tt15*tt5*lambda+1.0*tt8*tt9*tt15*lambda+0.5*mu*(2*tt8*tt9*tt15-2*f(1,2)*tt6));
    tt35=V*(-1.0*f(0,2)*tt6*tt5*lambda-1.0*tt8*tt10*tt15*tt5*lambda+1.0*tt8*tt10*tt15*lambda+0.5*mu*(2*f(0,2)*tt6+2*tt8*tt10*tt15));
    tt36=V*(-1.0*tt11*tt8*tt15*tt5*lambda+1.0*tt11*tt8*tt15*lambda+1.0*tt11*tt8*mu*tt15);
    tt37=V*(-1.0*f(1,1)*tt6*tt5*lambda-1.0*tt8*tt2*tt15*tt5*lambda+1.0*tt8*tt2*tt15*lambda+0.5*mu*(2*f(1,1)*tt6+2*tt8*tt2*tt15));
    tt38=V*(1.0*f(0,1)*tt6*tt5*lambda-1.0*tt8*tt12*tt15*tt5*lambda+1.0*tt8*tt12*tt15*lambda+0.5*mu*(2*tt8*tt12*tt15-2*f(0,1)*tt6));
    tt39=V*(-1.0*tt13*tt8*tt15*tt5*lambda+1.0*tt13*tt8*tt15*lambda+1.0*tt13*tt8*mu*tt15);
    tt40=pow(tt9,2);
    tt41=V*(-1.0*tt10*tt9*tt15*tt5*lambda+1.0*tt10*tt9*tt15*lambda+1.0*mu*tt10*tt9*tt15);
    tt42=V*(-1.0*tt11*tt9*tt15*tt5*lambda+1.0*tt11*tt9*tt15*lambda+1.0*tt11*mu*tt9*tt15);
    tt43=V*(-1.0*tt2*tt9*tt15*tt5*lambda+1.0*tt2*tt9*tt15*lambda+1.0*tt2*mu*tt9*tt15);
    tt44=V*(1.0*f(2,0)*tt6*tt5*lambda-1.0*tt12*tt9*tt15*tt5*lambda+1.0*tt12*tt9*tt15*lambda+0.5*mu*(2*tt12*tt9*tt15-2*f(2,0)*tt6));
    tt45=V*(-1.0*f(1,0)*tt6*tt5*lambda-1.0*tt13*tt9*tt15*tt5*lambda+1.0*tt13*tt9*tt15*lambda+0.5*mu*(2*f(1,0)*tt6+2*tt13*tt9*tt15));
    tt46=pow(tt10,2);
    tt47=V*(-1.0*tt11*tt10*tt15*tt5*lambda+1.0*tt11*tt10*tt15*lambda+1.0*tt11*mu*tt10*tt15);
    tt48=V*(-1.0*f(2,0)*tt6*tt5*lambda-1.0*tt2*tt10*tt15*tt5*lambda+1.0*tt2*tt10*tt15*lambda+0.5*mu*(2*f(2,0)*tt6+2*tt2*tt10*tt15));
    tt49=V*(-1.0*tt12*tt10*tt15*tt5*lambda+1.0*tt12*tt10*tt15*lambda+1.0*tt12*mu*tt10*tt15);
    tt50=V*(1.0*f(0,0)*tt6*tt5*lambda-1.0*tt13*tt10*tt15*tt5*lambda+1.0*tt13*tt10*tt15*lambda+0.5*mu*(2*tt13*tt10*tt15-2*f(0,0)*tt6));
    tt51=pow(tt11,2);
    tt52=V*(1.0*f(1,0)*tt6*tt5*lambda-1.0*tt11*tt2*tt15*tt5*lambda+1.0*tt11*tt2*tt15*lambda+0.5*mu*(2*tt11*tt2*tt15-2*f(1,0)*tt6));
    tt53=V*(-1.0*f(0,0)*tt6*tt5*lambda-1.0*tt11*tt12*tt15*tt5*lambda+1.0*tt11*tt12*tt15*lambda+0.5*mu*(2*f(0,0)*tt6+2*tt11*tt12*tt15));
    tt54=V*(-1.0*tt13*tt11*tt15*tt5*lambda+1.0*tt13*tt11*tt15*lambda+1.0*tt13*tt11*mu*tt15);
    tt55=pow(tt2,2);
    tt56=V*(-1.0*tt12*tt2*tt15*tt5*lambda+1.0*tt12*tt2*tt15*lambda+1.0*tt12*tt2*mu*tt15);
    tt57=V*(-1.0*tt13*tt2*tt15*tt5*lambda+1.0*tt13*tt2*tt15*lambda+1.0*tt13*tt2*mu*tt15);
    tt58=pow(tt12,2);
    tt59=V*(-1.0*tt13*tt12*tt15*tt5*lambda+1.0*tt13*tt12*tt15*lambda+1.0*tt13*tt12*mu*tt15);
    tt60=pow(tt13,2);
    (*hess)(0,0)=V*(-1.0*tt14*tt15*tt5*lambda+1.0*tt14*tt15*lambda+0.5*mu*(2*tt14*tt15+2));
    (*hess)(0,1)=tt16;
    (*hess)(0,2)=tt17;
    (*hess)(0,3)=tt18;
    (*hess)(0,4)=tt19;
    (*hess)(0,5)=tt20;
    (*hess)(0,6)=tt21;
    (*hess)(0,7)=tt22;
    (*hess)(0,8)=tt23;
    (*hess)(1,0)=tt16;
    (*hess)(1,1)=V*(-1.0*tt24*tt15*tt5*lambda+1.0*tt24*tt15*lambda+0.5*mu*(2*tt24*tt15+2));
    (*hess)(1,2)=tt25;
    (*hess)(1,3)=tt27;
    (*hess)(1,4)=tt28;
    (*hess)(1,5)=tt29;
    (*hess)(1,6)=tt30;
    (*hess)(1,7)=tt31;
    (*hess)(1,8)=tt32;
    (*hess)(2,0)=tt17;
    (*hess)(2,1)=tt25;
    (*hess)(2,2)=V*(-1.0*tt33*tt15*tt5*lambda+1.0*tt33*tt15*lambda+0.5*mu*(2*tt33*tt15+2));
    (*hess)(2,3)=tt34;
    (*hess)(2,4)=tt35;
    (*hess)(2,5)=tt36;
    (*hess)(2,6)=tt37;
    (*hess)(2,7)=tt38;
    (*hess)(2,8)=tt39;
    (*hess)(3,0)=tt18;
    (*hess)(3,1)=tt27;
    (*hess)(3,2)=tt34;
    (*hess)(3,3)=V*(-1.0*tt40*tt15*tt5*lambda+1.0*tt40*tt15*lambda+0.5*mu*(2*tt40*tt15+2));
    (*hess)(3,4)=tt41;
    (*hess)(3,5)=tt42;
    (*hess)(3,6)=tt43;
    (*hess)(3,7)=tt44;
    (*hess)(3,8)=tt45;
    (*hess)(4,0)=tt19;
    (*hess)(4,1)=tt28;
    (*hess)(4,2)=tt35;
    (*hess)(4,3)=tt41;
    (*hess)(4,4)=V*(-1.0*tt46*tt15*tt5*lambda+1.0*tt46*tt15*lambda+0.5*mu*(2*tt46*tt15+2));
    (*hess)(4,5)=tt47;
    (*hess)(4,6)=tt48;
    (*hess)(4,7)=tt49;
    (*hess)(4,8)=tt50;
    (*hess)(5,0)=tt20;
    (*hess)(5,1)=tt29;
    (*hess)(5,2)=tt36;
    (*hess)(5,3)=tt42;
    (*hess)(5,4)=tt47;
    (*hess)(5,5)=V*(-1.0*tt51*tt15*tt5*lambda+1.0*tt51*tt15*lambda+0.5*mu*(2*tt51*tt15+2));
    (*hess)(5,6)=tt52;
    (*hess)(5,7)=tt53;
    (*hess)(5,8)=tt54;
    (*hess)(6,0)=tt21;
    (*hess)(6,1)=tt30;
    (*hess)(6,2)=tt37;
    (*hess)(6,3)=tt43;
    (*hess)(6,4)=tt48;
    (*hess)(6,5)=tt52;
    (*hess)(6,6)=V*(-1.0*tt55*tt15*tt5*lambda+1.0*tt55*tt15*lambda+0.5*mu*(2*tt55*tt15+2));
    (*hess)(6,7)=tt56;
    (*hess)(6,8)=tt57;
    (*hess)(7,0)=tt22;
    (*hess)(7,1)=tt31;
    (*hess)(7,2)=tt38;
    (*hess)(7,3)=tt44;
    (*hess)(7,4)=tt49;
    (*hess)(7,5)=tt53;
    (*hess)(7,6)=tt56;
    (*hess)(7,7)=V*(-1.0*tt58*tt15*tt5*lambda+1.0*tt58*tt15*lambda+0.5*mu*(2*tt58*tt15+2));
    (*hess)(7,8)=tt59;
    (*hess)(8,0)=tt23;
    (*hess)(8,1)=tt32;
    (*hess)(8,2)=tt39;
    (*hess)(8,3)=tt45;
    (*hess)(8,4)=tt50;
    (*hess)(8,5)=tt54;
    (*hess)(8,6)=tt57;
    (*hess)(8,7)=tt59;
    (*hess)(8,8)=V*(-1.0*tt60*tt15*tt5*lambda+1.0*tt60*tt15*lambda+0.5*mu*(2*tt60*tt15+2));
    return E;
}
scalarD MaterialEnergy::evalFFung(const Mat3& f,scalarD lambda,scalarD mu,scalarD lambda2,scalarD mu2,scalarD c,scalarD V,scalarD dimMask,
                                  Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess)
{
    //input
    //scalarD V;
    //scalarD c;
    //scalarD dimMask;
    //Mat3d f;
    //scalarD lambda;
    //scalarD lambda2;
    //scalarD mu;
    //scalarD mu2;
    scalarD E;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;

    tt1=1.0*f(0,0);
    tt2=tt1-1;
    tt3=f(1,0)+f(0,1);
    tt4=pow(tt3,2);
    tt5=1.0*f(1,1);
    tt6=tt5-1;
    tt7=f(2,0)+f(0,2);
    tt8=pow(tt7,2);
    tt9=f(2,1)+f(1,2);
    tt10=pow(tt9,2);
    tt11=1.0*f(2,2);
    tt12=-dimMask;
    tt13=tt12+tt11;
    tt14=pow(tt13,2)+0.5*tt10+0.5*tt8+pow(tt6,2)+0.5*tt4+pow(tt2,2);
    tt15=tt12+tt11+tt5+tt1-2;
    tt16=pow(tt15,2);
    tt17=exp(lambda2*tt16+2.0*mu2*tt14);
    tt18=2.0*lambda2*tt15;
    tt19=tt18+4.0*tt2*mu2;
    tt20=2.0*tt15*lambda;
    tt21=0.5*(2.0*tt3*c*mu2*tt17+2.0*tt3*mu)*V;
    tt22=0.5*(2.0*tt7*c*mu2*tt17+2.0*tt7*mu)*V;
    tt23=tt18+4.0*tt6*mu2;
    tt24=0.5*(2.0*tt9*c*mu2*tt17+2.0*tt9*mu)*V;
    tt25=tt18+4.0*mu2*tt13;
    E=0.5*V*(tt16*lambda+c*(tt17-1)+2.0*mu*tt14);
    grad[0]=0.5*V*(tt20+c*tt19*tt17+4.0*tt2*mu);
    grad[1]=tt21;
    grad[2]=tt22;
    grad[3]=tt21;
    grad[4]=0.5*V*(tt20+c*tt23*tt17+4.0*tt6*mu);
    grad[5]=tt24;
    grad[6]=tt22;
    grad[7]=tt24;
    grad[8]=0.5*V*(tt20+c*tt25*tt17+4.0*mu*tt13);
    if(!hess)
        return E;

    tt26=4.0*mu;
    tt27=c*(4.0*mu2+2.0*lambda2)*tt17;
    tt28=2.0*lambda;
    tt29=1.0*tt3*c*mu2*tt19*tt17*V;
    tt30=1.0*tt7*c*mu2*tt19*tt17*V;
    tt31=2.0*c*lambda2*tt17;
    tt32=0.5*V*(tt28+c*tt19*tt23*tt17+tt31);
    tt33=1.0*tt9*c*mu2*tt19*tt17*V;
    tt34=0.5*V*(tt28+c*tt19*tt25*tt17+tt31);
    tt35=2.0*mu;
    tt36=2.0*c*mu2*tt17;
    tt37=pow(mu2,2);
    tt38=0.5*(4.0*tt4*c*tt37*tt17+tt36+tt35)*V;
    tt39=2.0*tt3*tt7*c*tt37*tt17*V;
    tt40=1.0*tt3*c*mu2*tt23*tt17*V;
    tt41=2.0*tt3*tt9*c*tt37*tt17*V;
    tt42=1.0*tt3*c*mu2*tt25*tt17*V;
    tt43=0.5*(4.0*tt8*c*tt37*tt17+tt36+tt35)*V;
    tt44=1.0*tt7*c*mu2*tt23*tt17*V;
    tt45=2.0*tt7*tt9*c*tt37*tt17*V;
    tt46=1.0*tt7*c*mu2*tt25*tt17*V;
    tt47=1.0*tt9*c*mu2*tt23*tt17*V;
    tt48=0.5*V*(tt28+c*tt23*tt25*tt17+tt31);
    tt49=0.5*(4.0*tt10*c*tt37*tt17+tt36+tt35)*V;
    tt50=1.0*tt9*c*mu2*tt25*tt17*V;
    (*hess)(0,0)=0.5*V*(tt28+c*pow(tt19,2)*tt17+tt27+tt26);
    (*hess)(0,1)=tt29;
    (*hess)(0,2)=tt30;
    (*hess)(0,3)=tt29;
    (*hess)(0,4)=tt32;
    (*hess)(0,5)=tt33;
    (*hess)(0,6)=tt30;
    (*hess)(0,7)=tt33;
    (*hess)(0,8)=tt34;
    (*hess)(1,0)=tt29;
    (*hess)(1,1)=tt38;
    (*hess)(1,2)=tt39;
    (*hess)(1,3)=tt38;
    (*hess)(1,4)=tt40;
    (*hess)(1,5)=tt41;
    (*hess)(1,6)=tt39;
    (*hess)(1,7)=tt41;
    (*hess)(1,8)=tt42;
    (*hess)(2,0)=tt30;
    (*hess)(2,1)=tt39;
    (*hess)(2,2)=tt43;
    (*hess)(2,3)=tt39;
    (*hess)(2,4)=tt44;
    (*hess)(2,5)=tt45;
    (*hess)(2,6)=tt43;
    (*hess)(2,7)=tt45;
    (*hess)(2,8)=tt46;
    (*hess)(3,0)=tt29;
    (*hess)(3,1)=tt38;
    (*hess)(3,2)=tt39;
    (*hess)(3,3)=tt38;
    (*hess)(3,4)=tt40;
    (*hess)(3,5)=tt41;
    (*hess)(3,6)=tt39;
    (*hess)(3,7)=tt41;
    (*hess)(3,8)=tt42;
    (*hess)(4,0)=tt32;
    (*hess)(4,1)=tt40;
    (*hess)(4,2)=tt44;
    (*hess)(4,3)=tt40;
    (*hess)(4,4)=0.5*V*(tt28+c*pow(tt23,2)*tt17+tt27+tt26);
    (*hess)(4,5)=tt47;
    (*hess)(4,6)=tt44;
    (*hess)(4,7)=tt47;
    (*hess)(4,8)=tt48;
    (*hess)(5,0)=tt33;
    (*hess)(5,1)=tt41;
    (*hess)(5,2)=tt45;
    (*hess)(5,3)=tt41;
    (*hess)(5,4)=tt47;
    (*hess)(5,5)=tt49;
    (*hess)(5,6)=tt45;
    (*hess)(5,7)=tt49;
    (*hess)(5,8)=tt50;
    (*hess)(6,0)=tt30;
    (*hess)(6,1)=tt39;
    (*hess)(6,2)=tt43;
    (*hess)(6,3)=tt39;
    (*hess)(6,4)=tt44;
    (*hess)(6,5)=tt45;
    (*hess)(6,6)=tt43;
    (*hess)(6,7)=tt45;
    (*hess)(6,8)=tt46;
    (*hess)(7,0)=tt33;
    (*hess)(7,1)=tt41;
    (*hess)(7,2)=tt45;
    (*hess)(7,3)=tt41;
    (*hess)(7,4)=tt47;
    (*hess)(7,5)=tt49;
    (*hess)(7,6)=tt45;
    (*hess)(7,7)=tt49;
    (*hess)(7,8)=tt50;
    (*hess)(8,0)=tt34;
    (*hess)(8,1)=tt42;
    (*hess)(8,2)=tt46;
    (*hess)(8,3)=tt42;
    (*hess)(8,4)=tt48;
    (*hess)(8,5)=tt50;
    (*hess)(8,6)=tt46;
    (*hess)(8,7)=tt50;
    (*hess)(8,8)=0.5*V*(tt28+c*pow(tt25,2)*tt17+tt27+tt26);
    return E;
}

EnergyPool::EnergyPool(const FEMBody& body):_body(body) {}
EnergyPool::ConstIterator EnergyPool::begin(const type_info& t) const
{
    const structE& elist=_energys.find(std::string(t.name()))->second;
    return elist._elist.begin()+elist._lock;
}
EnergyPool::ConstIterator EnergyPool::end(const type_info& t) const
{
    const structE& elist=_energys.find(std::string(t.name()))->second;
    return elist._elist.end();
}
void EnergyPool::addEnergy(boost::shared_ptr<FEMEnergy> e)
{
    if(!e)return;
    _energys[typeid(*e).name()]._elist.push_back(e);
}
void EnergyPool::delEnergy(boost::shared_ptr<FEMEnergy> e)
{
    if(!e)return;
    structE& elist=_energys[typeid(*e).name()];
    ELIST& ess=elist._elist;
    ess.erase(std::find(ess.begin()+elist._lock,ess.end(),e));
}
void EnergyPool::delEnergy(const std::string& name,sizeType id)
{
    structE& elist=_energys[name];
    ELIST& ess=elist._elist;
    vector<bool> flag(_body.getV(_body.nrV()-1)._index+1,false);
    for(ELIST::iterator beg=ess.begin()+elist._lock,end=ess.end(); beg!=end; beg++) {
        (*beg)->flagVertex(flag);
        if(flag[id]) {
            ess.erase(beg);
            return;
        }
    }
}
void EnergyPool::setMaterialType(MaterialEnergy::TYPE type)
{
    structE& elist=_energys[typeid(MaterialEnergy).name()];
    ELIST& ess=elist._elist;
    for(ELIST::iterator beg=ess.begin()+elist._lock,end=ess.end(); beg!=end; beg++)
        boost::dynamic_pointer_cast<MaterialEnergy>(*beg)->setType(type);
}
MaterialEnergy::TYPE EnergyPool::getMaterialType() const
{
    EMAP::const_iterator iter=_energys.find(typeid(MaterialEnergy).name());
    if(iter == _energys.end())
        return MaterialEnergy::LINEAR;

    sizeType type=-1;
    const structE& elist=iter->second;
    const ELIST& ess=elist._elist;
    for(ELIST::const_iterator beg=ess.begin(),end=ess.end(); beg!=end; beg++) {
        MaterialEnergy::TYPE t=boost::dynamic_pointer_cast<MaterialEnergy>(*beg)->getType();
        ASSERT_MSG(type == -1 || t == type,"User using a multimaterial body!")
        type=t;
    }
    return (MaterialEnergy::TYPE)type;
}
void EnergyPool::lock(bool setLock,const TypeChooser& c)
{
    for(EMAP::iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        structE& elist=beg->second;
        const ELIST& ess=elist._elist;
        if(ess.empty() || !c(typeid(*(ess[0]))))continue;
        elist._lock=setLock?(sizeType)elist._elist.size():0;
    }
}
void EnergyPool::clear()
{
    for(EMAP::iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        structE& elist=beg->second;
        elist._elist.resize(elist._lock);
    }
}
void EnergyPool::clear(const std::string& name)
{
    if(_energys.find(name) != _energys.end()) {
        structE& elist=_energys[name];
        elist._elist.resize(elist._lock);
    }
}
void EnergyPool::copy(const EnergyPool& other)
{
    _energys.clear();
    for(EMAP::const_iterator beg=other._energys.begin(),end=other._energys.end(); beg!=end; beg++) {
        const structE& elistO=beg->second;
        structE& elist=_energys[beg->first];

        elist._lock=elistO._lock;
        elist._elist=elistO._elist;
        for(ELIST::iterator ebeg=elist._elist.begin(),eend=elist._elist.end(); ebeg!=eend; ebeg++) {
            boost::shared_ptr<FEMEnergy> tmp=*ebeg;
            *ebeg=tmp->copy(_body);
        }
    }
}
scalarD EnergyPool::energyEval(Vec* f,TRIPS* H,scalar CF,scalar CH,const TypeChooser& c,bool lock) const
{
    scalarD ret=0.0f;
#define PARALLEL
#ifdef PARALLEL
    const OmpSettings& setting=OmpSettings::getOmpSettings();
    boost::shared_ptr<vector<Vec> > fParallel;
    if(f)fParallel.reset(new vector<Vec>(setting.nrThreads(),Vec::Zero(f->size())));
    boost::shared_ptr<vector<TRIPS> > HParallel;
    if(H)HParallel.reset(new vector<TRIPS>(setting.nrThreads(),TRIPS()));
    for(EMAP::const_iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        const structE& elist=beg->second;
        const ELIST& ess=elist._elist;
        if(ess.empty() || !c(typeid(*(ess[0]))))continue;

        sizeType nrE=(sizeType)ess.size();
#ifdef _MSC_VER
        OMP_PARALLEL_FOR_I(OMP_ADD(ret))
#endif
        for(sizeType i=lock?elist._lock:0; i<nrE; i++)
            ret+=ess[i]->eval(fParallel ? &((*fParallel)[omp_get_thread_num()]) : NULL,
                              HParallel ? &((*HParallel)[omp_get_thread_num()]) : NULL,CF,CH);
    }
    for(sizeType i=0; i<setting.nrThreads(); i++) {
        if(f)(*f)+=(*fParallel)[i];
        if(H)H->insert(H->end(),(*HParallel)[i].begin(),(*HParallel)[i].end());
    }
#else
    for(EMAP::const_iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        const structE& elist=beg->second;
        const ELIST& ess=elist._elist;
        if(ess.empty() || !c(typeid(*(ess[0]))))continue;
        sizeType nrE=(sizeType)ess.size();
        for(sizeType i=lock?elist._lock:0; i<nrE; i++)
            ret+=ess[i]->eval(f,H,CF,CH);
    }
#endif
    return ret;
}
scalarD EnergyPool::energyEval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH,const TypeChooser& c,bool lock) const
{
    scalarD ret=0.0f;
    for(EMAP::const_iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        const structE& elist=beg->second;
        const ELIST& ess=elist._elist;
        if(ess.empty() || !c(typeid(*(ess[0]))))continue;
        sizeType nrE=(sizeType)ess.size();
        for(sizeType i=lock?elist._lock:0; i<nrE; i++)
            ret+=ess[i]->eval(U,f,H,CF,CH);
    }
    return ret;
}
void EnergyPool::flagVertex(vector<bool>& flag,const TypeChooser& c) const
{
    for(EMAP::const_iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        const structE& elist=beg->second;
        const ELIST& ess=elist._elist;
        if(!ess.empty() && c(typeid(*(ess[0])))) {
            sizeType nrE=(sizeType)ess.size();
            OMP_PARALLEL_FOR_
            for(sizeType i=beg->second._lock; i<nrE; i++)
                ess[i]->flagVertex(flag);
        }
    }
}
void EnergyPool::debugEnergy()
{
    Vec dpos,dposTmp;
    const_cast<FEMBody&>(_body).getDPos(dpos);
    dposTmp=dpos;
    dposTmp.setRandom();
    const_cast<FEMBody&>(_body).setDPos(dposTmp);
    for(EMAP::const_iterator beg=_energys.begin(),end=_energys.end(); beg!=end; beg++) {
        structE& elist=_energys[beg->first];
        ELIST& ess=elist._elist;
        for(ELIST::iterator beg=ess.begin()+elist._lock,end=ess.end(); beg!=end; beg++)
            (*beg)->debugEnergy(_body);
    }
    const_cast<FEMBody&>(_body).setDPos(dpos);
}
