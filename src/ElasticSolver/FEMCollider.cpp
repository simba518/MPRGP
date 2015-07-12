#include "FEMCollider.h"
#include "FEMSystem.h"
#include "FEMMesh.h"
#include "FEMUtils.h"

USE_PRJ_NAMESPACE

DebugFEMCollider::DebugFEMCollider(const std::string& path,sizeType dim):_os("Collider",path,false),_dim(dim) {}
DebugFEMCollider::~DebugFEMCollider()
{
    _os.appendPoints(_vss.begin(),_vss.end());
    _os.appendCells(_issv.begin(),_issv.end(),_dim == 2 ? VTKWriter<scalar>::TRIANGLE : VTKWriter<scalar>::TETRA);
    _os.appendCells(_issp.begin(),_issp.end(),VTKWriter<scalar>::POINT);
    _os.appendCells(_issl.begin(),_issl.end(),VTKWriter<scalar>::LINE);
    _os.appendCustomPointData("Color",_css.begin(),_css.end());
}
void DebugFEMCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV)
{
    _issl.push_back(Vec2i(_vss.size(),_vss.size()+1));
    _vss.push_back(v[0]->_pos);
    Vec3 pos1=Vec3::Zero();
    for(sizeType i=1; i<nrV; i++)
        pos1-=v[i]->_pos*(scalar)(coef[i].dot(coef[0]));
    _vss.push_back(pos1);
    _css.push_back(1.0f);
    _css.push_back(1.0f);
}
void DebugFEMCollider::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary)
{
    //INFOV("On CV: %ld %ld",c->_index,v->_index)
    if(_dim == 2) {
        _issv.push_back(Vec4i(_vss.size()+0,_vss.size()+1,_vss.size()+2,-1));
        _issp.push_back(Vec4i::Constant(_vss.size()+3));
    } else {
        _issv.push_back(Vec4i(_vss.size()+0,_vss.size()+1,_vss.size()+2,_vss.size()+3));
        _issp.push_back(Vec4i::Constant(_vss.size()+4));
    }
    Vec3 pos=Vec3::Zero();
    for(sizeType vid=0; vid<_dim+1; vid++) {
        _vss.push_back((*c)[vid]);
        _css.push_back(1.0f);
        pos+=(*c)[vid]*bary[vid];
    }
    _vss.push_back(pos);
    _css.push_back(0.0f);
}
void DebugFEMCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n)
{
    //INFOV("On VN: %ld, (%f,%f,%f)",v->_index,n[0],n[1],n[2])
    _issl.push_back(Vec2i(_vss.size(),_vss.size()+1));
    _vss.push_back(v->_pos);
    _vss.push_back(v->_pos+n);
    _css.push_back(1.0f);
    _css.push_back(1.0f);
}

class SelfCollEnergy : public FEMEnergy
{
public:
    //dist=b4+(b1-b4,b2-b4,b3-b4)(v1-v4,v2-v4,v3-v4)^{-1}(vp-v4)
    //Energy: 0.5f*K*dist*dist
    SelfCollEnergy(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV,scalar K):_K(K) {
        _nrV=nrV;
        for(sizeType i=0; i<nrV; i++) {
            _b[i]=b[i];
            _v[i]=v[i];
            _coef[i]=coef[i];
        }
        _dist=0.0f;
        _dists.setZero();
    }
    SelfCollEnergy(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,scalar K):_K(K) {
        for(_nrV=0; _nrV < 4 && c->_v[_nrV]; _nrV++) {
            _b[_nrV]=bc;
            _v[_nrV]=c->_v[_nrV];
            _dists[_nrV]=c->_v[_nrV]->_matDist;
        }
        _b[_nrV]=bv;
        _v[_nrV++]=v;
        _dist=calcCoef();
    }
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const {
        scalarD dist=_dist;
        for(sizeType v=0; v<_nrV; v++)
            dist+=_coef[v].dot(_v[v]->_pos.cast<scalarD>());
        for(sizeType vi=0; vi<_nrV; vi++) {
            if(f)f->block<3,1>(_v[vi]->_index*3+_b[vi]->_offset,0)-=_coef[vi]*(CF*dist*_K);
            if(H)for(sizeType vj=0; vj<_nrV; vj++)
                    addBlock(*H,_v[vi]->_index*3+_b[vi]->_offset,
                             _v[vj]->_index*3+_b[vj]->_offset,
                             _coef[vi]*_coef[vj].transpose()*(_K*CH));
        }
        return 0.5f*dist*dist*_K;
    }
    void debug() {
#define DELTA 1E-7f
        scalarD f0=calcCoef();
        for(sizeType i=0; i<_nrV; i++)
            for(sizeType d=0; d<3; d++) {
                scalar& val=_v[i]->_pos[d];
                scalar valT=val;
                val+=DELTA;
                scalarD f1=calcCoef();
                scalarD dfN=(f1-f0)/DELTA;
                scalarD df=_coef[i][d];
                INFOV("A: %f, N: %f, ERR: %f",df,dfN,std::abs(df-dfN))
                val=valT;
            }
#undef DELTA
    }
protected:
    scalarD calcCoef() {
        for(sizeType v=0; v<5; v++)
            _coef[v].setZero();
        if(_nrV == 5)return calcCoef3D();
        else return calcCoef2D();
    }
    scalarD calcCoef3D() {
        Mat3d M,invM;
        M.col(0)=(_v[0]->_pos-_v[3]->_pos).cast<scalarD>();
        M.col(1)=(_v[1]->_pos-_v[3]->_pos).cast<scalarD>();
        M.col(2)=(_v[2]->_pos-_v[3]->_pos).cast<scalarD>();
        invM=M.inverse();

        Vec3d L=invM.transpose()*Vec3d(_dists[0]-_dists[3],_dists[1]-_dists[3],_dists[2]-_dists[3]);
        Vec3d R=invM*(_v[4]->_pos-_v[3]->_pos).cast<scalarD>();
        _coef[0]=-L*R[0];
        _coef[1]=-L*R[1];
        _coef[2]=-L*R[2];
        _coef[3]=L*R[0]+L*R[1]+L*R[2]-L;
        _coef[4]=L;
        return _dists[3]+L.dot((_v[4]->_pos-_v[3]->_pos).cast<scalarD>());
    }
    scalarD calcCoef2D() {
        Mat2d M,invM;
        M.col(0)=(_v[0]->_pos-_v[2]->_pos).block<2,1>(0,0).cast<scalarD>();
        M.col(1)=(_v[1]->_pos-_v[2]->_pos).block<2,1>(0,0).cast<scalarD>();
        invM=M.inverse();

        Vec2d L=invM.transpose()*Vec2d(_dists[0]-_dists[2],_dists[1]-_dists[2]);
        Vec2d R=invM*(_v[3]->_pos-_v[2]->_pos).block<2,1>(0,0).cast<scalarD>();
        _coef[0].block<2,1>(0,0)=-L*R[0];
        _coef[1].block<2,1>(0,0)=-L*R[1];
        _coef[2].block<2,1>(0,0)=L*R[0]+L*R[1]-L;
        _coef[3].block<2,1>(0,0)=L;
        return _dists[2]+L.dot((_v[3]->_pos-_v[2]->_pos).block<2,1>(0,0).cast<scalarD>());
    }
    //data
    boost::shared_ptr<FEMBody> _b[5];
    boost::shared_ptr<FEMVertex> _v[5];
    Vec3d _coef[5];
    Vec4d _dists;
    scalarD _dist;
    sizeType _nrV;
    const scalar _K;
};
class GeomCollEnergy : public FEMEnergy
{
public:
    //Energy: 0.5f*K*(n.dot(x0-x))^2
    GeomCollEnergy(boost::shared_ptr<FEMVertex> v,const Vec3& p0,scalar K):_v(v),_p0(p0.cast<scalarD>()),_K(K) {
        _n=_p0-v->_pos.cast<scalarD>();
        _n/=std::max<scalarD>(_n.norm(),1E-3f);
        _H=_n*_n.transpose()*_K;
    }
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const {
        Vec3d dx=_v->_pos.cast<scalarD>()-_p0;
        if(dx.dot(_n) > 0.0f)return 0.0f;
        if(f)f->block<3,1>(_v->_index*3,0)-=CF*(_H*dx);
        if(H)addBlock(*H,_v->_index*3,_v->_index*3,_H*CH);
        return 0.5f*dx.dot(_H*dx);
    }
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const {
        Vec3d dx=U[_v->_index]-_p0;
        if(dx.dot(_n) > 0.0f)return 0.0f;
        Matd UI=U(_v->_index);
        if(f)(*f)-=UI.transpose()*CF*(_H*dx);
        if(H)(*H)+=UI.transpose()*(_H*UI).eval()*CH;
        return 0.5f*dx.dot(_H*dx);
    }
protected:
    boost::shared_ptr<FEMVertex> _v;
    Vec3d _p0,_n;
    Mat3d _H;
    const scalar _K;
};
DefaultFEMCollider::DefaultFEMCollider(const FEMMesh& mesh,scalar K,scalar CF,scalar CH)
    :_mesh(mesh),_K(K),_CF(CF),_CH(CH)
{
    for(sizeType i=0; i<_mesh.nrB(); i++)
        _mesh.getB(i)._system->clearEnergy(typeid(GeomCollEnergy).name());
    _fE.setZero(_mesh.nrV()*3);
}
const DefaultFEMCollider::Vec& DefaultFEMCollider::getFE() const
{
    return _fE;
}
void DefaultFEMCollider::getHE(Eigen::SparseMatrix<scalarD,0,sizeType>& H)
{
    H.resize(_mesh.nrV()*3,_mesh.nrV()*3);
    H.setFromTriplets(_HTrips.begin(),_HTrips.end());
}
void DefaultFEMCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV)
{
    for(sizeType i=0; i<nrV; i++)
        b[i]->_system->onCollisionSelf(v[i]->_index);
    SelfCollEnergy(b,v,coef,nrV,_K).eval(&_fE,&_HTrips,_CF,_CH);
}
void DefaultFEMCollider::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary)
{
    for(sizeType i=0; i<4; i++)
        if(c->_v[i])bc->_system->onCollisionSelf(c->_v[i]->_index);
    bv->_system->onCollisionSelf(v->_index);
    SelfCollEnergy(bc,c,bv,v,_K).eval(&_fE,&_HTrips,_CF,_CH);
}
void DefaultFEMCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n)
{
    b->_system->onCollision(v->_index);
    boost::shared_ptr<FEMEnergy> e(new GeomCollEnergy(v,v->_pos+n,_K));
    b->_system->addEnergy(e);
}