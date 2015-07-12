#include <iomanip>
#include "FEMMesh.h"
#include "FEMCollision.h"
#include "FEMSparseReducedBasis.h"
#include "FEMUtils.h"
#include "FEMRotationUtil.h"
#include "FEMBVHBuilder.h"
#include "FEMSystem.h"

#include "GridOp.h"
#include "ParallelPoissonDiskSampling.h"
#include "CollisionDetection.h"

#include <boost/unordered_set.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/interprocess/streams/vectorstream.hpp>

USE_PRJ_NAMESPACE

//FEMMESH
//Vertex
FEMVertex::FEMVertex():Serializable(0),_index(-1),_surface(false)
{
    _pos.setZero();
    _pos0.setZero();
    _mass=0.0f;
    _matDist=0.0f;
}
bool FEMVertex::read(std::istream& is)
{
    readBinaryData(_pos,is);
    readBinaryData(_pos0,is);
    readBinaryData(_mass,is);
    readBinaryData(_matDist,is);
    readBinaryData(_index,is);
    readBinaryData(_surface,is);
    return is.good();
}
bool FEMVertex::write(std::ostream& os) const
{
    writeBinaryData(_pos,os);
    writeBinaryData(_pos0,os);
    writeBinaryData(_mass,os);
    writeBinaryData(_matDist,os);
    writeBinaryData(_index,os);
    writeBinaryData(_surface,os);
    return os.good();
}
//FEMCell
FEMCell::FEMCell():Serializable(2),_index(-1)
{
    _d.setZero();
    _mass=0.0f;
}
Vec3 FEMCell::operator[](sizeType i) const
{
    return _v[i]->_pos;
}
sizeType FEMCell::operator()(sizeType i) const
{
    return _v[i]->_index;
}
Vec3 FEMCell::get(sizeType i,const Vec* off) const
{
    Vec3 ret=_v[i]->_pos;
    if(off)ret+=off->block<3,1>(_v[i]->_index*3,0).cast<scalar>();
    return ret;
}
void FEMCell::addParticle(boost::shared_ptr<FEMInterp> I)
{
    I->_next=_pSet;
    _pSet=I;
}
BBox<scalar> FEMCell::getBB() const
{
    BBox<scalar> bb;
    int nrV=(_v[3] == NULL ? 3 : 4);
    for(int i=0; i<nrV; i++)
        bb.setUnion(_v[i]->_pos);
    return bb;
}
bool FEMCell::contain(const Vec3& pos,FEMInterp& I) const
{
    const FEMCell& c=*this;
    if(_v[3] == NULL) {
        TriangleTpl<scalar> tri(c[0],c[1],c[2]);
        I._coef.block<3,1>(0,0)=tri.bary(pos);
        return tri.isInside(pos);
    } else {
        TetrahedronTpl<scalar> tet(c[0],c[1],c[2],c[3]);
        ASSERT(!tet._swap)
        I._coef=tet.bary(pos);
        return tet.isInside(pos);
    }
}
void FEMCell::findWeight(const Vec3& pos,Vec4& bary,scalar& sqrDist) const
{
    Vec3 cp;
    bary.setConstant(0.0f);
    const FEMCell& c=*this;
    if(_v[3] == NULL) {
        TriangleTpl<scalar> tri(c[0],c[1],c[2]);
        tri.calcPointDist(pos,sqrDist,cp,bary);
    } else {
        TetrahedronTpl<scalar> tet(c[0],c[1],c[2],c[3]);
        ASSERT(!tet._swap)
        tet.calcPointDist(pos,sqrDist,cp,bary);
    }
}
Vec3 FEMCell::getVert(const FEMInterp& p,Vec3* N) const
{
    Vec3 pos=Vec3::Zero();
    for(sizeType i=0; i<4; i++)
        if(_v[i])pos+=p._coef[i]*_v[i]->_pos;
    if(N) {
        Mat3 F,FN,FR;
        buildF(F,FN,FR);
        (*N)=FN*(*N);
        (*N)/=std::max((*N).norm(),ScalarUtil<scalar>::scalar_eps);
    }
    return pos;
}
void FEMCell::buildF(Mat3& F,Mat3& FN,Mat3& FR,Eigen::Matrix<scalarD,9,9>* DRDF,bool forcePositive) const
{
    //compute F
    const FEMCell& c=*this;
    {
        Mat3 DX=Mat3::Zero();
        DX.col(0)=c[1]-c[0];
        DX.col(1)=c[2]-c[0];
        if(_v[3] != NULL)
            DX.col(2)=c[3]-c[0];
        F=DX*_d.cast<scalar>();
    }
    //compute SVD and force positive if required
    Mat3 U,V;
    Vec3 S;
    if(forcePositive) {
        if(_v[3])adjustF3D(F,U,V,S);
        else adjustF2D(F,U,V,S);
    } else {
        Eigen::JacobiSVD<Mat3> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
        U=svd.matrixU();
        V=svd.matrixV();
        S=svd.singularValues();
    }
    //compute FR
    FR=U*V.transpose();
    //compute FN
    Mat3 invSigma=Mat3::Zero();
    invSigma.diagonal()=(1.0f/S.array()).matrix();
    if(_v[3] == NULL)invSigma(2,2)=0.0f;
    Mat3 invS=V*invSigma*V.transpose();
    FN=FR*invS;
    //compute DRDF
    if(DRDF)derivRDF(*DRDF,F,U,S,V);
}
void FEMCell::makePositive()
{
    const FEMCell& c=*this;
    if(_v[3] != NULL) {
        TetrahedronTpl<scalar> tet(c[0],c[1],c[2],c[3]);
        if(tet._swap)std::swap(_v[2],_v[3]);
    }
}
void FEMCell::buildD()
{
    if(_v[3] == NULL) {
        //2D
        _d.setZero();
        Mat2d d22;
        d22.col(0)=(_v[1]->_pos0-_v[0]->_pos0).block<2,1>(0,0).cast<scalarD>();
        d22.col(1)=(_v[2]->_pos0-_v[0]->_pos0).block<2,1>(0,0).cast<scalarD>();
        _d.block<2,2>(0,0)=d22.inverse();
    } else {
        //3D
        _d.col(0)=(_v[1]->_pos0-_v[0]->_pos0).cast<scalarD>();
        _d.col(1)=(_v[2]->_pos0-_v[0]->_pos0).cast<scalarD>();
        _d.col(2)=(_v[3]->_pos0-_v[0]->_pos0).cast<scalarD>();
        Mat3d invD=_d.inverse();
        _d=invD;
    }
}
bool FEMCell::read(std::istream& is,IOData* dat)
{
    readBinaryData(_v[0],is,dat);
    readBinaryData(_v[1],is,dat);
    readBinaryData(_v[2],is,dat);
    readBinaryData(_v[3],is,dat);
    readBinaryData(_d,is);
    readBinaryData(_mass,is);
    readBinaryData(_index,is);
    readBinaryData(_pSet,is,dat);
    return is.good();
}
bool FEMCell::write(std::ostream& os,IOData* dat) const
{
    writeBinaryData(_v[0],os,dat);
    writeBinaryData(_v[1],os,dat);
    writeBinaryData(_v[2],os,dat);
    writeBinaryData(_v[3],os,dat);
    writeBinaryData(_d,os);
    writeBinaryData(_mass,os);
    writeBinaryData(_index,os);
    writeBinaryData(_pSet,os,dat);
    return os.good();
}
//FEMInterp
FEMInterp::FEMInterp():Serializable(1),_id(-1) {}
FEMInterp::FEMInterp(sizeType id):Serializable(1),_id(id) {}
FEMInterp::FEMInterp(sizeType id,const Vec4& coef):Serializable(1),_id(id),_coef(coef) {}
bool FEMInterp::read(std::istream& is,IOData* dat)
{
    readBinaryData(_cell,is,dat);
    readBinaryData(_next,is,dat);
    readBinaryData(_id,is);
    readBinaryData(_coef,is);
    return is.good();
}
bool FEMInterp::write(std::ostream& os,IOData* dat) const
{
    writeBinaryData(_cell,os,dat);
    writeBinaryData(_next,os,dat);
    writeBinaryData(_id,os);
    writeBinaryData(_coef,os);
    return os.good();
}
//FEMBody
FEMBody::FEMBody():Serializable(4)
{
    _X0.setZero();
    _W.setZero();
}
sizeType FEMBody::nrSV() const
{
    sizeType nrSV=0;
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        if(_vss[i]->_surface)nrSV++;
    return nrSV;
}
sizeType FEMBody::dim() const
{
    return _css[0]->_v[3] ? 3 : 2;
}
void FEMBody::applyTrans(const Mat4& R,bool pos,bool pos0)
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        if(pos) _vss[i]->_pos =R.block<3,3>(0,0)*_vss[i]->_pos +R.block<3,1>(0,3);
        if(pos0)_vss[i]->_pos0=R.block<3,3>(0,0)*_vss[i]->_pos0+R.block<3,1>(0,3);
    }
    if(pos0 && _basis)_basis->applyTrans(R);
    if(pos0 && _system)_system->applyTrans(R);
    assemble();
    updateMesh();
}
void FEMBody::getFace
(vector<std::pair<Vec3i,Vec2i> >* surfaceFaces,
 vector<std::pair<Vec3i,Vec2i> >* internalFaces) const
{
    sizeType Dim=dim();
    boost::unordered_map<Vec3i,Vec3i,Hash> faceMap;
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        const FEMCell& c=*(_css[i]);
        for(sizeType f=0; f<Dim+1; f++) {
            Vec3i id=Vec3i::Constant(-1);
            for(sizeType j=0,k=0; j<Dim+1; j++)
                if(j!=f)id[k++]=c._v[j]->_index;
            std::sort(id.data(),id.data()+Dim);

            boost::unordered_map<Vec3i,Vec3i,Hash>::iterator it=faceMap.find(id);
            if(it == faceMap.end())
                faceMap[id]=Vec3i(c._index,-1,c._v[f]->_index);
            else it->second[1]=c._index;
        }
    }
    for(boost::unordered_map<Vec3i,Vec3i,Hash>::const_iterator
		  beg=faceMap.begin(),end=faceMap.end(); beg!=end; beg++)
	  if(surfaceFaces && beg->second[1] == -1)
		surfaceFaces->push_back(make_pair(beg->first,Vec2i(beg->second[0],beg->second[2])));
	  else if(internalFaces && beg->second[1] >= 0)
		internalFaces->push_back(make_pair(beg->first,beg->second.block<2,1>(0,0)));
}
bool FEMBody::read(std::istream& is,IOData* dat)
{
    readVector(_vss,is,dat);
    readVector(_css,is,dat);
    readVector(_evss,is,dat);
    _pSet.read(is);
    readBinaryData(_basis,is,dat);
    readBinaryData(_cubatureKinetic,is,dat);
    readBinaryData(_cubaturePotential,is,dat);
    readBinaryData(_S,is);
    readBinaryData(_X0,is);
    readBinaryData(_W,is);
    return is.good();
}
bool FEMBody::write(std::ostream& os,IOData* dat) const
{
    writeVector(_vss,os,dat);
    writeVector(_css,os,dat);
    writeVector(_evss,os,dat);
    _pSet.write(os);
    writeBinaryData(_basis,os,dat);
    writeBinaryData(_cubatureKinetic,os,dat);
    writeBinaryData(_cubaturePotential,os,dat);
    writeBinaryData(_S,os);
    writeBinaryData(_X0,os);
    writeBinaryData(_W,os);
    return os.good();
}
bool FEMBody::writeABQ(const std::string& path) const
{
    boost::filesystem::ofstream os(path);
    os << "*NODE" << std::endl;
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
	  os<< std::setprecision(12) << (i+1) << ", " <<
           _vss[i]->_pos[0] << ", " <<
           _vss[i]->_pos[1] << ", " <<
           _vss[i]->_pos[2] << std::endl;
    os << "*ELEMENT, type=C3D4, ELSET=PART1" << std::endl;
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        os << (i+1);
        for(char v=0; v<4; v++)
            if(_css[i]->_v[v])
                os << ', ' << (_css[i]->_v[v]->_index+1);
            else os << ", 0";
        os << std::endl;
    }
    if(_basis) {
        boost::filesystem::path pathB=boost::filesystem::path(path).replace_extension(".BASIS");
        boost::filesystem::ofstream osB(pathB,ios::binary);
        _basis->write(osB);
    }
    if(_cubatureKinetic) {
        boost::filesystem::path pathC=boost::filesystem::path(path).replace_extension(".CUBK");
        boost::filesystem::ofstream osC(pathC,ios::binary);
        _cubatureKinetic->write(osC);
    }
    if(_cubaturePotential) {
        boost::filesystem::path pathC=boost::filesystem::path(path).replace_extension(".CUBP");
        boost::filesystem::ofstream osC(pathC,ios::binary);
        _cubaturePotential->write(osC);
    }
    os << "*ELSET, ELSET=EALL, GENERATE" << std::endl;
    os << "1, " << _css.size() << std::endl;
    return os.good();
}
void FEMBody::writeObj(ObjMesh& sMesh) const
{
    vector<int>& vMap=sMesh.getIG();
    if(vMap.size() != nrV()) {
        vector<std::pair<Vec3i,Vec2i> > sFace;
        getFace(&sFace,NULL);

        int nrSV=0;
        vMap.assign(nrV(),-1);
        for(sizeType f=0; f<(sizeType)sFace.size(); f++)
            for(char v=0; v<dim(); v++) {
                if(vMap[sFace[f].first[v]] == -1)
                    vMap[sFace[f].first[v]]=nrSV++;
                sFace[f].first[v]=vMap[sFace[f].first[v]];
            }
        sMesh.getV().resize(nrSV);
        sMesh.getI().resize(sFace.size());
        for(sizeType i=0; i<(sizeType)sFace.size(); i++)
            sMesh.getI()[i]=sFace[i].first;
    }

    for(sizeType i=0; i<(sizeType)vMap.size(); i++)
        if(vMap[i] >= 0)
            sMesh.getV()[vMap[i]]=getV(i)._pos;

    if(dim() == 3) {
        sMesh.smooth();
        sMesh.makeUniform();
        sMesh.smooth();
    }
}
void FEMBody::writeObj(std::ostream& os) const
{
    ObjMesh mesh;
    writeObj(mesh);
    mesh.getIG().clear();
    mesh.write(os);
}
void FEMBody::writeVTK(VTKWriter<scalar>& os,const Vec* off,const Vec* off0,sizeType bid,Vec* customV,Vec* customC) const
{
    //vss
    vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    vector<scalar> cFEMBodyIds,vSurfIds,vDists,cV,cC;
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        Vec3 v;
        if(off0)v=off0->block<3,1>(_vss[i]->_index*3,0).cast<scalar>();
        else v=_vss[i]->_pos;
        if(off)v+=off->block<3,1>(_vss[i]->_index*3,0).cast<scalar>();
        vss.push_back(v);
        vSurfIds.push_back(_vss[i]->_surface ? 1.0f : 0.0f);
        vDists.push_back(_vss[i]->_matDist);
        if(customV)
            cV.push_back((scalar)(*customV)[i]);
    }
    //iss
    vector<Vec4i,Eigen::aligned_allocator<Vec4i> > iss;
    for(sizeType t=0; t<(sizeType)_css.size(); t++) {
        const FEMCell& c=*(_css[t]);
        iss.push_back(Vec4i(c._v[0]->_index,c._v[1]->_index,
                            c._v[2]->_index,c._v[3]?c._v[3]->_index:-1));
        cFEMBodyIds.push_back((scalar)bid);
        if(customC)
            cC.push_back((scalar)(*customC)[t]);
    }
    os.setRelativeIndex();
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(iss.begin(),iss.end(),dim() == 2 ? VTKWriter<scalar>::TRIANGLE : VTKWriter<scalar>::TETRA,true);
    os.appendCustomData("cellFEMBodyId",cFEMBodyIds.begin(),cFEMBodyIds.end());
    os.appendCustomPointData("vertSurfId",vSurfIds.begin(),vSurfIds.end());
    os.appendCustomPointData("vertDists",vDists.begin(),vDists.end());
    if(customV)
        os.appendCustomPointData("customV",cV.begin(),cV.end());
    if(customC)
        os.appendCustomData("customC",cC.begin(),cC.end());
}
void FEMBody::parityCheck() const
{
    //check vss and css index
    for(sizeType i=0; i<(sizeType)_css.size(); i++)
        ASSERT(_css[i]->_index == i)
        for(sizeType i=0; i<(sizeType)_vss.size(); i++)
            ASSERT(_vss[i]->_index == i)
            //check particle
            ASSERT((sizeType)_evss.size() == _pSet.size())
            for(sizeType i=0; i<_pSet.size(); i++) {
                ASSERT(_evss[i]->_cell == _css[_evss[i]->_cell->_index])
                ASSERT(_evss[i]->_id == i)
                scalar err=(_evss[i]->_cell->getVert(*(_evss[i]))-_pSet[i]._pos).norm();
                if(err > EPS)INFOV("Particle %d err: %f",i,err)
                }
    sizeType nrP=0;
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        boost::shared_ptr<FEMInterp> I=_css[i]->_pSet;
        while(I) {
            ASSERT(I == _evss[I->_id])
            ASSERT(I->_cell == _css[i])
            ASSERT(I->_next != I)
            I=I->_next;
            nrP++;
        }
    }
    ASSERT(nrP == _pSet.size())
    //check face
    sizeType nrS=nrSV();
    for(sizeType i=0; i<nrS; i++)
        ASSERT(_vss[i]->_surface)
        for(sizeType i=nrS; i<(sizeType)_vss.size(); i++)
            ASSERT(!_vss[i]->_surface);
}
FEMBody& FEMBody::operator=(const FEMBody& other)
{
    boost::interprocess::basic_ovectorstream<vector<char> > os(ios::binary);
    {
        IOData dat;
        dat.registerType<FEMVertex>();
        dat.registerType<FEMCell>();
        dat.registerType<FEMInterp>();
        dat.registerType(copy());
        dat.registerType<SparseReducedBasis>();
        dat.registerType<Cubature>();
        other.write(os,&dat);
    }
    {
        boost::interprocess::basic_ivectorstream<vector<char> > is(os.vector(),ios::binary);
        IOData dat;
        dat.registerType<FEMVertex>();
        dat.registerType<FEMCell>();
        dat.registerType<FEMInterp>();
        dat.registerType(copy());
        dat.registerType<SparseReducedBasis>();
        dat.registerType<Cubature>();
        read(is,&dat);
    }
    _tree=other._tree;
    _system.reset((FEMSystem*)NULL);
    if(other._system)_system=other._system->copy(*this);
    parityCheck();
    return *this;
}
void FEMBody::assemble()
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        FEMVertex& v=*(_vss[i]);
        v._index=i;
        v._mass=0.0f;
        v._pos=v._pos0;
    }
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        FEMCell& c=*(_css[i]);
        c._index=i;
        c.makePositive();

        //build mass
        if(dim() == 2)
            c._mass=TriangleTpl<scalar>(c[0],c[1],c[2]).area();
        else c._mass=TetrahedronTpl<scalar>(c[0],c[1],c[2],c[3]).volume();

        //build deformation gradient
        c.buildD();
    }
    for(sizeType i=0; i<(sizeType)_css.size(); i++)
        for(sizeType j=0; j<(dim() == 3 ? 4 : 3); j++)
            _css[i]->_v[j]->_mass+=_css[i]->_mass/(dim() == 3 ? 4.0f : 3.0f);
}
void FEMBody::updateMesh()
{
    Mat3 F,FN,FR;
    sizeType nrP=_pSet.size();
    OMP_PARALLEL_FOR_I(OMP_PRI(F,FN,FR))
    for(sizeType i=0; i<nrP; i++) {
        const FEMInterp& I=*(_evss[i]);
        const FEMCell& c=*(_evss[i]->_cell);
        ParticleN<scalar>& p=_pSet[I._id];
        p._pos.setZero();
        for(int i=0; i<dim()+1; i++)
            p._pos+=I._coef[i]*c[i];
        c.buildF(F,FN,FR,NULL,false);
        p._normal=FN*p._normal;
        p._normal/=std::max(p._normal.norm(),ScalarUtil<scalar>::scalar_eps);
    }
}
void FEMBody::buildG(vector<Eigen::Triplet<scalarD,sizeType> >& trips) const
{
    Eigen::Matrix<scalarD,9,12> GComp;
    for(sizeType i=0; i<nrC(); i++) {
        const FEMCell& C=getC(i);
        if(dim() == 2)calcGComp2D(GComp,C._d);
        else calcGComp3D(GComp,C._d);
        for(sizeType r=0; r<9; r++)
            for(sizeType c=0; c<(dim()+1)*3; c++)
                if(GComp(r,c) != 0.0f)
                    trips.push_back(Eigen::Triplet<scalarD,sizeType>(r+i*9,C._v[c/3]->_index*3+c%3,GComp(r,c)));
    }
}
//getter
void FEMBody::getPos(Vec& X) const
{
    sizeType nrV=(sizeType)_vss.size();
    if(X.size() < nrV*3)
        X.resize(nrV*3);
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        X.block<3,1>(i*3,0)=_vss[i]->_pos.cast<scalarD>();
}
void FEMBody::getDPos(Vec& X) const
{
    sizeType nrV=(sizeType)_vss.size();
    if(X.size() < nrV*3)
        X.resize(nrV*3);
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        X.block<3,1>(i*3,0)=(_vss[i]->_pos-_vss[i]->_pos0).cast<scalarD>();
}
void FEMBody::getMass(Vec& M,bool surface) const
{
    sizeType nrV=(sizeType)_vss.size();
    if(M.size() < nrV*3)
        M.resize(nrV*3);
    M.block(0,0,nrV*3,1).setZero();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        if(!surface||_vss[i]->_surface)
            M.block<3,1>(i*3,0).setConstant(_vss[i]->_mass);
}
void FEMBody::getMass(Eigen::DiagonalMatrix<scalarD,-1>& M,bool surface) const
{
    Vec D;
    getMass(D,surface);
    M.resize(D.size());
    M.diagonal()=D;
}
void FEMBody::getMass(Eigen::SparseMatrix<scalarD,0,sizeType>& M,bool surface) const
{
    Vec D;
    getMass(D,surface);
    vector<Eigen::Triplet<scalarD,sizeType> > trips;
    addI(trips,0,0,D);

    M.resize(nrV()*3,nrV()*3);
    M.setFromTriplets(trips.begin(),trips.end());
}
//setter
void FEMBody::setPos(const Vec& X)
{
    sizeType nrV=(sizeType)_vss.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        _vss[i]->_pos=X.block<3,1>(i*3,0).cast<scalar>();
}
void FEMBody::setDPos(const Vec& X)
{
    sizeType nrV=(sizeType)_vss.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        _vss[i]->_pos=X.block<3,1>(i*3,0).cast<scalar>()+_vss[i]->_pos0;
}
void FEMBody::clearPos(char dim)
{
    sizeType nrV=(sizeType)_vss.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        for(char d=0; d<3; d++)
            if(dim == -1 || d == dim)
                _vss[i]->_pos[d]=_vss[i]->_pos0[d];
}

//add geometry: FEMVertex
boost::shared_ptr<FEMVertex> addVertex(const Vec3i& base,
                                       boost::unordered_map<Vec3i,boost::shared_ptr<FEMVertex>,Hash >& vmap,
                                       vector<boost::shared_ptr<FEMVertex> >& vss,scalar cellSz)
{
    if(vmap.find(base) == vmap.end()) {
        boost::shared_ptr<FEMVertex> v(new FEMVertex);
        v->_pos0=v->_pos=base.cast<scalar>()*cellSz;
        vss.push_back(v);
        vmap[base]=v;
        return v;
    } else return vmap[base];
}
//add geometry: FEMTet
boost::shared_ptr<FEMCell> addTet
(const Vec3i& v0,const Vec3i& v1,const Vec3i& v2,const Vec3i& v3,
 boost::unordered_map<Vec3i,vector<boost::shared_ptr<FEMCell> >,Hash >& cmap,
 boost::unordered_map<Vec3i,boost::shared_ptr<FEMVertex>,Hash >& vmap,
 vector<boost::shared_ptr<FEMVertex> >& vss,
 vector<boost::shared_ptr<FEMCell> >& css,
 scalar cellSz,sizeType dim)
{
    //build topology
    boost::shared_ptr<FEMCell> c(new FEMCell);
    c->_v[0]=addVertex(v0,vmap,vss,cellSz);
    c->_v[1]=addVertex(v1,vmap,vss,cellSz);
    c->_v[2]=addVertex(v2,vmap,vss,cellSz);
    if(dim == 3)
        c->_v[3]=addVertex(v3,vmap,vss,cellSz);
    c->makePositive();
    css.push_back(c);
    return c;
}
//add geometry: FEMCell
vector<boost::shared_ptr<FEMCell> > addCell(const Vec3i& base,
        boost::unordered_map<Vec3i,vector<boost::shared_ptr<FEMCell> >,Hash >& cmap,
        boost::unordered_map<Vec3i,boost::shared_ptr<FEMVertex>,Hash >& vmap,
        vector<boost::shared_ptr<FEMVertex> >& vss,
        vector<boost::shared_ptr<FEMCell> >& css,
        scalar cellSz,sizeType dim)
{
    if(cmap.find(base) != cmap.end())
        return cmap.find(base)->second;

    vector<boost::shared_ptr<FEMCell> > ret;
    if(dim == 3) {
        ret.push_back(addTet(base+Vec3i(0,0,0),
                             base+Vec3i(1,0,0),
                             base+Vec3i(0,0,1),
                             base+Vec3i(0,1,0),cmap,vmap, vss,css,cellSz,dim));
        ret.push_back(addTet(base+Vec3i(1,1,1),
                             base+Vec3i(0,1,1),
                             base+Vec3i(1,1,0),
                             base+Vec3i(1,0,1),cmap,vmap, vss,css,cellSz,dim));

        ret.push_back(addTet(base+Vec3i(1,0,1),
                             base+Vec3i(1,0,0),
                             base+Vec3i(0,0,1),
                             base+Vec3i(0,1,0),cmap,vmap, vss,css,cellSz,dim));
        ret.push_back(addTet(base+Vec3i(0,1,0),
                             base+Vec3i(0,1,1),
                             base+Vec3i(1,1,0),
                             base+Vec3i(1,0,1),cmap,vmap, vss,css,cellSz,dim));

        ret.push_back(addTet(base+Vec3i(0,0,1),
                             base+Vec3i(0,1,0),
                             base+Vec3i(0,1,1),
                             base+Vec3i(1,0,1),cmap,vmap, vss,css,cellSz,dim));
        ret.push_back(addTet(base+Vec3i(0,1,0),
                             base+Vec3i(1,0,0),
                             base+Vec3i(1,1,0),
                             base+Vec3i(1,0,1),cmap,vmap, vss,css,cellSz,dim));
    } else {
        ret.push_back(addTet(base+Vec3i(0,0,0),
                             base+Vec3i(1,0,0),
                             base+Vec3i(1,1,0),
                             base,cmap,vmap, vss,css,cellSz,dim));
        ret.push_back(addTet(base+Vec3i(0,0,0),
                             base+Vec3i(1,1,0),
                             base+Vec3i(0,1,0),
                             base,cmap,vmap, vss,css,cellSz,dim));
    }
    cmap[base]=ret;
    return ret;
}

//method
FEMMesh::FEMMesh(sizeType dim,boost::shared_ptr<FEMCollision> coll)
    :HasMagic(0xabcddcba),_dim(dim),_coll(coll)
{
    _coll->setMesh(*this);
}
FEMMesh::FEMMesh(const FEMMesh& other):HasMagic(0xabcddcba),_dim(other._dim)
{
    _coll=other._coll->copy();
    _coll->setMesh(*this);
    operator=(other);
}
void FEMMesh::updateMesh(scalar expand)
{
    for(sizeType i=0; i<(sizeType)_bss.size(); i++)
        _bss[i]->updateMesh();
    _coll->updateMesh();
}
void FEMMesh::reset(const std::string& path,scalar rad)
{
    //read data
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > nodes;
    std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> > tets;
    {
        boost::filesystem::ifstream is(path);
        sizeType id;
        char comma;
        string line;
        getline(is,line);
        while(is.good()) {
            if(beginsWith(line,"*NODE")) {
                while(getline(is,line).good() && !line.empty() && line[0] != '*') {
                    scalar x,y,z;
                    istringstream iss(line);
                    iss >> id >> comma >> x >> comma >> y >> comma >> z;
                    nodes.push_back(Vec3(x,y,z));
                }
            } else if(beginsWith(line,"*ELEMENT")) {
                while(getline(is,line).good() && !line.empty() && line[0] != '*') {
                    sizeType x,y,z,w;
                    istringstream iss(line);
                    iss >> id >> comma >> x >> comma >> y >> comma >> z >> comma >> w;
                    tets.push_back(Vec4i(x-1,y-1,z-1,w-1));
                }
            } else getline(is,line);
        }
    }
    //setup mesh
    ParticleSetN pSet;
    if(rad > 0.0f) {
        ParallelPoissonDiskSampling sampler(_dim);
        sampler.getRadius()=rad;
        sampler.sample(nodes,tets);
        pSet=sampler.getPSet();
    } else pSet.clear();
    bool inv=reset(nodes,tets,pSet);
    updateMesh();
    //setup basis
    if(inv && nrB() == 1) {
        boost::filesystem::path pathB=boost::filesystem::path(path).replace_extension(".BASIS");
        if(boost::filesystem::exists(pathB)) {
            boost::filesystem::ifstream isB(pathB,ios::binary);
            getB(0)._basis.reset(new SparseReducedBasis);
            getB(0)._basis->read(isB);
        }
        boost::filesystem::path pathC=boost::filesystem::path(path).replace_extension(".CUBK");
        if(boost::filesystem::exists(pathC)) {
            boost::filesystem::ifstream isC(pathC,ios::binary);
            getB(0)._cubatureKinetic.reset(new Cubature);
            getB(0)._cubatureKinetic->read(isC);
        }
        pathC=boost::filesystem::path(path).replace_extension(".CUBP");
        if(boost::filesystem::exists(pathC)) {
            boost::filesystem::ifstream isC(pathC,ios::binary);
            getB(0)._cubaturePotential.reset(new Cubature);
            getB(0)._cubaturePotential->read(isC);
        }
    }
    //write invariant version
    if(!inv)
        getB(0).writeABQ(path);
}
void FEMMesh::reset(BBox<scalar> bb,const ImplicitFunc<scalar>& f,scalar rad,scalar cellSz)
{
    if(_dim == 2)
        bb._maxC[2]=bb._minC[2]=0.0f;
    Vec3i nrFEMCell=ceil((Vec3)(bb.getExtent()/std::abs(rad)));
    ScalarField phi;
    phi.reset(nrFEMCell,bb,0.0f,true);
    GridOp<scalar,scalar>::copyFromImplictFunc(phi,f);

    ParallelPoissonDiskSampling sampler(_dim);
    sampler.getRadius()=std::abs(rad);
    sampler.sample(phi,true);
    ParticleSetN pSet=sampler.getPSet();
    for(sizeType i=0; i<pSet.size(); i++) {
        pSet[i]._normal=phi.sampleSafeGrad(pSet[i]._pos);
        pSet[i]._normal/=std::max(pSet[i]._normal.norm(),ScalarUtil<scalar>::scalar_eps);
    }
    reset(cellSz,pSet,rad>0.0f);
    updateMesh();
}
bool FEMMesh::reset(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& nodes,
                    const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& tets,const ParticleSetN& pset)
{
    _vss.resize(nodes.size());
    for(sizeType i=0; i<(sizeType)nodes.size(); i++) {
        _vss[i].reset(new FEMVertex);
        _vss[i]->_pos=_vss[i]->_pos0=nodes[i];
    }
    _css.resize(tets.size());
    for(sizeType i=0; i<(sizeType)tets.size(); i++) {
        _css[i].reset(new FEMCell);
        _css[i]->_v[0]=_vss[tets[i][0]];
        _css[i]->_v[1]=_vss[tets[i][1]];
        _css[i]->_v[2]=_vss[tets[i][2]];
        if(tets[i][3] == -1) {
            _dim=2;
        } else {
            _dim=3;
            _css[i]->_v[3]=_vss[tets[i][3]];
        }
    }
    _pSet=pset;
    _cellSz=0.0f;
    return assemble();
}
void FEMMesh::reset(scalar cellSz,const ParticleSetN& pset,bool usePSet)
{
    _vss.clear();
    _css.clear();
    _cellSz=cellSz;

    boost::unordered_map<Vec3i,vector<boost::shared_ptr<FEMCell> >,Hash > cmap;
    boost::unordered_map<Vec3i,boost::shared_ptr<FEMVertex>,Hash > vmap;
    for(sizeType i=0; i<pset.size(); i++)
        addCell(floor((Vec3)(pset[i]._pos/cellSz)),cmap,vmap,
                _vss,_css,_cellSz,_dim);
    if(usePSet)
        _pSet=pset;
    assemble();
}
void FEMMesh::applyTrans(const Mat4& R,sizeType bid,bool pos,bool pos0)
{
    for(sizeType i=0; i<(sizeType)_bss.size(); i++)
        if(i == bid || bid == -1)_bss[i]->applyTrans(R,pos,pos0);
}
FEMCollision& FEMMesh::getColl()
{
    return *_coll;
}
void FEMMesh::buildOffset()
{
    sizeType off=0;
    for(sizeType i=0; i<nrB(); i++) {
        _bss[i]->_offset=off;
        off+=_bss[i]->nrV()*3;
    }
}
sizeType FEMMesh::nrV() const
{
    sizeType nrV=0;
    for(sizeType i=0; i<nrB(); i++)
        nrV+=_bss[i]->nrV();
    return nrV;
}
//body operation
FEMMesh& FEMMesh::operator=(const FEMMesh& other)
{
    ASSERT_MSG(_dim == other._dim,"Copying mesh of different dimension!");
    if(this == &other)
        return *this;
    _bss.clear();
    return operator+=(other);
}
FEMMesh& FEMMesh::operator+=(const FEMMesh& other)
{
    if(this == &other) {
        FEMMesh tmp=other;
        return operator+=(tmp);
    }
    sizeType off=(sizeType)_bss.size();
    _bss.resize(off+other._bss.size());
    for(sizeType i=off; i<(sizeType)_bss.size(); i++) {
        _bss[i]=_coll->createBody();
        *(_bss[i])=*(other._bss[i-off]);
    }
    updateMesh();
    return *this;
}
FEMMesh FEMMesh::operator+(const FEMMesh& other) const
{
    FEMMesh tmp=*this;
    tmp+=other;
    return tmp;
}
FEMMesh& FEMMesh::operator-=(sizeType bid)
{
    _bss.erase(_bss.begin()+bid);
    updateMesh();
    return *this;
}
FEMMesh FEMMesh::operator-(sizeType bid) const
{
    FEMMesh tmp=*this;
    tmp-=bid;
    return tmp;
}
struct BaryCallback {
    BaryCallback(const vector<boost::shared_ptr<FEMCell> >& css,const Vec3& pos,FEMInterp& it)
        :_css(css),_pos(pos),_it(it) {}
    void updateDist(const Node<sizeType>& n,scalar& sqrDist) {
        FEMInterp ITmp;
        scalar sqrDistTmp;
        _css[n._cell]->findWeight(_pos,ITmp._coef,sqrDistTmp);
        if(sqrDistTmp < sqrDist) {
            sqrDist=sqrDistTmp;
            _it._coef=ITmp._coef;
            _it._id=n._cell;
        }
    }
    const vector<boost::shared_ptr<FEMCell> >& _css;
    const Vec3 _pos;
    FEMInterp& _it;
};
void FEMMesh::getBary(const vector<Vec3,Eigen::aligned_allocator<Vec3> >& ps,std::vector<FEMInterp>& evss,scalar sqrDistThres) const
{
    //assemble body
    vector<sizeType> bodyId;
    vector<boost::shared_ptr<FEMCell> > cells;
    for(sizeType i=0; i<(sizeType)_bss.size(); i++) {
        vector<boost::shared_ptr<FEMCell> >& css=_bss[i]->_css;
        for(sizeType c=0; c<(sizeType)css.size(); c++) {
            cells.push_back(css[c]);
            bodyId.push_back(i);
        }
    }

    //build BVH
    sizeType nrC=(sizeType)cells.size();
    std::vector<Node<sizeType> > bvh(nrC);
    for(sizeType i=0; i<nrC; i++) {
        bvh[i]._bb=cells[i]->getBB();
        bvh[i]._cell=i;
        bvh[i]._nrCell=1;
    }
    buildBVH<sizeType>(bvh,_dim,-1);

    //query each point
    ParticleSetN pSetErr;
    BVHQuery<sizeType> query(bvh,_dim,-1);
    for(sizeType i=0; i<(sizeType)ps.size(); i++) {
        BaryCallback cb(cells,ps[i],evss[i]);
        scalar eps=1E-3f,sqrDist;
        while(true) {
            sqrDist=numeric_limits<scalar>::max();
            query.pointQuery(ps[i],eps,cb,sqrDist);
            if(sqrDist < numeric_limits<scalar>::max()) {
                evss[i]._cell=cells[evss[i]._id];
                evss[i]._id=bodyId[evss[i]._id];
                break;
            }
            eps*=10.0f;
        }
        if(sqrDist != 0.0f) {
            WARNINGV("Particle %lu not within control mesh, distance=%f",i,sqrt(sqrDist));
            ParticleN<scalar> p;
            p._pos=ps[i];
            pSetErr.addParticle(p);
        }
        ASSERT_MSGV(sqrDist < sqrDistThres,"Particle %lu not in control mesh!",i)
    }
    pSetErr.writeVTK("./err.vtk");
}
//IO
void FEMMesh::writeVTK(const std::string& path) const
{
    VTKWriter<scalar> os("Mesh",path,true);
    for(sizeType i=0; i<(sizeType)_bss.size(); i++)
        _bss[i]->writeVTK(os,NULL,NULL,i);
}
void FEMMesh::writePSetVTK(const std::string& path) const
{
    VTKWriter<scalar> os("PSet",path,true);
    for(sizeType i=0; i<(sizeType)_bss.size(); i++)
        _bss[i]->getPSet().writeVTK(os);
}
bool FEMMesh::write(std::ostream& os) const
{
    if(!HasMagic::writeMagic(os))
        return false;

    IOData dat;
    dat.registerType<FEMVertex>();
    dat.registerType<FEMCell>();
    dat.registerType<FEMInterp>();
    dat.registerType(boost::dynamic_pointer_cast<Serializable>(_coll->createBody()));
    dat.registerType<SparseReducedBasis>();
    dat.registerType<Cubature>();

    writeBinaryData(_cellSz,os,&dat);
    writeBinaryData(_dim,os,&dat);
    writeVector(_bss,os,&dat);
    return os.good();
}
bool FEMMesh::read(std::istream& is)
{
    _bss.clear();
    if(!HasMagic::readMagic(is))
        return false;

    IOData dat;
    dat.registerType<FEMVertex>();
    dat.registerType<FEMCell>();
    dat.registerType<FEMInterp>();
    dat.registerType(boost::dynamic_pointer_cast<Serializable>(_coll->createBody()));
    dat.registerType<SparseReducedBasis>();
    dat.registerType<Cubature>();

    readBinaryData(_cellSz,is,&dat);
    readBinaryData(_dim,is,&dat);
    readVector(_bss,is,&dat);
    for(sizeType i=0; i<(sizeType)_bss.size(); i++)
        _bss[i]->parityCheck();
    updateMesh();
    return is.good();
}
//private
void FEMMesh::buildBody()
{
    buildIndex();
    //build body id using disjoint set
    vector<sizeType> bodyId(_vss.size(),-1);
    for(sizeType i=0; i<(sizeType)_css.size(); i++)
        for(sizeType v=1; v<4; v++)
            if(_css[i]->_v[v]) {
                sizeType ri=_css[i]->_v[0]->_index,riTmp=ri;
                while(bodyId[ri] >= 0)ri=bodyId[ri];
                //path compression
                if(riTmp != ri)bodyId[riTmp]=ri;

                sizeType rj=_css[i]->_v[v]->_index,rjTmp=rj;
                while(bodyId[rj] >= 0)rj=bodyId[rj];
                //path compression
                if(rjTmp != rj)bodyId[rjTmp]=rj;

                if(ri != rj)
                    bodyId[rj]=ri;
            }
    //assign body id
    typedef boost::unordered_map<boost::shared_ptr<FEMVertex>,boost::shared_ptr<FEMBody> > BMAP;
    BMAP bSet;
    _bss.clear();
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        sizeType ri=_vss[i]->_index;
        while(bodyId[ri] >= 0)ri=bodyId[ri];
        BMAP::iterator it=bSet.find(_vss[ri]);
        if(it == bSet.end()) {
            bSet[_vss[ri]]=_coll->createBody();
            it=bSet.find(_vss[ri]);
            _bss.push_back(it->second);
        }
        it->second->_vss.push_back(_vss[i]);
    }
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        sizeType ri=_css[i]->_v[0]->_index;
        while(bodyId[ri] >= 0)ri=bodyId[ri];
        bSet[_vss[ri]]->_css.push_back(_css[i]);
    }
    //reorder body
    buildIndex();
}
void FEMMesh::buildFace()
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_surface=false;
    //put face vertices first
    for(sizeType i=0; i<(sizeType)_bss.size(); i++) {
        std::vector<std::pair<Vec3i,Vec2i> > faces;
        _bss[i]->getFace(&faces,NULL);
        for(sizeType f=0; f<(sizeType)faces.size(); f++)
            for(sizeType j=0; j<3; j++)
                if(faces[f].first[j] >= 0)
                    _vss[faces[f].first[j]]->_surface=true;

        sizeType nrSV=0;
        std::vector<boost::shared_ptr<FEMVertex> >& vss=_bss[i]->_vss;
        for(sizeType v=0; v<(sizeType)vss.size(); v++)
            if(vss[v]->_surface)
                std::swap(vss[nrSV++],vss[v]);
    }
    //reorder according to body
    sizeType nrVTmp=(sizeType)_vss.size();
    sizeType nrCTmp=(sizeType)_css.size();
    _vss.clear();
    _css.clear();
    for(sizeType i=0; i<(sizeType)_bss.size(); i++) {
        FEMBody& body=*(_bss[i]);
        _vss.insert(_vss.end(),body._vss.begin(),body._vss.end());
        _css.insert(_css.end(),body._css.begin(),body._css.end());
    }
    ASSERT(nrVTmp==(sizeType)_vss.size())
    ASSERT(nrCTmp==(sizeType)_css.size())
    buildIndex();
}
void FEMMesh::buildEvss()
{
    if(_pSet.size() == 0)return;
    vector<FEMInterp> evss(_pSet.size());
    vector<Vec3,Eigen::aligned_allocator<Vec3> > ps(_pSet.size());
    for(sizeType i=0; i<_pSet.size(); i++)
        ps[i]=_pSet[i]._pos;
    getBary(ps,evss,1E-6f);
    for(sizeType i=0; i<_pSet.size(); i++) {
        boost::shared_ptr<FEMInterp> I(new FEMInterp(evss[i]));
        I->_id=_bss[evss[i]._id]->_pSet.size();

        _bss[evss[i]._id]->_evss.push_back(I);
        _bss[evss[i]._id]->_pSet.addParticle(_pSet[i]);
        I->_cell->addParticle(I);
    }
}
void FEMMesh::buildIndex()
{
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_index=i;
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_css.size(); i++)
        _css[i]->_index=i;
}
bool FEMMesh::assemble()
{
    //save order
    vector<boost::shared_ptr<FEMVertex> > vssTmp=_vss;
    vector<boost::shared_ptr<FEMCell> > cssTmp=_css;
    //assemble
    buildBody();
    buildFace();
    buildEvss();
    for(sizeType i=0; i<(sizeType)_bss.size(); i++) {
        _bss[i]->assemble();
        _bss[i]->parityCheck();
    }
    updateMesh();
    //check order
    bool inv=true;
    for(sizeType i=0; i<(sizeType)vssTmp.size(); i++)
        if(vssTmp[i] != _vss[i])inv=false;
    for(sizeType i=0; i<(sizeType)cssTmp.size(); i++)
        if(cssTmp[i] != _css[i])inv=false;
    //clear
    _vss.clear();
    _css.clear();
    _pSet.clear();
    return inv;
}
