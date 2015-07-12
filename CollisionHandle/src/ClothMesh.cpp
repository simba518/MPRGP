#include "ClothMesh.h"
#include <boost/interprocess/streams/vectorstream.hpp>

USE_PRJ_NAMESPACE

ClothMesh::ClothVertex::ClothVertex():Serializable(1) {}
ClothMesh::ClothVertex::ClothVertex(const Vec3d& pos,MESH_TYPE type):Serializable(1)
{
    _pos=pos;
    _lastPos=pos;
    _vel.setZero();
    _pos0=_pos;
    _mass=0.0f;
    _type=type;
}
bool ClothMesh::ClothVertex::write(std::ostream& os,IOData* dat) const
{
    writeBinaryData(_pos,os,dat);
    writeBinaryData(_lastPos,os,dat);
    writeBinaryData(_vel,os,dat);
    writeBinaryData(_pos0,os,dat);
    writeBinaryData(_mass,os,dat);
    writeBinaryData(_index,os,dat);
    writeBinaryData(_weight,os,dat);
    writeBinaryData(_type,os,dat);
    //writeVector(_oneRing,os,dat);
    return os.good();
}
bool ClothMesh::ClothVertex::read(std::istream& is,IOData* dat)
{
    readBinaryData(_pos,is,dat);
    readBinaryData(_lastPos,is,dat);
    readBinaryData(_vel,is,dat);
    readBinaryData(_pos0,is,dat);
    readBinaryData(_mass,is,dat);
    readBinaryData(_index,is,dat);
    readBinaryData(_weight,is,dat);
    readBinaryData(_type,is,dat);
    //readVector(_oneRing,is,dat);
    return is.good();
}

ClothMesh::ClothEdge::ClothEdge():Serializable(2) {}
ClothMesh::ClothEdge::ClothEdge(boost::shared_ptr<ClothVertex> v0,boost::shared_ptr<ClothVertex> v1,MESH_TYPE type):Serializable(2)
{
    _v[0]=v0;
    _v[1]=v1;
    _pos=(_v[0]->_pos+_v[1]->_pos)*0.5f;
    _vel.setZero();
    _pos0=_pos;
    _mass=0.0f;
    _type=type;
}
bool ClothMesh::ClothEdge::write(std::ostream& os,IOData* dat) const
{
    writeBinaryData(_v[0],os,dat);
    writeBinaryData(_v[1],os,dat);
    //writeBinaryData(_t[0],os,dat);
    //writeBinaryData(_t[1],os,dat);
    writeBinaryData(_pos,os,dat);
    writeBinaryData(_vel,os,dat);
    writeBinaryData(_pos0,os,dat);
    writeBinaryData(_mass,os,dat);
    writeBinaryData(_index,os,dat);
    writeBinaryData(_type,os,dat);
    return os.good();
}
bool ClothMesh::ClothEdge::read(std::istream& is,IOData* dat)
{
    readBinaryData(_v[0],is,dat);
    readBinaryData(_v[1],is,dat);
    //readBinaryData(_t[0],is,dat);
    //readBinaryData(_t[1],is,dat);
    readBinaryData(_pos,is,dat);
    readBinaryData(_vel,is,dat);
    readBinaryData(_pos0,is,dat);
    readBinaryData(_mass,is,dat);
    readBinaryData(_index,is,dat);
    readBinaryData(_type,is,dat);
    return is.good();
}

ClothMesh::ClothTriangle::ClothTriangle():Serializable(3) {}
ClothMesh::ClothTriangle::ClothTriangle(boost::shared_ptr<ClothEdge> e0,bool p0,boost::shared_ptr<ClothEdge> e1,bool p1,boost::shared_ptr<ClothEdge> e2,bool p2,MESH_TYPE type):Serializable(3),_edgeDir(p0,p1,p2)
{
    _e[0]=e0;
    _e[1]=e1;
    _e[2]=e2;
    _type=type;
}
boost::shared_ptr<ClothMesh::ClothVertex> ClothMesh::ClothTriangle::getV0() const
{
    return _edgeDir[0] ? _e[0]->_v[0] : _e[0]->_v[1];
}
boost::shared_ptr<ClothMesh::ClothVertex> ClothMesh::ClothTriangle::getV1() const
{
    return _edgeDir[1] ? _e[1]->_v[0] : _e[1]->_v[1];
}
boost::shared_ptr<ClothMesh::ClothVertex> ClothMesh::ClothTriangle::getV2() const
{
    return _edgeDir[2] ? _e[2]->_v[0] : _e[2]->_v[1];
}
bool ClothMesh::ClothTriangle::write(std::ostream& os,IOData* dat) const
{
    // writeBinaryData(_e[0],os,dat);
    // writeBinaryData(_e[1],os,dat);
    // writeBinaryData(_e[2],os,dat);
    // writeBinaryData(_edgeDir,os,dat);
    // writeBinaryData(_index,os,dat);
    // writeBinaryData(_type,os,dat);
    return os.good();
}
bool ClothMesh::ClothTriangle::read(std::istream& is,IOData* dat)
{
    // readBinaryData(_e[0],is,dat);
    // readBinaryData(_e[1],is,dat);
    // readBinaryData(_e[2],is,dat);
    // readBinaryData(_edgeDir,is,dat);
    // readBinaryData(_index,is,dat);
    // readBinaryData(_type,is,dat);
    return is.good();
}

ClothMesh::ClothMesh():Serializable(0) {}
ClothMesh::ClothMesh(const ObjMeshD& mesh,MESH_TYPE type):Serializable(0)
{
    reset(mesh,type);
}
ClothMesh::ClothMesh(const ClothMesh& other):Serializable(0)
{
    operator=(other);
}
ClothMesh& ClothMesh::operator=(const ClothMesh& other)
{
    boost::interprocess::basic_ovectorstream<vector<char> > os(ios::binary);
    other.write(os);
    boost::interprocess::basic_ivectorstream<vector<char> > is(os.vector(),ios::binary);
    read(is);
    assembleIndex();
    assembleA();
    return *this;
}
void ClothMesh::reset(const ObjMeshD& mesh,MESH_TYPE type)
{
    //vertex list
    _vss.resize(mesh.getV().size());
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++)
        _vss[i].reset(new ClothVertex(mesh.getV()[i],type));

    //build edge map
    typedef map<pair<int,int>,ObjMeshD::Edge,ObjMeshD::EdgeMap::LSS> EMAP;
    ObjMeshD::EdgeMap eMap;
    mesh.buildEdge(eMap);

    //edge list
    std::map<pair<int,int>,sizeType> eIdMap;
    sizeType ie=0;
    _ess.resize(eMap._ess.size());
    for(EMAP::const_iterator beg=eMap._ess.begin(),end=eMap._ess.end();beg!=end;beg++,ie++){
        _ess[ie].reset(new ClothEdge(_vss[beg->first.first],_vss[beg->first.second],type));
        eIdMap[beg->first]=ie;
    }

    //triangle list
    _tss.resize(mesh.getI().size());
    for(sizeType i=0; i<(sizeType)mesh.getI().size(); i++) {
        const Vec3i& it=mesh.getI()[i];
        pair<int,int> e0((int)it[0],(int)it[1]);
        pair<int,int> e1((int)it[1],(int)it[2]);
        pair<int,int> e2((int)it[2],(int)it[0]);
        bool p0=e0.first < e0.second;
        bool p1=e1.first < e1.second;
        bool p2=e2.first < e2.second;
        if(!p0)std::swap(e0.first,e0.second);
        if(!p1)std::swap(e1.first,e1.second);
        if(!p2)std::swap(e2.first,e2.second);
        _tss[i].reset(new ClothTriangle(_ess[eIdMap[e0]],p0,_ess[eIdMap[e1]],p1,_ess[eIdMap[e2]],p2,type));

        for(sizeType e=0; e<3; e++) {
            if(_tss[i]->_e[e]->_t[0])
                _tss[i]->_e[e]->_t[1]=_tss[i];
            else _tss[i]->_e[e]->_t[0]=_tss[i];
        }

        //face one ring
        _tss[i]->getV0()->_oneRing.push_back(_tss[i]);
        _tss[i]->getV1()->_oneRing.push_back(_tss[i]);
        _tss[i]->getV2()->_oneRing.push_back(_tss[i]);
    }

    //face one ring validify
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
        ClothMesh::ClothVertex& v=*(_vss[i]);
        std::sort(v._oneRing.begin(),v._oneRing.end());

        vector<boost::shared_ptr<ClothTriangle> >::iterator iter=
            std::unique(v._oneRing.begin(),v._oneRing.end());
        v._oneRing.erase(iter,v._oneRing.end());
    }
    assembleIndex();
    assembleA();
}
bool ClothMesh::write(std::ostream& os) const
{
    IOData dat;
    dat.registerType<ClothVertex>();
    dat.registerType<ClothEdge>();
    dat.registerType<ClothTriangle>();
    writeVector(_vss,os,&dat);
    writeVector(_ess,os,&dat);
    writeVector(_tss,os,&dat);
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        writeVector(_vss[i]->_oneRing,os,&dat);
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        writeBinaryData(_ess[i]->_t[0],os,&dat);
        writeBinaryData(_ess[i]->_t[1],os,&dat);
    }
    return os.good();
}
bool ClothMesh::read(std::istream& is)
{
    IOData dat;
    dat.registerType<ClothVertex>();
    dat.registerType<ClothEdge>();
    dat.registerType<ClothTriangle>();
    readVector(_vss,is,&dat);
    readVector(_ess,is,&dat);
    readVector(_tss,is,&dat);
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        readVector(_vss[i]->_oneRing,is,&dat);
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        readBinaryData(_ess[i]->_t[0],is,&dat);
        readBinaryData(_ess[i]->_t[1],is,&dat);
    }
    return is.good();
}
void ClothMesh::writeVTKC(const std::string& str,char defColor,const std::vector<char>* color,bool last) const
{
    VTKWriter<scalarD> os("Cloth",str,true);
    std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> > vss(_vss.size());
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss(_tss.size());
    std::vector<scalarD> css(vss.size());
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        vss[_vss[i]->_index]=last ? _vss[i]->_lastPos : _vss[i]->_pos;
        css[_vss[i]->_index]=(color && (sizeType)(color->size()) > i) ? (*color)[i] : defColor;
    }
    for(sizeType i=0; i<(sizeType)_tss.size(); i++)
        iss[i]=Vec3i(_tss[i]->getV0()->_index,
                     _tss[i]->getV1()->_index,
                     _tss[i]->getV2()->_index);
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(iss.begin(),iss.end(),VTKWriter<scalarD>::TRIANGLE);
    os.appendCustomPointData("Color",css.begin(),css.end());
}
void ClothMesh::writeVTKN(const std::string& str) const
{
    VTKWriter<scalarD> os("Cloth",str,true);
    std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> > vss(_tss.size()*3);
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss(_tss.size());
    for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
        vss[i*3+0]=_tss[i]->_e[0]->_pos+
                   _tss[i]->_e[2]->_pos-
                   _tss[i]->_e[1]->_pos;
        vss[i*3+1]=_tss[i]->_e[0]->_pos+
                   _tss[i]->_e[1]->_pos-
                   _tss[i]->_e[2]->_pos;
        vss[i*3+2]=_tss[i]->_e[1]->_pos+
                   _tss[i]->_e[2]->_pos-
                   _tss[i]->_e[0]->_pos;
        iss[i]=Vec3i(i*3+0,i*3+1,i*3+2);
    }
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(iss.begin(),iss.end(),VTKWriter<scalarD>::TRIANGLE);
}
void ClothMesh::assignMass(MASS_MODE mode)
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_mass=0.0f;
    for(sizeType i=0; i<(sizeType)_ess.size(); i++)
        _ess[i]->_mass=0.0f;
    for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
        ClothTriangle& tri=*(_tss[i]);
        TriangleTpl<scalarD> triGeom(tri.getV0()->_pos,tri.getV1()->_pos,tri.getV2()->_pos);

        if(mode == BARYCENTER) {
            scalarD A=triGeom.area();
            tri.getV0()->_mass+=A/3.0f;
            tri.getV1()->_mass+=A/3.0f;
            tri.getV2()->_mass+=A/3.0f;
            tri._e[0]->_mass+=A/3.0f;
            tri._e[1]->_mass+=A/3.0f;
            tri._e[2]->_mass+=A/3.0f;
        } else if(mode == CIRCUMCENTER) {
            Vec3d circum=triGeom.circumcenter();
            Vec3d cBary=triGeom.bary(circum);
            if(cBary[0] < 0.0f || cBary[1] < 0.0f || cBary[2] < 0.0f) {
                INFOV("Triangle %d Not Well-Centered (%f,%f,%f)!",i,cBary[0],cBary[1],cBary[2])
                ASSERT((triGeom._a*cBary[0]+triGeom._b*cBary[1]+triGeom._c*cBary[2]-circum).norm() < EPS)
            }

            scalarD ae0=abs(TriangleTpl<scalarD>(circum,triGeom._a,triGeom._b).area());
            scalarD ae1=abs(TriangleTpl<scalarD>(circum,triGeom._b,triGeom._c).area());
            scalarD ae2=abs(TriangleTpl<scalarD>(circum,triGeom._c,triGeom._a).area());
            tri.getV0()->_mass+=(ae0+ae2)/2.0f;
            tri.getV1()->_mass+=(ae0+ae1)/2.0f;
            tri.getV2()->_mass+=(ae1+ae2)/2.0f;
            tri._e[0]->_mass+=ae0;
            tri._e[1]->_mass+=ae1;
            tri._e[2]->_mass+=ae2;
        }
    }
}
void ClothMesh::saveLast()
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_lastPos=_vss[i]->_pos;
}
void ClothMesh::findBoundary(std::vector<std::vector<boost::shared_ptr<ClothTriangle> > >& boundary)
{
    //initalize
    assembleIndex();
    boundary.clear();
    boundary.resize(_vss.size());

    //find all boundary vertex
    std::vector<bool> bTag(_vss.size(),false);
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        if(!_ess[i]->_t[1]) {
            bTag[_ess[i]->_v[0]->_index]=true;
            bTag[_ess[i]->_v[1]->_index]=true;
        }
    }
    for(sizeType i=0,id; i<(sizeType)_tss.size(); i++) {
        //v0
        id=_tss[i]->getV0()->_index;
        if(bTag[id])boundary[id].push_back(_tss[i]);
        //v1
        id=_tss[i]->getV1()->_index;
        if(bTag[id])boundary[id].push_back(_tss[i]);
        //v2
        id=_tss[i]->getV2()->_index;
        if(bTag[id])boundary[id].push_back(_tss[i]);
    }
}
void ClothMesh::assembleIndex()
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_index=i;
    for(sizeType i=0; i<(sizeType)_ess.size(); i++)
        _ess[i]->_index=i;
    for(sizeType i=0; i<(sizeType)_tss.size(); i++)
        _tss[i]->_index=i;
}
void ClothMesh::assembleA()
{
    //build weight
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_weight=0.0f;
    for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
        _tss[i]->getV0()->_weight+=1.0f;
        _tss[i]->getV1()->_weight+=1.0f;
        _tss[i]->getV2()->_weight+=1.0f;
    }

    //build matrix
    std::vector<Eigen::Triplet<scalarD,sizeType> > trips;
    _A.resize((int)_vss.size()*3,(int)_ess.size()*3);
    for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
        ClothTriangle& tri=*(_tss[i]);
        scalarD deg;
        //contribute V0
        deg=tri.getV0()->_weight;
        addI3x3(trips,tri.getV0()->_index*3,tri._e[0]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV0()->_index*3,tri._e[2]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV0()->_index*3,tri._e[1]->_index*3,-1.0f/deg);
        //contribute V1
        deg=tri.getV1()->_weight;
        addI3x3(trips,tri.getV1()->_index*3,tri._e[0]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV1()->_index*3,tri._e[1]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV1()->_index*3,tri._e[2]->_index*3,-1.0f/deg);
        //contribute V2
        deg=tri.getV2()->_weight;
        addI3x3(trips,tri.getV2()->_index*3,tri._e[1]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV2()->_index*3,tri._e[2]->_index*3,1.0f/deg);
        addI3x3(trips,tri.getV2()->_index*3,tri._e[0]->_index*3,-1.0f/deg);
    }
    _A.setFromTriplets(trips.begin(),trips.end());

    //build solver
    _AAT=_A*_A.transpose();
    _AATSol.compute(_AAT);
}
void ClothMesh::assembleN(Vec* N,Vec* NV,Vec* M,bool pos0) const
{
    if(N)N->resize(_ess.size()*3);
    if(NV)NV->resize(_ess.size()*3);
    if(M)M->resize(_ess.size()*3);
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        if(N)N->block<3,1>(_ess[i]->_index*3,0)=pos0 ? _ess[i]->_pos0 : _ess[i]->_pos;
        if(NV)NV->block<3,1>(_ess[i]->_index*3,0)=_ess[i]->_vel;
        if(M)M->block<3,1>(_ess[i]->_index*3,0).setConstant(_ess[i]->_mass);
    }
}
void ClothMesh::assembleC(Vec* C,Vec* CV,Vec* M,bool pos0) const
{
    if(C)C->resize(_vss.size()*3);
    if(CV)CV->resize(_vss.size()*3);
    if(M)M->resize(_vss.size()*3);
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        if(C)C->block<3,1>(_vss[i]->_index*3,0)=pos0 ? _vss[i]->_pos0 : _vss[i]->_pos;
        if(CV)CV->block<3,1>(_vss[i]->_index*3,0)=_vss[i]->_vel;
        if(M)M->block<3,1>(_vss[i]->_index*3,0).setConstant(_vss[i]->_mass);
    }
}
void ClothMesh::assignN(const Vec* N,const Vec* NV)
{
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        if(N)_ess[i]->_pos=N->block<3,1>(_ess[i]->_index*3,0);
        if(NV)_ess[i]->_vel=NV->block<3,1>(_ess[i]->_index*3,0);
    }
}
void ClothMesh::assignC(const Vec* C,const Vec* CV)
{
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        if(C)_vss[i]->_pos=C->block<3,1>(_vss[i]->_index*3,0);
        if(CV)_vss[i]->_vel=CV->block<3,1>(_vss[i]->_index*3,0);
    }
}
void ClothMesh::convertN2C()
{
    Vec N,C;
    assembleN(&N);
    C=_A*N;
    assignC(&C);
}
void ClothMesh::convertC2N()
{
    Vec C,N,RHS,L;
    assembleN(&N);
    assembleC(&C);
    RHS=_A*N-C;
    L=_AATSol.solve(RHS);
    N-=_A.transpose()*L;
    assignN(&N);
}
void ClothMesh::convertC2Obj(ObjMeshD& mesh)
{
    mesh.getV().resize(_vss.size());
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        mesh.getV()[_vss[i]->_index]=_vss[i]->_pos;

    mesh.getI().resize(_tss.size());
    for(sizeType i=0; i<(sizeType)_tss.size(); i++)
        mesh.getI()[i]=Vec3i(_tss[i]->getV0()->_index,
                             _tss[i]->getV1()->_index,
                             _tss[i]->getV2()->_index);

    mesh.smooth();
}
Vec3d ClothMesh::perp(const Vec3d& in)
{
    if(abs(in[0]) <= abs(in[1]) && abs(in[0]) <= abs(in[2]))
        return Vec3d(0.0f,in[2],-in[1]).normalized();
    else if(abs(in[1]) <= abs(in[0]) && abs(in[1]) <= abs(in[2]))
        return Vec3d(in[2],0.0f,-in[0]).normalized();
    else return Vec3d(in[1],-in[0],0.0f).normalized();
}
void ClothMesh::addI3x3(std::vector<Eigen::Triplet<scalarD,sizeType> >& trips,sizeType R,sizeType C,const scalarD& coef)
{
    trips.push_back(Eigen::Triplet<scalarD,sizeType>(R+0,C+0,coef));
    trips.push_back(Eigen::Triplet<scalarD,sizeType>(R+1,C+1,coef));
    trips.push_back(Eigen::Triplet<scalarD,sizeType>(R+2,C+2,coef));
}
void ClothMesh::parityCheck()
{
    for(sizeType i=0; i<(sizeType)_ess.size(); i++) {
        for(sizeType j=0; j<2; j++)
            if(_ess[i]->_t[j]) {
                ASSERT (_ess[i]->_t[j]->_e[0] == _ess[i] ||
                        _ess[i]->_t[j]->_e[1] == _ess[i] ||
                        _ess[i]->_t[j]->_e[2] == _ess[i]);
            }
    }
    for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
        for(sizeType j=0; j<3; j++) {
            boost::shared_ptr<ClothVertex> va=
                _tss[i]->_edgeDir[j] ?
                _tss[i]->_e[j]->_v[0] :
                _tss[i]->_e[j]->_v[1];

            sizeType lastJ=(j+2)%3;
            boost::shared_ptr<ClothVertex> vb=
                _tss[i]->_edgeDir[lastJ] ?
                _tss[i]->_e[lastJ]->_v[1] :
                _tss[i]->_e[lastJ]->_v[0];

            ASSERT(va == vb);
        }
    }

    {
        assignMass(BARYCENTER);
        scalarD mv=0.0f,me=0.0f,mt=0.0f;
        for(sizeType i=0; i<(sizeType)_vss.size(); i++)
            mv+=_vss[i]->_mass;
        for(sizeType i=0; i<(sizeType)_ess.size(); i++)
            me+=_ess[i]->_mass;
        for(sizeType i=0; i<(sizeType)_tss.size(); i++)
            mt+=TriangleTpl<scalarD>
                (_tss[i]->getV0()->_pos,
                 _tss[i]->getV1()->_pos,
                 _tss[i]->getV2()->_pos).area();
        ASSERT(abs(mv-me) < EPS && abs(mv-mt) < EPS)
    }

    {
        assignMass(CIRCUMCENTER);
        scalarD mv=0.0f,me=0.0f;
        for(sizeType i=0; i<(sizeType)_vss.size(); i++)
            mv+=_vss[i]->_mass;
        for(sizeType i=0; i<(sizeType)_ess.size(); i++)
            me+=_ess[i]->_mass;
        ASSERT(abs(mv-me) < EPS)
    }
}
void ClothMesh::parityCheckVss()
{
    //assemble N
    Vec N,C;
    N.resize(_ess.size()*3);
    for(sizeType i=0; i<(sizeType)_ess.size(); i++)
        N.block<3,1>(_ess[i]->_index*3,0)=_ess[i]->_pos;

    //assemble C
    C=_A*N;
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        ASSERT((_vss[i]->_pos-C.block<3,1>(_vss[i]->_index*3,0)).norm() < EPS)
    }
Mat4d ClothMesh::getTransform() const
{
    Mat4d LHS=Mat4d::Zero(),RHS=Mat4d::Zero();
    for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
        Vec4d pos0;
        pos0.block<3,1>(0,0)=_vss[i]->_pos0;
        pos0[3]=1.0f;

        Vec4d pos1;
        pos1.block<3,1>(0,0)=_vss[i]->_pos;
        pos1[3]=1.0f;

        LHS+=pos0*pos0.transpose();
        RHS+=pos0*pos1.transpose();
    }
    return (LHS.ldlt().solve(RHS)).eval().transpose();
}
void ClothMesh::transform(const Mat3d& R,const Vec3d& T)
{
    saveLast();
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
        _vss[i]->_pos=R*_vss[i]->_pos+T;
}
