#include "geom/staticGeom.h"
#include "geom/BVHBuilder.h"
#include "MakeMesh.h"
#include "CollisionDetection.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

struct BoxGeomCell : public StaticGeomCell {
    BoxGeomCell():StaticGeomCell(0) {}
    BoxGeomCell(const Mat4& T,sizeType dim,const Vec3& ext):StaticGeomCell(T,dim,0),_ext(ext) {}
    virtual void getMesh(ObjMeshTpl<scalar>& mesh) const {
        if(_dim == 2)MakeMesh::makeBox2D(mesh,_ext);
        else MakeMesh::makeBox3D(mesh,_ext);
        StaticGeomCell::getMesh(mesh);
    }
    virtual bool read(std::istream& is) {
        StaticGeomCell::read(is);
        readBinaryData(_ext,is);
        return is.good();
    }
    virtual bool write(std::ostream& os) const {
        StaticGeomCell::write(os);
        writeBinaryData(_ext,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new BoxGeomCell);
    }
protected:
    virtual BBox<scalar> getBBInner() const {
        return BBox<scalar>(-_ext,_ext);
    }
    virtual bool distInner(const Vec3& pt,Vec3& n) const {
        scalar minDist=numeric_limits<scalar>::max(),dist;
        for(sizeType d=0; d<_dim; d++) {
            dist=_ext[d]+pt[d];
            if(dist < minDist) {
                minDist=dist;
                n=-Vec3::Unit(d)*dist;
            }
            dist=_ext[d]-pt[d];
            if(dist < minDist) {
                minDist=dist;
                n=Vec3::Unit(d)*dist;
            }
        }
        return minDist > 0.0f;
    }
    Vec3 _ext;
};
struct SphereGeomCell : public StaticGeomCell {
    SphereGeomCell():StaticGeomCell(1) {}
    SphereGeomCell(const Mat4& T,sizeType dim,scalar rad):StaticGeomCell(T,dim,1),_rad(rad) {}
    virtual void getMesh(ObjMeshTpl<scalar>& mesh) const {
        if(_dim == 2)MakeMesh::makeSphere2D(mesh,_rad,16);
        else MakeMesh::makeSphere3D(mesh,_rad,16);
        StaticGeomCell::getMesh(mesh);
    }
    virtual bool read(std::istream& is) {
        StaticGeomCell::read(is);
        readBinaryData(_rad,is);
        return is.good();
    }
    virtual bool write(std::ostream& os) const {
        StaticGeomCell::write(os);
        writeBinaryData(_rad,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new SphereGeomCell);
    }
protected:
    virtual BBox<scalar> getBBInner() const {
        Vec3 rad=Vec3::Zero();
        rad.block(0,0,_dim,1).setConstant(_rad);
        return BBox<scalar>(-rad,rad);
    }
    virtual bool distInner(const Vec3& pt,Vec3& n) const {
        scalar len=pt.norm();
        if(len > _rad)
            return false;
        n=pt*(_rad-len)/std::max<scalar>(len,1E-6f);
        return true;
    }
    scalar _rad;
};
struct CylinderGeomCell : public StaticGeomCell {
    CylinderGeomCell():StaticGeomCell(2) {}
    CylinderGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y):StaticGeomCell(T,dim,2),_rad(rad),_y(y) {}
    virtual void getMesh(ObjMeshTpl<scalar>& mesh) const {
        MakeMesh::makeCylinder3D(mesh,_rad,_y,64,64);
        StaticGeomCell::getMesh(mesh);
    }
    virtual bool read(std::istream& is) {
        StaticGeomCell::read(is);
        readBinaryData(_rad,is);
        readBinaryData(_y,is);
        return is.good();
    }
    virtual bool write(std::ostream& os) const {
        StaticGeomCell::write(os);
        writeBinaryData(_rad,os);
        writeBinaryData(_y,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new CylinderGeomCell);
    }
protected:
    virtual BBox<scalar> getBBInner() const {
        Vec3 cor(_rad,_y,_rad);
        if(_dim == 2)cor[2]=0.0f;
        return BBox<scalar>(-cor,cor);
    }
    virtual bool distInner(const Vec3& pt,Vec3& n) const {
        scalar len=Vec3(pt[0],0.0f,pt[2]).norm();
        //boundary
        scalar dist=_rad-len;
        if(dist < 0.0f)return false;
        n=Vec3(pt[0],0.0f,pt[2])*dist/std::max<scalar>(len,1E-6f);
        //bottom
        scalar dist2=_y+pt[1];
        if(dist2 < 0.0f)return false;
        if(dist2 < dist) {
            dist=dist2;
            n=-Vec3::Unit(1)*dist2;
        }
        //top
        scalar dist3=_y-pt[1];
        if(dist3 < 0.0f)return false;
        if(dist3 < dist) {
            dist=dist3;
            n=Vec3::Unit(1)*dist3;
        }
        return true;
    }
    scalar _rad,_y;
};
struct ObjMeshGeomCell : public StaticGeomCell {
    ObjMeshGeomCell():StaticGeomCell(3),_query(_bvh,_dim,-1),_insideOut(false) {}
    ObjMeshGeomCell(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth,bool insideOut)
        :StaticGeomCell(trans,mesh.getDim(),3),_depth(depth),_query(_bvh,_dim,-1),_insideOut(insideOut) {
        _vss=mesh.getV();
        _iss=mesh.getI();
        _bvh.resize(_iss.size());
        for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
            Node<sizeType>& n=_bvh[i];
            n._nrCell=1;
            n._cell=i;
            n._bb.reset();
            for(sizeType v=0; v<_dim; v++)
                n._bb.setUnion(_vss[_iss[i][v]]);
        }
        buildBVH<sizeType>(_bvh,_dim,-1);
        // writeBVHByLevel<sizeType>(_bvh,-1);
        if(_depth == 0.0f)
            _depth=mesh.getBB().getExtent().maxCoeff()*0.1f;
    }
    void getMesh(ObjMeshTpl<scalar>& mesh) const {
        mesh.getV()=_vss;
        mesh.getI()=_iss;
        mesh.setDim((int)_dim);
        mesh.smooth();
        StaticGeomCell::getMesh(mesh);
    }
    bool read(std::istream& is) {
        StaticGeomCell::read(is);
        readVector(_vss,is);
        readVector(_iss,is);
        readVector(_bvh,is);
        readBinaryData(_depth,is);
        readBinaryData(_insideOut,is);
        return is.good();
    }
    bool write(std::ostream& os) const {
        StaticGeomCell::write(os);
        writeVector(_vss,os);
        writeVector(_iss,os);
        writeVector(_bvh,os);
        writeBinaryData(_depth,os);
        writeBinaryData(_insideOut,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new ObjMeshGeomCell);
    }
    scalar depth() const {
        return _depth;
    }
    void updateDist(const Node<sizeType>& node,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist) const {
        if(_dim == 2)calcDist2D(_iss[node._cell],pt,cp,n,dist);
        else calcDist3D(_iss[node._cell],pt,cp,n,dist);
    }
protected:
    virtual BBox<scalar> getBBInner() const {
        BBox<scalar> bb=_bvh.back()._bb;
        if(_insideOut)
            return bb.enlarge(_depth,_dim);
        else return bb;
    }
    virtual bool distInner(const Vec3& pt,Vec3& n) const {
        Vec3 cp;
        scalar dist=numeric_limits<scalar>::max();
        _query.pointDistQuery(pt,*this,cp,n,dist);
        cp-=pt;
        if(_insideOut)
            n*=-1.0f;
        if(cp.dot(n) < 0.0f)
            return false;
        n=cp;
        return dist < numeric_limits<scalar>::max();
    }
	virtual void closestInner(const Vec3& pt,Vec3& n) const
	{
		Vec3 cp;
        scalar dist=numeric_limits<scalar>::max();
        _query.pointDistQuery(pt,*this,cp,n,dist);
        n=cp-pt;
		ASSERT(dist < numeric_limits<scalar>::max());
	}
    void calcDist2D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist) const {
        Vec3 cpTmp,b;
        scalar distNew;
        LineSegTpl<scalar> l(_vss[I[0]],_vss[I[1]]);

        l.calcPointDist(pt,distNew,cpTmp,b);
        distNew=sqrt(distNew);
        if(distNew < dist) {
            dist=distNew;
            cp=cpTmp;
            n=l.normal();
        }
    }
    void calcDist3D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist) const {
        Vec3 cpTmp,b;
        scalar distNew;
        Triangle t(_vss[I[0]],_vss[I[1]],_vss[I[2]]);

        t.calcPointDist(pt,distNew,cpTmp,b);
        distNew=sqrt(distNew);
        if(distNew < dist) {
            dist=distNew;
            cp=cpTmp;
            n=t.normal();
        }
    }
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > _vss;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _iss;
    std::vector<Node<sizeType> > _bvh;
    BVHQuery<sizeType> _query;
    scalar _depth;
    char _insideOut;
};

//Geom
StaticGeomCell::StaticGeomCell(sizeType type):Serializable(type) {}
StaticGeomCell::StaticGeomCell(const Mat4& T,sizeType dim,sizeType type):Serializable(type),_T(T),_invT(T.inverse()),_dim(dim) {}
void StaticGeomCell::getMesh(ObjMeshTpl<scalar>& mesh) const
{
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
        Vec3& pos=mesh.getV()[i];
        pos=_T.block<3,3>(0,0)*pos+_T.block<3,1>(0,3);
    }
    mesh.smooth();
}
BBox<scalar> StaticGeomCell::getBB() const
{
    Vec3 pt;
    BBox<scalar> tmp=getBBInner(),ret;
    for(sizeType x=0; x<2; x++)
        for(sizeType y=0; y<2; y++)
            for(sizeType z=0; z<2; z++) {
                pt[0]=(x==0) ? tmp._minC[0] : tmp._maxC[0];
                pt[1]=(y==0) ? tmp._minC[1] : tmp._maxC[1];
                pt[2]=(z==0) ? tmp._minC[2] : tmp._maxC[2];
                ret.setUnion(_T.block<3,3>(0,0)*pt+_T.block<3,1>(0,3));
            }
    return ret;
}
bool StaticGeomCell::dist(const Vec3& pt,Vec3& n) const
{
    Vec3 pt0=_invT.block<3,3>(0,0)*pt+_invT.block<3,1>(0,3);
	if(distInner(pt0,n)){
		n=_T.block<3,3>(0,0)*n;
		return true;
	}
	return false;
}
void StaticGeomCell::closest(const Vec3& pt,Vec3& n) const
{
    Vec3 pt0=_invT.block<3,3>(0,0)*pt+_invT.block<3,1>(0,3);
	closestInner(pt0,n);
	n=_T.block<3,3>(0,0)*n;
}
bool StaticGeomCell::read(std::istream& is)
{
    readBinaryData(_T,is);
    readBinaryData(_invT,is);
    readBinaryData(_dim,is);
    readBinaryData(_index,is);
    return is.good();
}
bool StaticGeomCell::write(std::ostream& os) const
{
    writeBinaryData(_T,os);
    writeBinaryData(_invT,os);
    writeBinaryData(_dim,os);
    writeBinaryData(_index,os);
    return os.good();
}
void StaticGeomCell::closestInner(const Vec3& pt,Vec3& n) const
{
	ASSERT_MSG(false,"Not supported: function closestInner!")
}

//Static Geometry
StaticGeom::StaticGeom(sizeType dim):HasMagic(0xaabbccdd),_dim(dim)
{
    _bvh.reset(new vector<Node<boost::shared_ptr<StaticGeomCell> > >);
}
const vector<Node<boost::shared_ptr<StaticGeomCell> > >& StaticGeom::getBVH() const
{
    return *_bvh;
}
void StaticGeom::assemble()
{
    _bvh->clear();
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        _css[i]->_index=i;
        Node<boost::shared_ptr<StaticGeomCell> > n;
        n._l=n._r=n._parent=-1;
        n._nrCell=1;
        n._cell=_css[i];
        n._bb=n._cell->getBB();
        _bvh->push_back(n);
    }
    buildBVH(*_bvh,_dim,boost::shared_ptr<StaticGeomCell>());
}
void StaticGeom::addGeomCell(boost::shared_ptr<StaticGeomCell> c)
{
    _css.push_back(c);
}
void StaticGeom::addGeomBox(const Mat4& trans,const BBox<scalar>& bb)
{
    Mat4 T=Mat4::Identity();
    T.block<3,1>(0,3)=(bb._maxC+bb._minC)/2.0f;
    addGeomBox(trans*T,bb.getExtent()*0.5f);
}
void StaticGeom::addGeomBox(const Mat4& trans,const Vec3& ext)
{
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new BoxGeomCell(trans,_dim,ext)));
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,2>& obb)
{
    ASSERT(_dim == 2)
    Mat4 m=Mat4::Identity();
    m.block<2,2>(0,0)=obb._rot;
    m.block<2,1>(0,3)=obb._trans;
    addGeomBox(m,Vec3(obb._ext[0],obb._ext[1],0.0f));
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,3>& obb)
{
    ASSERT(_dim == 3)
    Mat4 m=Mat4::Identity();
    m.block<3,3>(0,0)=obb._rot;
    m.block<3,1>(0,3)=obb._trans;
    addGeomBox(m,obb._ext);
}
void StaticGeom::addGeomCylinder(const Mat4& trans,scalar rad,scalar y)
{
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans,_dim,rad,y)));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext)
{
    scalar alpha=-plane[3]/plane.block<3,1>(0,0).squaredNorm();
    Vec3 p0=plane.block<3,1>(0,0)*alpha;
    Quat q;
    q.setFromTwoVectors(Vec3::Unit(1),plane.block<3,1>(0,0).normalized());

    Mat4 T=Mat4::Identity();
    T.block<3,1>(0,3)=p0-plane.block<3,1>(0,0).normalized()*ext;
    T.block<3,3>(0,0)=q.toRotationMatrix();
    if(_dim == 3)
        _css.push_back(boost::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans*T,_dim,ext,ext)));
    else _css.push_back(boost::shared_ptr<StaticGeomCell>(new BoxGeomCell(trans*T,_dim,Vec3::Constant(ext))));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext)
{
    Vec4 p;
    p.block<3,1>(0,0)=plane._n;
    p[3]=-plane._n.dot(plane._x0);
    addGeomPlane(trans,p,ext);
}
void StaticGeom::addGeomSphere(const Vec3& ctr,scalar rad)
{
    Mat4 T=Mat4::Identity();
    T.block<3,1>(0,3)=ctr;
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new SphereGeomCell(T,_dim,rad)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const ObjMeshTpl<scalar>& mesh,scalar depth,bool insideOut)
{
    _css.push_back(boost::shared_ptr<StaticGeomCell>(new ObjMeshGeomCell(trans,mesh,depth,insideOut)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const std::string& path,scalar depth,bool insideOut)
{
    ObjMeshTpl<scalar> mesh;
    boost::filesystem::ifstream is(path);
    mesh.read(is,false,false);
    mesh.smooth();
    mesh.makeUniform();
    addGeomMesh(trans,mesh,depth,insideOut);
}
//IO
void StaticGeom::writeVTK(const std::string& path) const
{
    ObjMeshTpl<scalar> mesh;
    VTKWriter<scalar> os("Geom",path,true);
    boost::filesystem::path components=boost::filesystem::path(path).parent_path()/"geomComponents/";
    boost::filesystem::create_directory(components);
    for(sizeType i=0; i<(sizeType)_css.size(); i++) {
        _css[i]->getMesh(mesh);
        ostringstream oss;
        oss << components.string() << "/comp" << i << ".obj";
        mesh.write(oss.str());
        mesh.writePov(boost::filesystem::ofstream(boost::filesystem::path(oss.str()).replace_extension(".pov")),false);
        mesh.writePov(boost::filesystem::ofstream(boost::filesystem::path(oss.str()).replace_extension(".spov")),true);
        mesh.writeVTK(os,false,false);
    }
}
void StaticGeom::writeBVH() const
{
    writeBVHByLevel<boost::shared_ptr<StaticGeomCell> >(*_bvh,boost::shared_ptr<StaticGeomCell>());
}
bool StaticGeom::write(std::ostream& os) const
{
    if(!HasMagic::writeMagic(os))
        return false;
    IOData dat;
    dat.registerType<BoxGeomCell>();
    dat.registerType<SphereGeomCell>();
    dat.registerType<CylinderGeomCell>();
    dat.registerType<ObjMeshGeomCell>();
    writeBinaryData(_dim,os);
    writeVector(_css,os,&dat);
    writeVector(*_bvh,os,&dat);
    return os.good();
}
bool StaticGeom::read(std::istream& is)
{
    _css.clear();
    _bvh->clear();
    if(!HasMagic::readMagic(is))
        return false;

    IOData dat;
    dat.registerType<BoxGeomCell>();
    dat.registerType<SphereGeomCell>();
    dat.registerType<CylinderGeomCell>();
    dat.registerType<ObjMeshGeomCell>();
    readBinaryData(_dim,is);
    readVector(_css,is,&dat);
    readVector(*_bvh,is,&dat);
    return is.good();
}