#include "StaticDraw.h"
#include "MakeMesh.h"
#include <GL/freeglut.h>

USE_PRJ_NAMESPACE

StaticDraw::StaticDraw(const std::string& group):PVSMScene::PVSMAction(group)
{
    setDebugNormal(false,0.002f);
    _displayList=false;
    _list=numeric_limits<GLuint>::max();
}
StaticDraw::~StaticDraw()
{
    if(glIsList(_list))
        glDeleteLists((GLuint)_list,1);
}
bool StaticDraw::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    if(_vertices.size() == 0)
        return false;
    if(_displayList && _list == -1) {
        _list=glGenLists(1);
        glNewList((GLuint)_list,GL_COMPILE);
        _displayList=false;
        draw(s,src);
        _displayList=true;
        glEndList();
    }
    return false;
}
bool StaticDraw::draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    if(_vertices.size() == 0)
        return false;
    if(_displayList) {
        src->_mat._ref->setOpenGLMaterial();
        glCallList((GLuint)_list);
    } else {
        src->_mat._ref->setOpenGLMaterial();
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3,GL_FLOAT,0,_vertices.data());
        if(_normals.size() > 0) {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT,0,_normals.data());
        }
        glDrawElements(_type,(GLsizei)_indices.size(),GL_UNSIGNED_INT,_indices.data());
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    return false;
}
bool StaticDraw::postDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    if(_debugNormal && _normals.size() == _vertices.size()) {
        glPushAttrib(GL_ENABLE_BIT|GL_LIGHTING_BIT|GL_LINE_BIT);
        glDisable(GL_LIGHTING);
        glLineWidth(1.0f);
        glBegin(GL_LINES);
        for(sizeType i=0; i<_normals.size(); i+=3) {
            glColor3f(1.0f,1.0f,1.0f);
            glVertex3d(_vertices[i+0],
                       _vertices[i+1],
                       _vertices[i+2]);
            glVertex3d(_vertices[i+0]+_normals[i+0]*_normalLen,
                       _vertices[i+1]+_normals[i+1]*_normalLen,
                       _vertices[i+2]+_normals[i+2]*_normalLen);
        }
        glEnd();
        glEnable(GL_LIGHTING);
        glPopAttrib();
    }
    return false;
}
void StaticDraw::setDebugNormal(bool debug,scalar len)
{
    _debugNormal=debug;
    _normalLen=len;
}
void StaticDraw::toArray
(ObjMesh& mesh,boost::shared_ptr<PVSMScene::PVSMSource> src,Mat4 T0,
 Eigen::Matrix<GLfloat,-1,1>& vertices,Eigen::Matrix<GLfloat,-1,1>& normals,
 Eigen::Matrix<GLuint,-1,1>& indices,GLenum& type)
{
    Mat4 T=PVSMScene::getTrans(src)*T0;
    vertices.setOnes(mesh.getV().size()*3);
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
        mesh.getV()[i]=PVSMScene::transform(T,mesh.getV()[i]);
        vertices.block<3,1>(i*3,0)=mesh.getV()[i].cast<GLfloat>();
    }

    Mat3 N=getNMat(T);
    normals.setOnes(mesh.getN().size()*3);
    for(sizeType i=0; i<(sizeType)mesh.getN().size(); i++)
        normals.block<3,1>(i*3,0)=(N*mesh.getN()[i]).normalized().cast<GLfloat>();

    indices.resize(mesh.getI().size()*3);
    for(sizeType i=0; i<(sizeType)mesh.getI().size(); i++)
        indices.block<3,1>(i*3,0)=mesh.getI()[i].cast<GLuint>();

    type=GL_TRIANGLES;
}
void StaticDraw::toArray(ObjMesh& mesh,boost::shared_ptr<PVSMScene::PVSMSource> src,Mat4 T0)
{
    toArray(mesh,src,T0,_vertices,_normals,_indices,_type);
}
Mat3 StaticDraw::getNMat(const Mat4& T)
{
    Mat3 ret=T.block<3,3>(0,0).inverse();
    return ret.transpose();
}
//StaticObjMesh
StaticObjMesh::StaticObjMesh():StaticDraw("staticObjMesh")
{
    _names["staticObjMesh"]=PVSMScene::OBJMESH;
    _displayList=true;
}
bool StaticObjMesh::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMObjMesh> msh=
        boost::dynamic_pointer_cast<PVSMScene::PVSMObjMesh>(src);
    if(msh) {
        ObjMesh mesh;
        boost::filesystem::ifstream is(msh->_file[0]);
        mesh.read(is,false,false);
        toArray(mesh,src);
    }
    return false;
}
//StaticBox
StaticBox::StaticBox():StaticDraw("staticBox")
{
    _names["staticBox"]=PVSMScene::BOX;
    _displayList=true;
}
bool StaticBox::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMBox> box=
        boost::dynamic_pointer_cast<PVSMScene::PVSMBox>(src);
    if(box) {
        ObjMesh mesh;
        MakeMesh::makeDiscreteBox3D(mesh,Vec3(box->_X[0],box->_Y[0],box->_Z[0])/2.0f);
        Mat4 T0=Mat4::Identity();
        T0.block<3,1>(0,3)=box->_Ctr.get().cast<scalar>();
        toArray(mesh,src,T0);
    }
    return false;
}
//StaticPlane
StaticPlane::StaticPlane():StaticDraw("staticPlane")
{
    _names["staticPlane"]=PVSMScene::PLANE;
    _displayList=true;
}
bool StaticPlane::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMPlane> plane=
        boost::dynamic_pointer_cast<PVSMScene::PVSMPlane>(src);
    if(plane) {
        ObjMesh mesh;
        mesh.getV().push_back(plane->_Ctr.get());
        mesh.getV().push_back(plane->_P1.get());
        mesh.getV().push_back(plane->_P2.get());
        mesh.getV().push_back(plane->_P1.get()+plane->_P2.get()-plane->_Ctr.get());
        mesh.getI().push_back(Vec3i(0,1,3));
        mesh.getI().push_back(Vec3i(0,3,2));
        mesh.smooth();
        toArray(mesh,src);
    }
    return false;
}
//StaticSphere
StaticSphere::StaticSphere():StaticDraw("staticSphere")
{
    _names["staticSphere"]=PVSMScene::SPHERE;
    _displayList=true;
}
bool StaticSphere::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMSphere> sphere=
        boost::dynamic_pointer_cast<PVSMScene::PVSMSphere>(src);
    if(sphere) {
        ObjMesh mesh;
        MakeMesh::makeSphere3D(mesh,sphere->_Rad[0],32);
        Mat4 T0=Mat4::Identity();
        T0.block<3,1>(0,3)=sphere->_Ctr.get().cast<scalar>();
        toArray(mesh,src,T0);
    }
    return false;
}
//StaticLine
StaticLine::StaticLine():StaticDraw("staticLine")
{
    _names["staticLine"]=PVSMScene::LINE;
    _displayList=true;
}
bool StaticLine::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Mat4 T=PVSMScene::getTrans(src);
    boost::shared_ptr<PVSMScene::PVSMLine> line=
        boost::dynamic_pointer_cast<PVSMScene::PVSMLine>(src);
    if(line) {
        _vertices.resize(6);
        _vertices.block<3,1>(0,0)=PVSMScene::transform(T,line->_P1.get()).cast<GLfloat>();
        _vertices.block<3,1>(3,0)=PVSMScene::transform(T,line->_P2.get()).cast<GLfloat>();

        _indices.resize(2);
        _indices=Vec2i(0,1).cast<GLuint>();
        _type=GL_LINES;
    }
    return false;
}
//StaticPoint
StaticPoint::StaticPoint():StaticDraw("staticPoint")
{
    _names["staticPoint"]=PVSMScene::POINT;
    _displayList=true;
}
bool StaticPoint::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Mat4 T=PVSMScene::getTrans(src);
    boost::shared_ptr<PVSMScene::PVSMPoint> point=
        boost::dynamic_pointer_cast<PVSMScene::PVSMPoint>(src);
    if(point) {
        _vertices.resize(3);
        _vertices=PVSMScene::transform(T,point->_Ctr.get()).cast<GLfloat>();

        _indices.resize(1);
        _indices[0]=0;
        _type=GL_POINTS;
    }
    return false;
}
//Light
Light::Light(boost::shared_ptr<LightManager> mgr):PVSMScene::PVSMAction("light"),_mgr(mgr)
{
    _names["lightDir"]=PVSMScene::LINE;
}
bool Light::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    src->_mat._ref->setOpenGLLight(s);
    return false;
}
bool Light::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Mat4 T=PVSMScene::getTrans(src);
    boost::shared_ptr<PVSMScene::PVSMLine> line=
        boost::dynamic_pointer_cast<PVSMScene::PVSMLine>(src);
    if(line) {
        _pos=PVSMScene::transformH(T,line->_P1.get());
        _dir=PVSMScene::transformH(T,line->_P2.get())-_pos;
    }

    src->_mat._ref->_lightId=_mgr->grabId();
    src->_mat._ref->_LPos.set<GLfloat>(_pos.cast<GLfloat>());
    src->_mat._ref->_LDir.set<GLfloat>(_dir.cast<GLfloat>());
    return false;
}
