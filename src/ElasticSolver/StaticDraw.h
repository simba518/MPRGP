#ifndef STATIC_DRAW_H
#define STATIC_DRAW_H

#include "MathBasic.h"
#include "ObjMesh.h"
#include "PVSMScene.h"

PRJ_BEGIN

class StaticDraw : public PVSMScene::PVSMAction
{
public:
    StaticDraw(const std::string& group);
    ~StaticDraw();
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool postDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    void setDebugNormal(bool debug,scalar len);
    static void toArray(ObjMesh& mesh,boost::shared_ptr<PVSMScene::PVSMSource> src,Mat4 T0,
                        Eigen::Matrix<GLfloat,-1,1>& vertices,Eigen::Matrix<GLfloat,-1,1>& normals,
                        Eigen::Matrix<GLuint,-1,1>& indices,GLenum& type);
protected:
    void toArray(ObjMesh& mesh,boost::shared_ptr<PVSMScene::PVSMSource> src,Mat4 T0=Mat4::Identity());
    static Mat3 getNMat(const Mat4& T);
    Eigen::Matrix<GLfloat,-1,1> _vertices;
    Eigen::Matrix<GLfloat,-1,1> _normals;
    Eigen::Matrix<GLuint,-1,1> _indices;
    GLenum _type;
    GLuint _list;
    bool _debugNormal,_displayList;
    scalar _normalLen;
};
class StaticObjMesh : public StaticDraw
{
public:
    StaticObjMesh();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticObjMesh);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class StaticBox : public StaticDraw
{
public:
    StaticBox();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticBox);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class StaticPlane : public StaticDraw
{
public:
    StaticPlane();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticPlane);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class StaticSphere : public StaticDraw
{
public:
    StaticSphere();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticSphere);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class StaticLine : public StaticDraw
{
public:
    StaticLine();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticLine);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class StaticPoint : public StaticDraw
{
public:
    StaticPoint();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new StaticPoint);
    }
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
};
class Light : public PVSMScene::PVSMAction
{
public:
    struct LightManager {
        LightManager() {
            _id=0;
        }
        GLuint grabId() {
            return _id++;
        }
        GLuint _id;
    };
    Light(boost::shared_ptr<LightManager> mgr=boost::shared_ptr<LightManager>(new LightManager));
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return boost::shared_ptr<PVSMAction>(new Light(_mgr));
    }
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
private:
    Vec4 _pos,_dir;
    boost::shared_ptr<LightManager> _mgr;
};

PRJ_END

#endif