#ifndef FEM_EMBEDDED_MESH_H
#define FEM_EMBEDDED_MESH_H

#include "FEMMesh.h"
#include "StaticDraw.h"

PRJ_BEGIN

class FEMEmbeddedMesh : public StaticDraw
{
public:
    FEMEmbeddedMesh(const FEMMesh& mesh);
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self);
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
protected:
    const FEMMesh& _mesh;
    vector<FEMInterp> _evss;
    Eigen::Matrix<GLfloat,-1,1> _normals0;
};
class FEMGPUEmbeddedMesh : public PVSMScene::PVSMAction
{
public:
    FEMGPUEmbeddedMesh(const FEMMesh& mesh);
    ~FEMGPUEmbeddedMesh();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self);
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
private:
    const FEMMesh& _mesh;
    boost::shared_ptr<FEMSystem> _sys;
    sizeType _szSys;
    //texture buffer object
    GLuint _ctbuf,_ctid;
    GLuint _btbuf,_btid;
    //transform feedback
    cwc::glShader *_shaderRecon;
    GLuint _vbuf0,_vbuf1;
    GLuint _nbuf0,_nbuf1;
    GLuint _ibuf;
};
class FEMComputeMesh : public StaticDraw
{
public:
    FEMComputeMesh(const FEMMesh& mesh);
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return self;
    }
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
private:
    const FEMMesh& _mesh;
    vector<boost::shared_ptr<FEMVertex> > _verts;
};

PRJ_END

#endif