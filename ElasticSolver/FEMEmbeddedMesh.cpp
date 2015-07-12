#include "glsl.h"
#include "FEMEmbeddedMesh.h"
#include "FEMSystem.h"
#include "FEMUtils.h"
#include <boost/unordered_set.hpp>

USE_PRJ_NAMESPACE

FEMEmbeddedMesh::FEMEmbeddedMesh(const FEMMesh& mesh):StaticDraw("embeddedMesh"),_mesh(mesh)
{
    _names["embeddedMesh"]=PVSMScene::OBJMESH;
}
boost::shared_ptr<PVSMScene::PVSMAction> FEMEmbeddedMesh::copy(boost::shared_ptr<PVSMAction> self)
{
    return boost::shared_ptr<PVSMAction>(new FEMEmbeddedMesh(_mesh));
}
bool FEMEmbeddedMesh::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<_mesh.nrB(); i++) {
        FEMSystem& sys=*(_mesh.getB(i)._system);
        Cold L;
        sys.getPosL(L);
        sys.body().setDPos(sys.LtoF(L));
    }
    sizeType nrV=(sizeType)_vertices.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i+=3) {
        Vec3 N=_normals0.block<3,1>(i,0).cast<scalar>();
        _vertices.block<3,1>(i,0)=_evss[i/3]._cell->getVert(_evss[i/3],&N).cast<GLfloat>();
        _normals.block<3,1>(i,0)=N.cast<GLfloat>();
    }
    return false;
}
bool FEMEmbeddedMesh::initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMObjMesh> msh=
        boost::dynamic_pointer_cast<PVSMScene::PVSMObjMesh>(src);
    if(msh) {
        ObjMesh mesh;
        boost::filesystem::ifstream is(msh->_file[0]);
        mesh.read(is,false,false);
        toArray(mesh,src);

        _normals0=_normals;
        _evss.resize(mesh.getV().size());
        _mesh.getBary(mesh.getV(),_evss,numeric_limits<scalar>::max());
    }
    return false;
}

template<typename T>
void genBuffer(GLenum type,GLenum mode,GLuint& bid,const T& vec)
{
    glGenBuffersARB(1,&bid);
    glBindBufferARB(type,bid);
    glBufferDataARB(type,vec.size()*sizeof(typename T::Scalar),vec.data(),mode);
    glBindBufferARB(type,0);
}
template<typename T>
void genTBuffer(GLenum type,GLuint& tbo,GLuint& tex,const T& vec)
{
    glGenBuffers(1,&tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT,tbo);
    glBufferData(GL_TEXTURE_BUFFER_EXT,vec.size()*sizeof(typename T::Scalar),vec.data(),GL_STATIC_DRAW);

    glGenTextures(1,&tex);
    glBindTexture(GL_TEXTURE_BUFFER_EXT,tex);
    glTexBufferEXT(GL_TEXTURE_BUFFER_EXT,type,tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT,0);
}
FEMGPUEmbeddedMesh::FEMGPUEmbeddedMesh(const FEMMesh& mesh):_mesh(mesh)
{
    _names["GPUEmbeddedMesh"]=PVSMScene::OBJMESH;
    _ctbuf=_ctid=
               _btbuf=_btid=
                          _vbuf0=_vbuf1=
                                     _nbuf0=_nbuf1=
                                             _ibuf=numeric_limits<GLuint>::max();
}
FEMGPUEmbeddedMesh::~FEMGPUEmbeddedMesh()
{
#define DEL_TEX(TEX) if(glIsTexture(TEX))glDeleteTextures(1,&TEX);
#define DEL_BUF(BUF) if(glIsBufferARB(BUF))glDeleteBuffersARB(1,&BUF);
    DEL_BUF(_ctbuf)
    DEL_TEX(_ctid)
    DEL_BUF(_btbuf)
    DEL_TEX(_btid)
    DEL_BUF(_vbuf0)
    DEL_BUF(_vbuf1)
    DEL_BUF(_nbuf0)
    DEL_BUF(_nbuf1)
    DEL_BUF(_ibuf)
#undef DEL_TEX
#undef DEL_BUF
}
boost::shared_ptr<PVSMScene::PVSMAction> FEMGPUEmbeddedMesh::copy(boost::shared_ptr<PVSMAction> self)
{
    return boost::shared_ptr<PVSMAction>(new FEMGPUEmbeddedMesh(_mesh));
}
bool FEMGPUEmbeddedMesh::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    {
        //assemble coord
        Cold coordD;
        _sys->getPosL(coordD);

        Colf coord(_szSys*4);
        for(sizeType i=0; i<_szSys; i++)
            coord.block<4,1>(i*4,0).setConstant((GLfloat)coordD[i]);

        glBindBuffer(GL_TEXTURE_BUFFER_EXT,_ctbuf);
        GLfloat* data=(GLfloat*)glMapBuffer(GL_TEXTURE_BUFFER_EXT,GL_READ_WRITE);
        Eigen::Map<Colf>(data,_szSys*4)=coord.cast<GLfloat>();
        glUnmapBuffer(GL_TEXTURE_BUFFER_EXT);
    }
    GLint nrP;
    {
        //find no. vertex
        glBindBuffer(GL_ARRAY_BUFFER_ARB,_vbuf0);
        glGetBufferParameteriv(GL_ARRAY_BUFFER_ARB,GL_BUFFER_SIZE,&nrP);
        nrP/=sizeof(GLfloat)*3;
    }
    //draw
    _shaderRecon->begin();
    {
        //setup shader
        _shaderRecon->setUniform1i("nrBasis",(GLint)_szSys);
        _shaderRecon->setUniform1i("colOff",(GLint)nrP*4);

        glActiveTextureARB(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER_EXT,_ctid);
        _shaderRecon->setUniform1i("coordR",0);

        glActiveTextureARB(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER_EXT,_btid);
        _shaderRecon->setUniform1i("basisR",1);
    }
    {
        //render to vbo
        glEnable(GL_RASTERIZER_DISCARD);
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,0,_vbuf1);
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER,1,_nbuf1);
        glBeginTransformFeedback(GL_POINTS);
        {
            glEnableClientState(GL_VERTEX_ARRAY);
            glBindBufferARB(GL_ARRAY_BUFFER_ARB,_vbuf0);
            glVertexPointer(3,GL_FLOAT,0,0);

            glEnableClientState(GL_NORMAL_ARRAY);
            glBindBufferARB(GL_ARRAY_BUFFER_ARB,_nbuf0);
            glNormalPointer(GL_FLOAT,0,0);

            glDrawArrays(GL_POINTS,0,nrP);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
            glBindBufferARB(GL_TEXTURE_BUFFER_EXT,0);
        }
        glEndTransformFeedback();
        glDisable(GL_RASTERIZER_DISCARD);
    }
    _shaderRecon->end();
    return false;
}
bool FEMGPUEmbeddedMesh::draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    src->_mat._ref->setOpenGLMaterial();
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB,_vbuf1);
    glVertexPointer(3,GL_FLOAT,0,0);

    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB,_nbuf1);
    glNormalPointer(GL_FLOAT,0,0);

    GLint nrI;
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, _ibuf);
    glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER_ARB,GL_BUFFER_SIZE,&nrI);
    glDrawElements(GL_TRIANGLES,(GLsizei)nrI/sizeof(GLuint),GL_UNSIGNED_INT,NULL);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
    return false;
}
bool FEMGPUEmbeddedMesh::initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    boost::shared_ptr<PVSMScene::PVSMObjMesh> msh=
        boost::dynamic_pointer_cast<PVSMScene::PVSMObjMesh>(src);
    if(!msh)
        return false;

    //build VBO
    ObjMesh mesh;
    boost::filesystem::ifstream is(msh->_file[0]);
    mesh.read(is,false,false);
    Eigen::Matrix<GLfloat,-1,1> vertices,normals;
    Eigen::Matrix<GLuint,-1,1> indices;
    GLenum type;
    StaticDraw::toArray(mesh,src,Mat4::Identity(),vertices,normals,indices,type);

    //interpolator
    vector<FEMInterp> evss(mesh.getV().size());
    _mesh.getBary(mesh.getV(),evss,numeric_limits<scalar>::max());
    sizeType bodyId=-1;
    for(sizeType i=0; i<(sizeType)evss.size(); i++)
        if(bodyId!=-1) ASSERT_MSG(bodyId == evss[i]._id,"Interbody Mesh Detected!")
            else bodyId=evss[i]._id;

    //assemble offset
    _sys=_mesh.getB(bodyId)._system;
    _szSys=_sys->sizeL();
    _sys->body().setDPos(_sys->LtoF(Cold::Zero(_szSys)));

    //build shader
    const char* vars[2]= {"pos","nor"};
    _shaderRecon=s.getRenderer().getManager().loadfromFile("ReducedRecon.txt",NULL,2,vars,GL_SEPARATE_ATTRIBS);

    //build buffer
    genBuffer(GL_ARRAY_BUFFER_ARB,GL_STATIC_DRAW,_vbuf0,vertices);
    genBuffer(GL_ARRAY_BUFFER_ARB,GL_DYNAMIC_DRAW,_vbuf1,vertices);
    genBuffer(GL_ARRAY_BUFFER_ARB,GL_STATIC_DRAW,_nbuf0,normals);
    genBuffer(GL_ARRAY_BUFFER_ARB,GL_DYNAMIC_DRAW,_nbuf1,normals);
    genBuffer(GL_ELEMENT_ARRAY_BUFFER_ARB,GL_DYNAMIC_DRAW,_ibuf,indices);

    //assemble basis
    sizeType colOff=16*evss.size();
    Colf basis=Colf::Zero(_szSys*colOff);
    for(sizeType i=0,off=0; i<_szSys; i++,off+=colOff) {
        _sys->body().setDPos(_sys->LtoF(Cold::Unit(_szSys,i)));
        for(sizeType v=0; v<(sizeType)evss.size(); v++) {
            Vec3 pos=evss[v]._cell->getVert(evss[v]);
            basis.block<3,1>(off+v*16+0,0)=pos.cast<scalarF>()-vertices.block<3,1>(v*3,0);

            Mat3 F,FN,FR;
            evss[v]._cell->buildF(F,FN,FR);
            F-=Mat3::Identity();
            basis.block<3,1>(off+v*16+4,0)=F.col(0).cast<scalarF>();
            basis.block<3,1>(off+v*16+8,0)=F.col(1).cast<scalarF>();
            basis.block<3,1>(off+v*16+12,0)=F.col(2).cast<scalarF>();
        }
    }
    Colf coord(_szSys*4);
    genTBuffer(GL_RGBA32F_ARB,_ctbuf,_ctid,coord);
    genTBuffer(GL_RGBA32F_ARB,_btbuf,_btid,basis);
    return false;
}

FEMComputeMesh::FEMComputeMesh(const FEMMesh& mesh):StaticDraw("computeMesh"),_mesh(mesh)
{
    _names["computeMesh"]=PVSMScene::OBJMESH;
}
bool FEMComputeMesh::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<_mesh.nrB(); i++) {
        FEMSystem& sys=*(_mesh.getB(i)._system);
        Cold L;
        sys.getPosL(L);
        sys.body().setDPos(sys.LtoF(L));
    }
    sizeType nrV=(sizeType)_verts.size();
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrV; i++)
        _vertices.block<3,1>(i*3,0)=_verts[i]->_pos.cast<GLfloat>();
    return false;
}
bool FEMComputeMesh::initStage4(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    //create vertex map
    sizeType j=0;
    boost::unordered_map<boost::shared_ptr<FEMVertex>,sizeType> vertMap;
    for(sizeType b=0; b<(sizeType)_mesh.nrB(); b++)
        for(sizeType i=0; i<_mesh.getB(b).nrV(); i++)
            if(_mesh.getB(b).getV(i)._surface)vertMap[_mesh.getB(b).getVPtr(i)]=j++;
    //assemble vertex
    _verts.resize(j);
    for(boost::unordered_map<boost::shared_ptr<FEMVertex>,sizeType>::const_iterator
            beg=vertMap.begin(),end=vertMap.end(); beg!=end; beg++)
        _verts[beg->second]=beg->first;
    _vertices.resize(j*3);
    //assemble index
    boost::unordered_set<Vec2i,Hash> indices;
    for(sizeType b=0; b<(sizeType)_mesh.nrB(); b++)
        for(sizeType i=0; i<_mesh.getB(b).nrC(); i++) {
            const FEMCell& c=_mesh.getB(b).getC(i);
#define ADDLINE(I,J)	\
if(c._v[I]->_surface && c._v[J]->_surface)	\
{	\
	Vec2i id(vertMap[c._v[I]],vertMap[c._v[J]]);	\
	if(id[0] > id[1])std::swap(id[0],id[1]);	\
	indices.insert(id);	\
}
            ADDLINE(0,1)
            ADDLINE(0,2)
            ADDLINE(1,2)
            ADDLINE(0,3)
            ADDLINE(1,3)
            ADDLINE(2,3)
        }

    j=0;
    _indices.resize(indices.size()*2);
    for(boost::unordered_set<Vec2i,Hash>::const_iterator
            beg=indices.begin(),end=indices.end(); beg!=end; beg++) {
        _indices[j++]=(GLuint)((*beg)[0]);
        _indices[j++]=(GLuint)((*beg)[1]);
    }
    _type=GL_LINES;
    return false;
}
