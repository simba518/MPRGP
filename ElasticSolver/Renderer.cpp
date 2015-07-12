#include "glsl.h"
#include "Renderer.h"
#include "MathBasic.h"
#include <fstream>
#include <GL/freeglut.h>

USE_PRJ_NAMESPACE

//Renderable
void Renderable::drawObject()
{
    {
        GLfloat material_Ka[] = {0.5f, 0.0f, 0.0f, 1.0f};
        GLfloat material_Kd[] = {0.4f, 0.4f, 0.5f, 1.0f};
        GLfloat material_Ks[] = {0.8f, 0.8f, 0.0f, 1.0f};
        GLfloat material_Ke[] = {0.1f, 0.0f, 0.0f, 0.0f};
        GLfloat material_Se = 20.0f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_Ka);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_Kd);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_Ks);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, material_Ke);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material_Se);

        glutSolidSphere(0.5,16,16);
    }

    {
        GLfloat material_Ka[] = {0.5f, 0.5f, 0.5f, 1.0f};
        GLfloat material_Kd[] = {0.4f, 0.4f, 0.4f, 1.0f};
        GLfloat material_Ks[] = {0.8f, 0.8f, 0.0f, 1.0f};
        GLfloat material_Ke[] = {0.1f, 0.0f, 0.0f, 0.0f};
        GLfloat material_Se = 20.0f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_Ka);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_Kd);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_Ks);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, material_Ke);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material_Se);

        GLfloat y=-0.8f;
        glBegin(GL_TRIANGLES);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d(-10.0f,y,-10.0f);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d( 10.0f,y,-10.0f);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d( 10.0f,y, 10.0f);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d(-10.0f,y,-10.0f);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d( 10.0f,y, 10.0f);
        glNormal3d(0.0f,1.0f,0.0f);
        glVertex3d(-10.0f,y, 10.0f);
        glEnd();
    }
}
//DefaultRenderer
DefaultRenderer::~DefaultRenderer() {}
void DefaultRenderer::init()
{
    _SM.reset(new cwc::glShaderManager);
}
void DefaultRenderer::render(Renderable& rend)
{
    glViewport(0,0,_w,_h);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    rend.drawObject();
}
void DefaultRenderer::reshape(int w,int h)
{
    _w=w;
    _h=h;
}
void DefaultRenderer::drawLights()
{
    GLfloat pos[4];
    glPushMatrix();
    for(int i=0; i<MAX_LIGHTS; i++) {
        glColor3f(0.0f,0.0f,0.0f);
        if(glIsEnabled(GL_LIGHT0+i)) {
            glDisable(GL_LIGHTING);
            glGetLightfv(GL_LIGHT0+i,GL_POSITION,pos);
            glTranslatef(pos[0],pos[1],pos[2]);
            glutSolidSphere(_rad,16,16);
            glEnable(GL_LIGHTING);
        }
    }
    glPopMatrix();
}
void DefaultRenderer::setShowLights(bool show,GLfloat rad)
{
    _showLight=show;
    _rad=rad;
}
cwc::glShaderManager& DefaultRenderer::getManager()
{
    return *_SM;
}
//utility
GLuint DefaultRenderer::getMaxNumLight()
{
    return MAX_LIGHTS;
}
void DefaultRenderer::checkGLError()
{
    GLenum glErr;
    //int retCode=0;
    glErr=glGetError();
    while(glErr != GL_NO_ERROR) {
        const GLubyte* sError=gluErrorString(glErr);
        if(sError) std::cout << "GL Error #%d (%s)" << sError << std::endl;
        else std::cout << "GL Error #" << glErr << " (no message available)" << std::endl;
        system("pause");
        //retCode=1;
        glErr=glGetError();
    }
}
void DefaultRenderer::drawQuad()
{
    //Drawing quad
    glBegin(GL_QUADS);
    glTexCoord2d(0,0);
    glVertex3f(0,0,0);
    glTexCoord2d(1,0);
    glVertex3f(1,0,0);
    glTexCoord2d(1,1);
    glVertex3f(1,1,0);
    glTexCoord2d(0,1);
    glVertex3f(0,1,0);
    glEnd();
}
void DefaultRenderer::checkFBO()
{
    // check FBO status
    GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
        WARNING("GL_FRAMEBUFFER_COMPLETE_EXT failed, CANNOT use FBO");
}
void DefaultRenderer::genTexture(GLuint& tid,GLsizei w,GLsizei h,GLint internalFormat,GLenum format,GLenum type,GLint filter,GLint wrap)
{
    glGenTextures(1,&tid);
    glBindTexture(GL_TEXTURE_2D,tid);
    // GL_LINEAR does not make sense for depth texture. However, next tutorial shows usage of GL_LINEAR and PCF
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,filter);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,filter);
    // Remove artefact on the edges of the shadowmap
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,wrap);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,wrap);
    //glTexParameterfv( GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor );
    // No need to force GL_DEPTH_COMPONENT24, drivers usually give you the max precision if available
    glTexImage2D(GL_TEXTURE_2D,0,internalFormat,w,h,0,format,type,0);
    glBindTexture(GL_TEXTURE_2D,0);
}
void DefaultRenderer::writeTexture(GLuint tid,GLenum format,std::string path)
{
    GLint w,h;
    glBindTexture(GL_TEXTURE_2D,tid);
    glGetTexLevelParameteriv(GL_TEXTURE_2D,0,GL_TEXTURE_WIDTH,&w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D,0,GL_TEXTURE_HEIGHT,&h);

    std::vector<GLubyte> pixels(w*h,0);
    glGetTexImage(GL_TEXTURE_2D,0,format,GL_UNSIGNED_BYTE,&(pixels[0]));
    std::ofstream os(path.c_str());
    os << "P5\n" << w << " " << h << "\n255\n";
    os.write((char*)&(pixels[0]),pixels.size());
    checkGLError();
}
void DefaultRenderer::writeColorTexture(GLuint tid,std::string path)
{
    writeTexture(tid,GL_RED,std::string(path).append("Red.ppm"));
    writeTexture(tid,GL_BLUE,std::string(path).append("Blue.ppm"));
    writeTexture(tid,GL_GREEN,std::string(path).append("Green.ppm"));
}
//PerPixelRenderer
PerPixelRenderer::~PerPixelRenderer()
{
    release();
}
void PerPixelRenderer::init()
{
    DefaultRenderer::init();
    _shaderDraw=_SM->loadfromFile("PerPixelVert.txt","PerPixelFrag.txt");
    _showLight=true;
    _released=true;
    _shadowRatio=1;
}
void PerPixelRenderer::render(Renderable& rend)
{
    if(_released) {
        DefaultRenderer::render(rend);
        return;
    }

    glEnable(GL_TEXTURE_2D);
    //render shadow map
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    for(int i=0; i<MAX_LIGHTS; i++)
        if(glIsEnabled(GL_LIGHT0+i)) {
            glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_fboId[i]);
            glViewport(0,0,_w*_shadowRatio,_h*_shadowRatio);
            setupMatrices(i);
            setTextureMatrix(i);
            glCullFace(GL_FRONT);
            beforeShadowDraw(i);
            rend.drawObject();
            glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
            afterShadowDraw(i);
        }
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    //render scene
    glViewport(0,0,_w,_h);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    if(_showLight)drawLights();
    if(_shaderDraw)_shaderDraw->begin();
    beforeShadowRender(_shaderDraw);
    glCullFace(GL_BACK);
    rend.drawObject();
    if(_shaderDraw)_shaderDraw->end();
}
void PerPixelRenderer::reshape(int w,int h)
{
    DefaultRenderer::reshape(w,h);
    release();
    for(int i=0; i<MAX_LIGHTS; i++) {
        genTexture(_tid[i],w,h,GL_DEPTH_COMPONENT,GL_DEPTH_COMPONENT,GL_UNSIGNED_BYTE);
        //create a framebuffer object
        glGenFramebuffersEXT(1,&(_fboId[i]));
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_fboId[i]);
        //Instruct openGL that we won't bind a color texture with the currently binded FBO
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
        //attach the texture to FBO depth attachment point
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_DEPTH_ATTACHMENT_EXT,GL_TEXTURE_2D,_tid[i],0);
        checkFBO();
        //switch back to window-system-provided framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    }
    _released=false;
}
void PerPixelRenderer::release()
{
    if(!_released)
        for(int i=0; i<MAX_LIGHTS; i++) {
            glDeleteTextures(1,&_tid[i]);
            glDeleteFramebuffersEXT(1,&_fboId[i]);
        }
    _released=true;
}
void PerPixelRenderer::setupMatrices(int L) const
{
    GLfloat pos[4],dir[3];
    glGetLightfv(GL_LIGHT0+L,GL_POSITION,pos);
    glGetLightfv(GL_LIGHT0+L,GL_SPOT_DIRECTION,dir);
    dir[0]+=pos[0];
    dir[1]+=pos[1];
    dir[2]+=pos[2];

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(100.0f,(GLdouble)_w/(GLdouble)_h,0.1f,1000.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(pos[0],pos[1],pos[2],
              dir[0],dir[1],dir[2],
              0,1,0);
}
void PerPixelRenderer::setTextureMatrix(int L) const
{
    static GLdouble modelView[16];
    static GLdouble projection[16];

    // This is matrix transform every coordinate x,y,z
    // x = x* 0.5 + 0.5
    // y = y* 0.5 + 0.5
    // z = z* 0.5 + 0.5
    // Moving from unit cube [-1,1] to [0,1]
    const GLdouble bias[16] = {
        0.5, 0.0, 0.0, 0.0,
        0.0, 0.5, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
        0.5, 0.5, 0.5, 1.0
    };

    // Grab modelview and transformation matrices
    glGetDoublev(GL_MODELVIEW_MATRIX,modelView);
    glGetDoublev(GL_PROJECTION_MATRIX,projection);

    glMatrixMode(GL_TEXTURE);
    glActiveTextureARB(GL_TEXTURE0+L);

    glLoadIdentity();
    glLoadMatrixd(bias);

    // concatating all matrice into one.
    glMultMatrixd(projection);
    glMultMatrixd(modelView);

    // Go back to normal matrix mode
    glMatrixMode(GL_MODELVIEW);
}
void PerPixelRenderer::beforeShadowDraw(int L)
{
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
    glClear(GL_DEPTH_BUFFER_BIT);
}
void PerPixelRenderer::beforeShadowRender(cwc::glShader* shader)
{
    GLint tids[MAX_LIGHTS];
    for(int i=0; i<MAX_LIGHTS; i++) {
        tids[i]=i;
        glActiveTextureARB(GL_TEXTURE0+i);
        glBindTexture(GL_TEXTURE_2D,_tid[i]);
    }
    if(shader)shader->setUniform1iv("STex",MAX_LIGHTS,tids);
    if(shader)shader->setUniform1f("zOff",0.001f);
}
//VSMPerPixelRenderer
void VSMPerPixelRenderer::init()
{
    PerPixelRenderer::init();
    _shaderDepth=_SM->loadfromFile("VSMDepthVert.txt","VSMDepthFrag.txt");
    _shaderBlur=_SM->loadfromFile("VSMDepthVert.txt","VSMBlurFrag.txt");
    generateBlurKernel(3);
}
void VSMPerPixelRenderer::reshape(int w,int h)
{
    DefaultRenderer::reshape(w,h);
    release();
    for(int i=0; i<MAX_LIGHTS; i++) {
        genTexture(_tid[i],w,h,GL_DEPTH_COMPONENT,GL_DEPTH_COMPONENT,GL_UNSIGNED_BYTE);
        genTexture(_ctid[i],w,h,GL_RGB32F_ARB,GL_RGB,GL_FLOAT,GL_LINEAR);
        //create a framebuffer object
        glGenFramebuffersEXT(1,&(_fboId[i]));
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_fboId[i]);
        //Instruct openGL that we won't bind a color texture with the currently binded FBO
        //glDrawBuffer(GL_NONE);
        //glReadBuffer(GL_NONE);
        //attach the texture to FBO depth attachment point
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_DEPTH_ATTACHMENT_EXT,GL_TEXTURE_2D,_tid[i],0);
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,_ctid[i],0);
        checkFBO();
        //switch back to window-system-provided framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    }
    {
        genTexture(_ctTmp,w,h,GL_RGB32F_ARB,GL_RGB,GL_FLOAT,GL_LINEAR);
        //create a framebuffer object
        glGenFramebuffersEXT(1,&_fboTmp);
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_fboTmp);
        //Instruct openGL that we won't bind a color texture with the currently binded FBO
        //glDrawBuffer(GL_NONE);
        //glReadBuffer(GL_NONE);
        //attach the texture to FBO depth attachment point
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,_ctTmp,0);
        checkFBO();
        //switch back to window-system-provided framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    }
    _released=false;
}
void VSMPerPixelRenderer::release()
{
    if(!_released)
        for(int i=0; i<MAX_LIGHTS; i++) {
            glDeleteTextures(1,&_tid[i]);
            glDeleteTextures(1,&_ctid[i]);
            glDeleteFramebuffersEXT(1,&_fboId[i]);
        }
    glDeleteTextures(1,&_ctTmp);
    glDeleteFramebuffersEXT(1,&_fboTmp);
    _released=true;
}
void VSMPerPixelRenderer::beforeShadowDraw(int L)
{
    if(_shaderDepth)_shaderDepth->begin();
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_FALSE);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}
void VSMPerPixelRenderer::afterShadowDraw(int L)
{
    if(_shaderDepth)_shaderDepth->end();
    //begin blur
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0f,1.0f,0.0f,1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //blur horizontally
    blurPass(1.0f/(_w*_shadowRatio),0.0f,_fboTmp,_ctid[L]);
    //blur vertically
    blurPass(0.0f,1.0f/(_h*_shadowRatio),_fboId[L],_ctTmp);
    //end blur
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    checkGLError();
}
void VSMPerPixelRenderer::blurPass(GLfloat dx,GLfloat dy,GLuint fbo,GLuint color) const
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo);
    glViewport(0,0,_w*_shadowRatio,_h*_shadowRatio);
    if(_shaderBlur)_shaderBlur->begin();
    if(_shaderBlur)_shaderBlur->setUniform2f("scaleU",dx,dy);
    if(_shaderBlur)_shaderBlur->setUniform1fv("offB",(GLsizei)_off.size(),(GLfloat*)&(_off[0]));
    if(_shaderBlur)_shaderBlur->setUniform1fv("coefB",(GLsizei)_coef.size(),(GLfloat*)&(_coef[0]));
    if(_shaderBlur)_shaderBlur->setUniform1i("nrB",(GLint)_off.size());
    if(_shaderBlur)_shaderBlur->setUniform1i("motions",0);
    glActiveTextureARB(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D,color);
    drawQuad();
    if(_shaderBlur)_shaderBlur->end();
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
}
void VSMPerPixelRenderer::generateBlurKernel(GLint nr)
{
    GLfloat sigma=(GLfloat)nr*3.0f,total=0.0f;
    _off.clear();
    _coef.clear();
    for(GLint i=-nr; i<=nr; i++) {
        GLfloat fi=(GLfloat)i;
        _off.push_back(fi);
        _coef.push_back(std::exp(-0.5f*(fi*fi)/(sigma*sigma)));
        total+=_coef.back();
    }
    for(GLint i=0; i<(GLint)_coef.size(); i++)
        _coef[i]/=total;
    /*#define NR_BLUR_MAX 20
    	_off.assign(NR_BLUR_MAX,0.0f);
    	_coef.assign(NR_BLUR_MAX,0.0f);
    	_off[0]=-3.0f;_coef[0]=0.015625f;
    	_off[1]=-2.0f;_coef[1]=0.09375f;
    	_off[2]=-1.0f;_coef[2]=0.234375f;
    	_off[3]= 0.0f;_coef[3]=0.3125f;
    	_off[4]= 1.0f;_coef[4]=0.234375f;
    	_off[5]= 2.0f;_coef[5]=0.09375f;
    	_off[6]= 3.0f;_coef[6]=0.015625f;
    #undef NR_BLUR_MAX*/
}
void VSMPerPixelRenderer::beforeShadowRender(cwc::glShader* shader)
{
    GLint tids[MAX_LIGHTS];
    for(int i=0; i<MAX_LIGHTS; i++) {
        tids[i]=i;
        glActiveTextureARB(GL_TEXTURE0+i);
        glBindTexture(GL_TEXTURE_2D,_ctid[i]);
    }
    if(shader)shader->setUniform1iv("STex",MAX_LIGHTS,tids);
    if(shader)shader->setUniform1f("zOff",0.0f);
}
