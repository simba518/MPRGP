#include "Console.h"
#include <time.h>
#include <GL/freeglut.h>

USE_PRJ_NAMESPACE

DefaultConsole::DefaultConsole()
{
    setColor(1.0f,0.5f,0.5f,1.0f);
}
void DefaultConsole::setColor(float r,float g,float b,float a)
{
    _r=r;
    _g=g;
    _b=b;
    _a=a;
}
void DefaultConsole::addMsg(const std::string& msg,sizeType dur)
{
    if(dur > 0)
        _msgs.push(make_pair(msg,(sizeType)clock()+dur));
    else if(dur < 0)
        _cmsgs[dur]=msg;
}
void DefaultConsole::rmvMsg(sizeType dur)
{
    if(dur < 0)
        _cmsgs.erase(dur);
}
void DefaultConsole::reshape(int w,int h)
{
    _w=w;
    _h=h;
}
void DefaultConsole::render()
{
    initDrawMsg();
    sizeType curr=(sizeType)clock();
    for(std::map<sizeType,std::string>::const_iterator
            beg=_cmsgs.begin(),end=_cmsgs.end(); beg!=end; beg++)
        drawMsg(beg->second);

    std::queue<std::pair<std::string,sizeType> > bk;
    while(!_msgs.empty()) {
        std::pair<std::string,sizeType> top=_msgs.front();
        _msgs.pop();
        if(top.second > curr) {
            drawMsg(top.first);
            bk.push(top);
        }
    }
    bk.swap(_msgs);
    finishDrawMsg();
}
void DefaultConsole::initDrawMsg()
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0f,1.0f,0.0f,1.0f);

    glDisable(GL_LIGHTING);
    glColor4f(_r,_g,_b,_a);

    GLfloat delta=0.01f;
    _posx=delta;
    _posy=1.0f-delta*5.0f;
}
void DefaultConsole::finishDrawMsg()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}
void DefaultConsole::drawMsg(const std::string& str)
{
    glRasterPos2f(_posx,_posy);
    for(sizeType i=0; i<(sizeType)str.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,str[i]);
    _posy-=24.0f*1.2f/(GLfloat)_h;
}