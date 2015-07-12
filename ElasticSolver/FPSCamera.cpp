#include "FPSCamera.h"
#include <GL/freeglut.h>

USE_PRJ_NAMESPACE

FPSCamera::FPSCamera():PVSMScene::PVSMAction("camera")
{
    _names["cameraLine"]=PVSMScene::LINE;
    reset();
}
bool FPSCamera::preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Vec3 dir=calcDir()+_pos;
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(_pos[0],_pos[1],_pos[2],
              dir[0], dir[1], dir[2],
              _frame.col(2)[0],_frame.col(2)[1],_frame.col(2)[2]);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
    return false;
}
bool FPSCamera::frame(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    if(_forward) {
        _pos+=s.getFrameTime()*calcDir()*_moveSpd;
        _forward=false;
    }
    if(_backward) {
        _pos-=s.getFrameTime()*calcDir()*_moveSpd;
        _backward=false;
    }
    if(_left) {
        _pos+=s.getFrameTime()*calcLeft()*_moveSpd;
        _left=false;
    }
    if(_right) {
        _pos-=s.getFrameTime()*calcLeft()*_moveSpd;
        _right=false;
    }
    return false;
}
bool FPSCamera::reshape(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int w,int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0f,(GLdouble)w/(GLdouble)h,0.01f,1000.0f);
    return false;
}
bool FPSCamera::mouse(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int button,int state,int x,int y)
{
    _lastX=x;
    _lastY=y;
    checkMode(button,state,_DIR_BUT,_motionMode);
    checkMode(button,state,_SCALE_BUT,_scaleMode);
    checkMode(button,state,_TRANS_BUT,_transMode);
    return false;
}
bool FPSCamera::mouseWheel(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int wheel,int dir,int x,int y)
{
    if(dir > 0)_moveSpd*=2.0f;
    else _moveSpd=std::max<scalar>(_moveSpd*0.5f,1.0f);
    ostringstream oss;
    oss << "Set moveSpd: " << _moveSpd;
    s.getConsole().addMsg(oss.str(),1000);
    return false;
}
bool FPSCamera::motion(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int x,int y)
{
    if(_motionMode) {
        _theta-=(x-_lastX)*_rotSpd;
        while(_theta < -180.0f)_theta+=360.0f;
        while(_theta > 180.0f)_theta-=360.0f;

        _phi-=(y-_lastY)*_rotSpd;
        _phi=std::min<scalar>(std::max<scalar>(_phi,-89.0f),89.0f);
    }
    if(_scaleMode) {
        Vec3 ctr=(_pos+calcDir());
        _pos=(_pos-ctr)*std::exp((y-_lastY)*_scaleSpd)+ctr;
    }
    if(_transMode) {
        _pos-=calcLeft()*((x-_lastX)*_transSpd);
        _pos-=calcTop()*((y-_lastY)*_transSpd);
    }
    _lastX=x;
    _lastY=y;
    return false;
}
bool FPSCamera::keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y)
{
    if(key == _FORWARD)_forward=true;
    if(key == _BACKWARD)_backward=true;
    if(key == _LEFT)_left=true;
    if(key == _RIGHT)_right=true;
    return false;
}
void FPSCamera::deactivate()
{
    _motionMode=false;
    _scaleMode=false;
    _transMode=false;
    _forward=false;
    _backward=false;
    _left=false;
    _right=false;
}
bool FPSCamera::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    reset();
    boost::shared_ptr<PVSMScene::PVSMLine> line=
        boost::dynamic_pointer_cast<PVSMScene::PVSMLine>(src);
    if(line) {
        setPos(line->_P1.get());
        setDir(line->_P2.get()-line->_P1.get());
    }
    return false;
}
//setting
void FPSCamera::reset()
{
    _frame.setIdentity();
    setPos(Vec3::Zero());
    setDir(Vec3::Unit(0));
    setMoveSpd(5.0f);
    setRotSpd(0.1f);
    setScaleSpd(0.01f);
    setTransSpd(0.01f);

    _FORWARD='w';
    _BACKWARD='s';
    _LEFT='a';
    _RIGHT='d';
    _DIR_BUT=GLUT_LEFT_BUTTON;
    _SCALE_BUT=GLUT_RIGHT_BUTTON;
    _TRANS_BUT=GLUT_MIDDLE_BUTTON;
    deactivate();
}
void FPSCamera::setUp(const Vec3& up)
{
    _frame.col(2)=up.normalized();
    int maxC;
    up.cwiseAbs().maxCoeff(&maxC);
    if(maxC == 0)_frame.col(0)=Vec3::Unit(1);
    else _frame.col(0)=Vec3::Unit(0);
    _frame.col(1)=_frame.col(2).cross(_frame.col(0)).normalized();
    _frame.col(0)=_frame.col(1).cross(_frame.col(2)).normalized();
}
void FPSCamera::setPos(const Vec3& pos)
{
    _pos=pos;
}
void FPSCamera::setDir(const Vec3& dir)
{
    Vec3 f=_frame.transpose()*dir;
    _theta=std::atan2(f[1],f[0])*180.0f/M_PI;
    _phi=std::atan(dir[2]/dir.block<2,1>(0,0).norm())*180.0f/M_PI;
}
void FPSCamera::setMoveSpd(scalar spd)
{
    _moveSpd=spd;
}
void FPSCamera::setRotSpd(scalar spd)
{
    _rotSpd=spd;
}
void FPSCamera::setScaleSpd(scalar spd)
{
    _scaleSpd=spd;
}
void FPSCamera::setTransSpd(scalar spd)
{
    _transSpd=spd;
}
Vec3 FPSCamera::calcDir() const
{
    scalar ct=cos(_theta*M_PI/180.0f);
    scalar st=sin(_theta*M_PI/180.0f);
    scalar cp=cos(_phi*M_PI/180.0f);
    scalar sp=sin(_phi*M_PI/180.0f);
    return _frame*Vec3(ct*cp,st*cp,sp);
}
Vec3 FPSCamera::calcLeft() const
{
    return _frame.col(2).cross(calcDir()).normalized();
}
Vec3 FPSCamera::calcTop() const
{
    return calcDir().cross(calcLeft());
}
