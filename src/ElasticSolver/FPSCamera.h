#ifndef FPS_CAMERA_H
#define FPS_CAMERA_H

#include "MathBasic.h"
#include "PVSMScene.h"

PRJ_BEGIN

class FPSCamera : public PVSMScene::PVSMAction
{
public:
    FPSCamera();
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
        return self;
    }
    virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool frame(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool reshape(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int w,int h);
    virtual bool mouse(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int button,int state,int x,int y);
    virtual bool motion(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int x,int y);
    virtual bool mouseWheel(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int wheel,int dir,int x,int y);
    virtual bool keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y);
    virtual void deactivate();
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    //setting
    void reset();
    void setUp(const Vec3& up);
    void setPos(const Vec3& pos);
    void setDir(const Vec3& dir);
    void setMoveSpd(scalar spd);
    void setRotSpd(scalar spd);
    void setScaleSpd(scalar spd);
    void setTransSpd(scalar spd);
private:
    Vec3 calcDir() const;
    Vec3 calcLeft() const;
    Vec3 calcTop() const;
    Vec3 _pos;
    Mat3 _frame;
    scalar _theta,_phi;
    scalar _moveSpd,_rotSpd,_scaleSpd,_transSpd;
    //key binding
    int _lastX,_lastY;
    char _FORWARD,_BACKWARD,_LEFT,_RIGHT;
    unsigned int _DIR_BUT,_SCALE_BUT,_TRANS_BUT;
    bool _motionMode,_scaleMode,_transMode;
    bool _forward,_backward,_left,_right;
};

PRJ_END

#endif