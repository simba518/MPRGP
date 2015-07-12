#ifndef DEFORM_EDITOR_H
#define DEFORM_EDITOR_H

#include "MathBasic.h"
#include "FEMSystem.h"
#include "PVSMScene.h"
#include <GL/freeglut.h>

PRJ_BEGIN

class DeformEditor : public PVSMScene::PVSMAction
{
public:
    typedef Eigen::Vector2i POS;
    DeformEditor(FEMSolver& sol,const boost::property_tree::ptree& bt);
    virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self);
    virtual bool draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool frame(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool reshape(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int w,int h);
    virtual bool mouse(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int button,int state,int x,int y);
    virtual bool motion(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int x,int y);
    virtual bool mouseWheel(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int wheel,int dir,int x,int y);
    virtual bool keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y);
    virtual void deactivate();
    virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage2(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    virtual bool initStage3(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src);
    void setRadius(scalar rad) {
        _rad=rad;
    }
    void setConstraintK(scalar K) {
        _K=K;
    }
    void reset();
protected:
    void buildFrame();
    void setDist(Vec3 pt);
    void updateEnergy(bool dragMode,bool fixMode,bool delMode);
    //solver
    FEMSolver& _sol;
    FEMMesh& _mesh;
    boost::shared_ptr<FEMGeom> _geom;
    //key binding
    int _DRAG_BUT,_FIX_BUT,_DEL_BUT;
    scalar _rad,_K;
    //temporary
    POS _sz,_curr;
    Vec3d _c,_x,_y,_z;
    vector<boost::shared_ptr<FixPointEnergy> > _fixEnergy;
    boost::shared_ptr<FixPointEnergy> _currEnergy;
    boost::shared_ptr<PVSMScene::PVSMMaterial> _mat;
    boost::property_tree::ptree _bt;
    bool _isFix,_doSim;
};

PRJ_END

#endif