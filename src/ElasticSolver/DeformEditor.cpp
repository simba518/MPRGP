#include "DeformEditor.h"
#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMCollision.h"
#include "FEMReducedSystem.h"
#include<iomanip>

USE_PRJ_NAMESPACE

DeformEditor::DeformEditor(FEMSolver& sol,const boost::property_tree::ptree& bt)
    :PVSMScene::PVSMAction("deformer"),_sol(sol),_mesh(sol.getMesh()),_geom(new FEMGeom(3))
{
    _bt=bt;
    _sol._geom=_geom;
    _names["deformer"]=PVSMScene::OBJMESH;
    _names["constraint"]=PVSMScene::BOX;
    _names["gravity"]=PVSMScene::LINE;

    _names["geomPlane"]=PVSMScene::PLANE;
    _names["geomBox"]=PVSMScene::BOX;
    _names["geomSphere"]=PVSMScene::SPHERE;
    _names["geomObjMesh"]=PVSMScene::OBJMESH;
    reset();
}
boost::shared_ptr<PVSMScene::PVSMAction> DeformEditor::copy(boost::shared_ptr<PVSMAction> self)
{
    return self;
}
bool DeformEditor::draw(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    GLsizei nrE=(GLsizei)_fixEnergy.size();
    if(_mat && nrE > 0) {
        _mat->setOpenGLMaterial();
        Vec3 pos;
        vector<GLfloat> vss(nrE*6);
        for(sizeType i=0,j=0; i<nrE; i++) {
            const FixPointEnergy& e=*(_fixEnergy[i]);
            pos=e.getCell()->getVert(e.getI());
            vss[j++]=(float)pos[0];
            vss[j++]=(float)pos[1];
            vss[j++]=(float)pos[2];

            pos=e.getPos();
            vss[j++]=(float)pos[0];
            vss[j++]=(float)pos[1];
            vss[j++]=(float)pos[2];
        }
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3,GL_FLOAT,0,&(vss[0]));
        glDrawArrays(GL_LINES,0,nrE*2);
        glDrawArrays(GL_POINTS,0,nrE*2);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    return false;
}
bool DeformEditor::frame(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    ostringstream oss;
    oss << "FPS: ";
    if(_doSim) {
        boost::timer::cpu_timer timer;
        timer.start();
        _sol.advance(s.getFrameTime());
        timer.stop();
        oss << setprecision(5) << (1000000000.0L/timer.elapsed().wall);
    } else oss << setprecision(5) << (1000.0f/s.getFPS());
    if(_sol._tree.get<bool>("showFPS",true))
        s.getConsole().addMsg(oss.str(),-1);
    else s.getConsole().rmvMsg(-1);
    return false;
}
bool DeformEditor::reshape(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int w,int h)
{
    _sz=POS(w,h);
    return false;
}
bool DeformEditor::mouse(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int button,int state,int x,int y)
{
    bool dragMode,fixMode,delMode;
    checkMode(button,state,_DRAG_BUT,dragMode);
    checkMode(button,state,_FIX_BUT,fixMode);
    checkMode(button,state,_DEL_BUT,delMode);
    _curr=POS(x,y);

    //if delete mode
    buildFrame();
    updateEnergy(dragMode,fixMode,delMode);
    return false;
}
bool DeformEditor::motion(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int x,int y)
{
    if(_currEnergy)
        _currEnergy->setPos(((_c+_z-_x+_y)+
                             (_x*2)*(x/(scalar)_sz[0])-
                             (_y*2)*(y/(scalar)_sz[1])).cast<scalar>());
    return false;
}
bool DeformEditor::mouseWheel(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,int wheel,int dir,int x,int y)
{
    if(dir > 0)_K*=2.0f;
    else _K=std::max<scalar>(_K*0.5f,1.0f);
    ostringstream oss;
    oss << "Current Constraint K: " << _K;
    s.getConsole().addMsg(oss.str(),1000);
    return false;
}
bool DeformEditor::keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y)
{
    if(key == 's') {
        _doSim=!_doSim;
        s.getConsole().addMsg(_doSim?"Start Simulation":"Stop Simulation",1000);
    } else if(key == 'r') {
        for(sizeType i=0; i<_sol.getMesh().nrB(); i++)
            _sol.resetBody(i);
        s.getConsole().addMsg("Reset All Body",1000);
    }
    return false;
}
void DeformEditor::deactivate()
{
    updateEnergy(false,false,false);
}
bool DeformEditor::initStage1(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Mat4 T=PVSMScene::getTrans(src);
    if(isInGroup(*src,"gravity")) {
        boost::shared_ptr<PVSMScene::PVSMLine> line=
            boost::dynamic_pointer_cast<PVSMScene::PVSMLine>(src);
        Vec3 grav=PVSMScene::transform(T,line->_P2.get())-
                  PVSMScene::transform(T,line->_P1.get());
        _bt.put<scalar>("gravX",grav[0]);
        _bt.put<scalar>("gravY",grav[1]);
        _bt.put<scalar>("gravZ",grav[2]);
        _mat=line->_mat._ref;
    }
    if(isInGroup(*src,"deformer")) {
        boost::shared_ptr<PVSMScene::PVSMObjMesh> mesh=
            boost::dynamic_pointer_cast<PVSMScene::PVSMObjMesh>(src);
        boost::filesystem::path path(mesh->_file[0]);
        path.replace_extension(".abq");

        FEMMesh tmpMesh(3,_sol.getMesh().getColl().copy());
        tmpMesh.reset(path.string(),0.0f);
        tmpMesh.applyTrans(T,0,true,true);
        tmpMesh.getB(0)._tree=_bt;
        tmpMesh.getB(0)._system.reset(new FEMReducedSystem(tmpMesh.getB(0)));
        boost::dynamic_pointer_cast<FEMReducedSystem>(tmpMesh.getB(0)._system)->readEnergy(path.replace_extension(".elastic").string());
        _mesh+=tmpMesh;
        //_mesh.writeBVH();
    }
    if(isInGroup(*src,"geomPlane")) {
        boost::shared_ptr<PVSMScene::PVSMPlane> plane=
            boost::dynamic_pointer_cast<PVSMScene::PVSMPlane>(src);
        _geom->addGeomPlane(T,plane->getPlane());
    }
    if(isInGroup(*src,"geomBox")) {
        boost::shared_ptr<PVSMScene::PVSMBox> box=
            boost::dynamic_pointer_cast<PVSMScene::PVSMBox>(src);
        _geom->addGeomBox(T,box->getBox());
    }
    if(isInGroup(*src,"geomSphere")) {
        boost::shared_ptr<PVSMScene::PVSMSphere> sphere=
            boost::dynamic_pointer_cast<PVSMScene::PVSMSphere>(src);
        _geom->addGeomSphere(PVSMScene::transform(T,sphere->_Ctr.get()),sphere->_Rad[0]);
    }
    if(isInGroup(*src,"geomObjMesh")) {
        boost::shared_ptr<PVSMScene::PVSMObjMesh> objMesh=
            boost::dynamic_pointer_cast<PVSMScene::PVSMObjMesh>(src);
        _geom->addGeomMesh(T,objMesh->_file[0]);
    }
    return true;
}
bool DeformEditor::initStage2(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    Mat4 T=PVSMScene::getTrans(src);
    if(isInGroup(*src,"constraint")) {
        boost::shared_ptr<PVSMScene::PVSMBox> box=
            boost::dynamic_pointer_cast<PVSMScene::PVSMBox>(src);
        FEMGeom geom(3);
        geom.addGeomBox(T,box->getBox());
        geom.assemble();
        _sol.addConstraintPoint(geom);
    }
    return true;
}
bool DeformEditor::initStage3(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src)
{
    _geom->assemble();
    _sol.resetImplicitEuler(1E-4f,1);
    _sol.setSelfColl(false);
    _sol.setCollK(1E6f);
    for(sizeType i=0; i<_sol.getMesh().nrB(); i++) {
        _sol.getSystem<FEMReducedSystem>(i).addEnergyMass();
        _sol.getSystem<FEMReducedSystem>(i).buildU();
        _sol.getSystem<FEMReducedSystem>(i).onDirty();
    }
    return false;
}
void DeformEditor::reset()
{
    _DRAG_BUT=GLUT_LEFT_BUTTON;
    _FIX_BUT=GLUT_MIDDLE_BUTTON;
    _DEL_BUT=GLUT_RIGHT_BUTTON;

    setRadius(0.01f);
    setConstraintK(1E3f);
    deactivate();
    _doSim=true;
}
void DeformEditor::buildFrame()
{
    static GLdouble modelView[16];
    static GLdouble projection[16];
    glGetDoublev(GL_MODELVIEW_MATRIX,modelView);
    glGetDoublev(GL_PROJECTION_MATRIX,projection);
    getViewMatrixFrame<scalarD>(Eigen::Map<Mat4d>(modelView),
                                Eigen::Map<Mat4d>(projection),
                                _c,_x,_y,_z);
}
void DeformEditor::setDist(Vec3 pt)
{
    Vec3d ptd=pt.cast<scalarD>();
    scalarD coef=(ptd-_c).dot(_z)/_z.dot(_z);
    _x*=coef;
    _y*=coef;
    _z*=coef;
}
void DeformEditor::updateEnergy(bool dragMode,bool fixMode,bool delMode)
{
    LineSeg l(_c.cast<scalar>(),
              ((_c+_z-_x+_y)+
               (_x*2)*(_curr[0]/(scalar)_sz[0])-
               (_y*2)*(_curr[1]/(scalar)_sz[1])).cast<scalar>());
    l._y=(l._y-l._x).normalized()*10000.0f+l._x;

    if(delMode) {
        scalar sqrDist;
        Vec3 cp,b;
        for(sizeType i=0; i<(sizeType)_fixEnergy.size();) {
            l.calcPointDist(_fixEnergy[i]->getPos(),sqrDist,cp,b);
            if(sqrDist < _rad*_rad) {
                _sol.delEnergy(_fixEnergy[i]);
                if(_fixEnergy[i] == _currEnergy)
                    _currEnergy.reset((FixPointEnergy*)NULL);
                _fixEnergy[i]=_fixEnergy.back();
                _fixEnergy.pop_back();
            } else i++;
        }
    }
    //if in editing mode, reset current energy
    if(dragMode||fixMode) {
        FEMInterp cd;
        scalar dist=numeric_limits<scalar>::max();
        _mesh.getColl().rayCastMesh(l,cd,dist);
        if(dist < numeric_limits<scalar>::max()) {
            boost::shared_ptr<FEMCell> c=cd._cell;
            _currEnergy.reset(new FixPointEnergy(3,c,cd,_K,c->getVert(cd)));
            _fixEnergy.push_back(_currEnergy);
            _sol.addEnergy(_currEnergy);

            _isFix=fixMode;
            setDist(l._x+(l._y-l._x).normalized()*dist);
        }
    }
    //else delete current energy
    else if(_currEnergy) {
        if(!_isFix) {
            _fixEnergy.erase(std::find(_fixEnergy.begin(),_fixEnergy.end(),_currEnergy));
            _sol.delEnergy(_currEnergy);
        }
        _currEnergy.reset((FixPointEnergy*)NULL);
    }
}
