#include "PVSMScene.h"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/tokenizer.hpp>

#include <GL/freeglut.h>

USE_PRJ_NAMESPACE

//PVSMTransform3
Mat4 PVSMScene::PVSMTransform3::get() const
{
    Mat4 T=Mat4::Identity();
    T.block<3,1>(0,3)=getTranslate();
    T.block<3,3>(0,0).diagonal()=getScale();
    T.block<3,3>(0,0)=getRotate()*T.block<3,3>(0,0);
    return T;
}
Vec3 PVSMScene::PVSMTransform3::getScale() const
{
    return _S.get();
}
Mat3 PVSMScene::PVSMTransform3::getRotate() const
{
    Mat3 ret=Mat3::Identity();
    ret=getQuat(_R[2]*M_PI/180.0f,ret*Vec3::Unit(2))*ret;
    ret=getQuat(_R[0]*M_PI/180.0f,ret*Vec3::Unit(0))*ret;
    ret=getQuat(_R[1]*M_PI/180.0f,ret*Vec3::Unit(1))*ret;
    return ret;
}
Vec3 PVSMScene::PVSMTransform3::getTranslate() const
{
    return _P.get();
}
Mat3 PVSMScene::PVSMTransform3::getQuat(scalar rot,const Vec3& axis)
{
    scalar c=cos(0.5f*rot);
    scalar s=sin(0.5f*rot);
    scalar t=s/axis.norm();
    return Eigen::Quaternion<scalar>(c,axis[0]*t,axis[1]*t,axis[2]*t).toRotationMatrix();
}

//PVSMMaterial
#define PGL(T,A) #T#A,A
PVSMScene::PVSMPSize::PVSMPSize():PVSMMat<1>(PGL(MATERIAL,GL_POINT_SIZE)) {}
void PVSMScene::PVSMPSize::setOpenGL() const
{
    glPointSize(_value[0]);
}
PVSMScene::PVSMLWidth::PVSMLWidth():PVSMMat<1>(PGL(MATERIAL,GL_LINE_WIDTH)) {}
void PVSMScene::PVSMLWidth::setOpenGL() const
{
    glLineWidth(_value[0]);
}
PVSMScene::PVSMMaterial::PVSMMaterial():_type("type"),
//material
    _MAmb(PGL(MATERIAL,GL_AMBIENT)),
    _MDiff(PGL(MATERIAL,GL_DIFFUSE)),
    _MSpec(PGL(MATERIAL,GL_SPECULAR)),
    _MEmi(PGL(MATERIAL,GL_EMISSION)),
    _MShi(PGL(MATERIAL,GL_SHININESS)),
//light
    _LPos(PGL(LIGHT,GL_POSITION)),
    _LDir(PGL(LIGHT,GL_SPOT_DIRECTION)),
    _LAmb(PGL(LIGHT,GL_AMBIENT)),
    _LDiff(PGL(LIGHT,GL_DIFFUSE)),
    _LSpec(PGL(LIGHT,GL_SPECULAR)),
    _LSpotE(PGL(LIGHT,GL_SPOT_EXPONENT)),
    _LCutOff(PGL(LIGHT,GL_SPOT_CUTOFF)),
    _LAC0(PGL(LIGHT,GL_CONSTANT_ATTENUATION)),
    _LAC1(PGL(LIGHT,GL_LINEAR_ATTENUATION)),
    _LAC2(PGL(LIGHT,GL_QUADRATIC_ATTENUATION))
{
    //randomInit();
}
void PVSMScene::PVSMMaterial::randomInit(const PVSMSource* src)
{
    _type._value="Material";
    if(dynamic_cast<const PVSMLine*>(src) || dynamic_cast<const PVSMPoint*>(src)) {
        _MAmb.set<scalarF>(Vec4f(0.44f,0.72f,0.72f,1.0f));
        _MDiff.set<scalarF>(Vec4f(0.44f,0.72f,0.72f,1.0f));
    } else if(dynamic_cast<const PVSMObjMesh*>(src)) {
        _MAmb.set<scalarF>(Vec4f(0.6f,0.46f,0.3f,1.0f));
        _MDiff.set<scalarF>(Vec4f(0.6f,0.46f,0.3f,1.0f));
    } else {
        _MAmb.set<scalarF>(Vec4f(0.1f,0.1f,0.1f,1.0f));
        _MDiff.set<scalarF>(Vec4f(0.4f,0.4f,0.4f,1.0f));
    }
    _MSpec.set<scalarF>(Vec4f(1.0f,1.0f,1.0f,1.0f));
    _MEmi.set<scalarF>(Vec4f(0.3f,0.3f,0.5f,0.0f));
    _MShi.set(50.0f);
    _MLWidth.set(3.0f);
    _MPSize.set(8.0f);

    _LAmb.set<scalarF>(Vec4f(0.2f,0.2f,0.2f,1.0f));
    _LDiff.set<scalarF>(Vec4f(0.6f,0.6f,0.6f,1.0f));
    _LSpec.set<scalarF>(Vec4f(1.0f,0.6f,0.6f,1.0f));

    _LSpotE.set(2.0f);
    _LCutOff.set(45.0f);
    _LAC0.set(1.0f);
    _LAC1.set(0.0f);
    _LAC2.set(0.0f);
}
bool PVSMScene::PVSMMaterial::read(const std::string& name,const ELEMENT& e)
{
    return PVSMObject::read(name,e) && _type.read(e,"Material") &&

           _MAmb.read(e,_id._value) && _MDiff.read(e,_id._value) &&
           _MSpec.read(e,_id._value) && _MEmi.read(e,_id._value) &&
           _MShi.read(e,_id._value) && _MLWidth.read(e,_id._value) && _MPSize.read(e,_id._value) &&

           _LAmb.read(e,_id._value) && _LDiff.read(e,_id._value) && _LSpec.read(e,_id._value) &&
           _LSpotE.read(e,_id._value) && _LCutOff.read(e,_id._value) &&
           _LAC0.read(e,_id._value) && _LAC1.read(e,_id._value) && _LAC2.read(e,_id._value);
}
std::string PVSMScene::PVSMMaterial::write(ELEMENT& e) const
{
    _type.write(e);
    _MAmb.write(e,_id._value);
    _MDiff.write(e,_id._value);
    _MSpec.write(e,_id._value);
    _MEmi.write(e,_id._value);
    _MShi.write(e,_id._value);
    _MLWidth.write(e,_id._value);
    _MPSize.write(e,_id._value);

    _LAmb.write(e,_id._value);
    _LDiff.write(e,_id._value);
    _LSpec.write(e,_id._value);
    _LSpotE.write(e,_id._value);
    _LCutOff.write(e,_id._value);
    _LAC0.write(e,_id._value);
    _LAC1.write(e,_id._value);
    _LAC2.write(e,_id._value);

    PVSMObject::write(e);
    return "Proxy";
}
boost::shared_ptr<PVSMScene::PVSMObject> PVSMScene::PVSMMaterial::copy() const
{
    return boost::shared_ptr<PVSMObject>(new PVSMMaterial);
}
void PVSMScene::PVSMMaterial::setOpenGLMaterial()
{
    _MAmb.setOpenGL();
    _MDiff.setOpenGL();
    _MSpec.setOpenGL();
    _MEmi.setOpenGL();
    _MShi.setOpenGL();
    _MLWidth.setOpenGL();
    _MPSize.setOpenGL();
}
void PVSMScene::PVSMMaterial::setOpenGLLight(PVSMScene& s)
{
    if(_lightId < s.getRenderer().getMaxNumLight()) {
        glEnable(GL_LIGHT0+_lightId);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        _LPos.setOpenGL(_lightId);
        _LDir.setOpenGL(_lightId);
        glPopMatrix();

        _LAmb.setOpenGL(_lightId);
        _LDiff.setOpenGL(_lightId);
        _LSpec.setOpenGL(_lightId);
        _LSpotE.setOpenGL(_lightId);
        _LCutOff.setOpenGL(_lightId);
        _LAC0.setOpenGL(_lightId);
        _LAC1.setOpenGL(_lightId);
        _LAC2.setOpenGL(_lightId);
    }
}
#undef PGL

//PVSMAction
bool PVSMScene::PVSMAction::isInGroup(const PVSMSource& src,const std::string& forceGroup) const
{
    typedef boost::tokenizer<boost::char_separator<char> > Token;
    boost::char_separator<char> sep(",");
    Token tok(src._name._value,sep);
    for(Token::iterator it=tok.begin(); it!=tok.end(); ++it) {
        if(!forceGroup.empty()) {
            if(forceGroup == *it)
                return true;
        } else {
            boost::unordered_map<std::string,SOURCE_TYPE>::const_iterator iter=_names.find(*it);
            if(iter != _names.end() && iter->second == src.getType())
                return true;
        }
    }
    return false;
}
void PVSMScene::PVSMAction::checkMode(int button,int state,int buttonMode,bool& mode)
{
    if(button == buttonMode) {
        if(state == GLUT_DOWN)mode=true;
        if(state == GLUT_UP)mode=false;
    } else mode=false;
}

//IO
bool PVSMScene::read(const std::string& path)
{
    _objs.clear();
    _names.clear();
    _gSetting=GlobalSetting();

    //fill candidate
    setCandidates();

    //read pass
    ELEMENT pt;
    boost::filesystem::ifstream is(path);
    boost::property_tree::read_xml(is,pt);
    if(pt.empty())
        return false;
    readElement("",pt);

    //link pass
    INFOV("Read %lu elements from paraview",_objs.size())
    for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator
            beg=_objs.begin(),end=_objs.end(); beg!=end; beg++)
        if(!beg->second->link(*this))
            return false;

    //assign name pass
    for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator
            beg=_objs.begin(),end=_objs.end(); beg!=end; beg++)
        beg->second->_name._value=findName(beg->second->_id._value);

    updateScene();
    return true;
}
bool PVSMScene::write(const std::string& path) const
{
    //write all items
    ELEMENT pt;
    ELEMENT ids;
    PVSMAttribute<std::string> idsName("name");
    idsName._value="sources";
    idsName.write(ids);
    for(boost::unordered_map<int,std::string>::const_iterator
            beg=_names.begin(),end=_names.end(); beg!=end; beg++) {
        PVSMItem item;
        item._id._value=beg->first;
        item._name._value=beg->second;

        ELEMENT id;
        ids.add_child(item.write(id),id);
    }
    pt.add_child("ProxyCollection",ids);

    //write objects
    for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator
            beg=_objs.begin(),end=_objs.end(); beg!=end; beg++) {
        ELEMENT obj;
        pt.add_child(beg->second->write(obj),obj);
    }

    //write global setting
    ELEMENT gs;
    pt.add_child(_gSetting.write(gs),gs);

    boost::property_tree::xml_writer_settings<char> settings('\t',1);
    boost::filesystem::ofstream os(path);
    boost::property_tree::write_xml(os,pt,settings);
    return true;
}
void PVSMScene::clearUnused()
{
    //discard unused
    {
        boost::unordered_map<int,boost::shared_ptr<PVSMObject> > validObjs;
        for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator
                beg=_objs.begin(),end=_objs.end(); beg!=end; beg++) {
            boost::shared_ptr<PVSMSource> src=
                boost::dynamic_pointer_cast<PVSMSource>(beg->second);
            if(!src || src->_name._value != "")
                validObjs[beg->first]=beg->second;
        }
        _objs.swap(validObjs);
    }
    while(true) {
        boost::unordered_map<int,boost::shared_ptr<PVSMObject> > validObjs;
        for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator
                beg=_objs.begin(),end=_objs.end(); beg!=end; beg++) {
            boost::shared_ptr<PVSMSource> src=
                boost::dynamic_pointer_cast<PVSMSource>(beg->second);
            if(src || beg->second.use_count() > 0)
                validObjs[beg->first]=beg->second;
        }
        if(_objs.size() == validObjs.size())
            break;
        _objs.swap(validObjs);
    }
    INFOV("Read %lu valid objects",_objs.size())
}
PVSMScene::GlobalSetting PVSMScene::getGlobalSetting()
{
    return _gSetting;
}
//group management
void PVSMScene::delAction(boost::shared_ptr<PVSMAction> action)
{
    _candidateActions.erase(action);
}
void PVSMScene::addAction(boost::shared_ptr<PVSMAction> action)
{
    _candidateActions.insert(action);
}
void PVSMScene::delSceneAction(boost::shared_ptr<PVSMAction> action)
{
    _sceneActions.erase(action);
}
void PVSMScene::addSceneAction(boost::shared_ptr<PVSMAction> action)
{
    _sceneActions.insert(action);
}
void PVSMScene::updateScene()
{
    for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::iterator
            beg=_objs.begin(),end=_objs.end(); beg!=end; beg++) {
        boost::shared_ptr<PVSMSource> src=
            boost::dynamic_pointer_cast<PVSMSource>(beg->second);
        if(src) {
            src->_actions.clear();
            for(boost::unordered_set<boost::shared_ptr<PVSMAction> >::const_iterator
                    begA=_candidateActions.begin(),endA=_candidateActions.end(); begA!=endA; begA++)
                if((*begA)->isInGroup(*src))
                    src->_actions.push_back((*begA)->copy(*begA));
        }
    }
}
//events
void PVSMScene::render()
{
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);

    visit(PRE_DRAW);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glClearColor(_gSetting._clearColor[0],_gSetting._clearColor[1],
                 _gSetting._clearColor[2],_gSetting._clearColor[3]);
    getRenderer().render(*this);

    _console->setColor(_gSetting._textColor.get()[0],_gSetting._textColor.get()[1],
                       _gSetting._textColor.get()[2],_gSetting._textColor.get()[3]);
    _console->render();
    visit(POST_DRAW);
    glutSwapBuffers();
}
void PVSMScene::idle()
{
    visit(IDLE);
}
void PVSMScene::reshape(int w,int h)
{
    getRenderer().reshape(w,h);
    getConsole().reshape(w,h);
    Info f;
    f._w=w;
    f._h=h;
    visit(RESHAPE,&f);
}
void PVSMScene::mouse(int button,int state,int x,int y)
{
    Info f;
    f._button=button;
    f._state=state;
    f._x=x;
    f._y=y;
    visit(MOUSE,&f);
}
void PVSMScene::motion(int x,int y)
{
    Info f;
    f._x=x;
    f._y=y;
    visit(MOTION,&f);
}
void PVSMScene::passiveMotion(int x,int y)
{
    Info f;
    f._x=x;
    f._y=y;
    visit(PASSIVE_MOTION,&f);
}
void PVSMScene::mouseWheel(int wheel,int dir,int x,int y)
{
    Info f;
    f._wheel=wheel;
    f._dir=dir;
    f._x=x;
    f._y=y;
    visit(MOUSE_WHEEL,&f);
}
void PVSMScene::keyboard(unsigned char key,int x,int y)
{
    Info f;
    f._key=key;
    f._x=x;
    f._y=y;
    visit(KEYBOARD,&f);
}
void PVSMScene::keyboardUp(unsigned char key,int x,int y)
{
    Info f;
    f._key=key;
    f._x=x;
    f._y=y;
    visit(KEYBOARD_UP,&f);
}
void PVSMScene::init()
{
    if(!_renderer)
        _renderer.reset(new DefaultRenderer);
    if(!_console)
        _console.reset(new DefaultConsole);
    getRenderer().init();
    visit(INIT_STAGE1);
    visit(INIT_STAGE2);
    visit(INIT_STAGE3);
    visit(INIT_STAGE4);
}
DefaultRenderer& PVSMScene::getRenderer()
{
    return *_renderer;
}
const DefaultRenderer& PVSMScene::getRenderer() const
{
    return *_renderer;
}
void PVSMScene::setRenderer(boost::shared_ptr<DefaultRenderer> rend)
{
    _renderer=rend;
}
DefaultConsole& PVSMScene::getConsole()
{
    return *_console;
}
const DefaultConsole& PVSMScene::getConsole() const
{
    return *_console;
}
int PVSMScene::getFPS() const
{
    return 30;
}
scalar PVSMScene::getFrameTime() const
{
    return 1.0f/getFPS();
}
void PVSMScene::frame()
{
    visit(FRAME);
    glutPostRedisplay();
}
sizeType PVSMScene::getStamp() const
{
    return _stamp;
}

//IO
void PVSMScene::readElement(const std::string& name,const ELEMENT& e)
{
    if(name == "<xmlattr>" || name == "Property")
        return;
    //check for source group
    if(name == "ProxyCollection" && PVSMAttribute<std::string>("name").read(e,"sources")) {
        for(ELEMENT::const_iterator beg=e.begin(),end=e.end(); beg!=end; beg++) {
            PVSMItem it;
            if(it.read(beg->first,beg->second))
                registerName(it._id._value,it._name._value);
        }
    } else if(name == "GlobalSetting") {
        _gSetting.read(name,e);
    } else {
        //check if this is a valid element
        for(sizeType i=0; i<(sizeType)_candidateObjs.size(); i++)
            if(_candidateObjs[i]->read(name,e)) {
                _objs[_candidateObjs[i]->_id._value]=_candidateObjs[i];
                _candidateObjs[i]=_candidateObjs[i]->copy();
            }
        //walk down
        for(ELEMENT::const_iterator beg=e.begin(),end=e.end(); beg!=end; beg++)
            readElement(beg->first,beg->second);
    }
}
void PVSMScene::setCandidates()
{
    _candidateObjs.clear();
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMTransform3));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMTransformFilter));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMMaterial));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMObjMesh));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMBox));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMPlane));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMSphere));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMLine));
    _candidateObjs.push_back(boost::shared_ptr<PVSMObject>(new PVSMPoint));
}
template <typename T>
boost::shared_ptr<T> PVSMScene::findObject(int id) const
{
    boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::const_iterator iter;
    if((iter=_objs.find(id)) == _objs.end())
        return boost::shared_ptr<T>((T*)NULL);
    return boost::dynamic_pointer_cast<T>(iter->second);
}
template <typename T>
boost::shared_ptr<T> PVSMScene::findOrCreateObject(int& id)
{
    boost::shared_ptr<T> ret=findObject<T>(id);
    if(!ret) {
        id=getNewId();
        ret.reset(new T);
        ret->_id._value=id;
        _objs[id]=ret;
    }
    return ret;
}
boost::shared_ptr<PVSMScene::PVSMObject> PVSMScene::findPtr(int id) const
{
    return _objs.find(id)->second;
}
int PVSMScene::getNewId() const
{
    int i=0;
    while(_objs.find(i) != _objs.end())
        i++;
    return i;
}
std::string PVSMScene::findName(int id) const
{
    boost::unordered_map<int,std::string>::const_iterator iter;
    if((iter=_names.find(id)) == _names.end())
        return "";
    return iter->second;
}
void PVSMScene::registerName(int id,const std::string& name)
{
    _names[id]=name;
}

//group management
void PVSMScene::visit(ACTION_TYPE type,const Info* f)
{
    _stamp++;
    for(boost::unordered_map<int,boost::shared_ptr<PVSMObject> >::iterator
            beg=_objs.begin(),end=_objs.end(); beg!=end; beg++) {
        boost::shared_ptr<PVSMSource> src=
            boost::dynamic_pointer_cast<PVSMSource>(beg->second);
        if(src) {
            vector<boost::shared_ptr<PVSMAction> >& actions=src->_actions;
            for(sizeType i=0; i<(sizeType)actions.size(); i++)
                performAction(src,actions[i].get(),type,f);
        }
    }
    for(boost::unordered_set<boost::shared_ptr<PVSMAction> >::iterator
            beg=_sceneActions.begin(),end=_sceneActions.end(); beg!=end; beg++) {
        boost::shared_ptr<PVSMSource> temp;
        performAction(temp,beg->get(),type,f);
    }

}
void PVSMScene::performAction(boost::shared_ptr<PVSMSource>& src,PVSMAction* action,ACTION_TYPE type,const Info* f)
{
    if(action->getStamp() == _stamp)return;
    switch(type) {
    case PRE_DRAW:
        if(!action->preDraw(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case DRAW:
        if(!action->draw(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case POST_DRAW:
        if(!action->postDraw(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case IDLE:
        if(!action->idle(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case FRAME:
        if(!action->frame(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case RESHAPE:
        if(!action->reshape(*this,src,f->_w,f->_h)) {
            action->setStamp(_stamp);
        }
        break;

    case INIT_STAGE1:
        if(!action->initStage1(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case INIT_STAGE2:
        if(!action->initStage2(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case INIT_STAGE3:
        if(!action->initStage3(*this,src)) {
            action->setStamp(_stamp);
        }
        break;
    case INIT_STAGE4:
        if(!action->initStage4(*this,src)) {
            action->setStamp(_stamp);
        }
        break;

    case MOUSE:
        if(action->isActive() && !action->mouse(*this,src,f->_button,f->_state,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    case MOTION:
        if(action->isActive() && !action->motion(*this,src,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    case PASSIVE_MOTION:
        if(action->isActive() && !action->passiveMotion(*this,src,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    case MOUSE_WHEEL:
        if(action->isActive() && !action->mouseWheel(*this,src,f->_wheel,f->_dir,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    case KEYBOARD:
        if(action->isActive() && !action->keyboard(*this,src,f->_key,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    case KEYBOARD_UP:
        if(action->isActive() && !action->keyboardUp(*this,src,f->_key,f->_x,f->_y)) {
            action->setStamp(_stamp);
        }
        break;
    }
}
void PVSMScene::drawObject()
{
    visit(DRAW);
}

PVSMSwitch::PVSMSwitch(boost::shared_ptr<PVSMAction> A,boost::shared_ptr<PVSMAction> B)
    :PVSMAction("Switch"),_A(A),_B(B),_switchA(true)
{
    _A->setActive(true);
    _B->setActive(false);
    _key='1';
}
bool PVSMSwitch::keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y)
{
    if(key==_key) {
        _switchA=!_switchA;
        _A->setActive(_switchA);
        _B->setActive(!_switchA);
    }
    return false;
}
void PVSMSwitch::setKey(unsigned char key)
{
    _key=key;
}