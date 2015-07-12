#ifndef PVSM_SCENE_H
#define PVSM_SCENE_H

#include "MathBasic.h"
#include "Renderer.h"
#include "Console.h"
#include "CollisionDetection.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

PRJ_BEGIN

class PVSMScene : public Renderable
{
public:
    typedef boost::property_tree::ptree ELEMENT;
    enum ACTION_TYPE {
        PRE_DRAW,
        DRAW,
        POST_DRAW,
        IDLE,
        FRAME,
        RESHAPE,

        MOUSE,
        MOTION,
        PASSIVE_MOTION,
        MOUSE_WHEEL,

        KEYBOARD,
        KEYBOARD_UP,

        INIT_STAGE1,
        INIT_STAGE2,
        INIT_STAGE3,
        INIT_STAGE4,
    };
    enum FILTER_TYPE {
        TRANSFORM,
    };
    enum SOURCE_TYPE {
        OBJMESH,
        BOX,
        PLANE,
        SPHERE,
        LINE,
        POINT,
    };
    //PVSM elements
    template <typename T>
    struct PVSMAttribute {
        PVSMAttribute(const std::string& name):_name(name) {}
        bool read(const ELEMENT& e) {
            std::string text=e.get<std::string>("<xmlattr>."+_name,"");
            if(text.empty())return false;
            istringstream iss(text);
            iss >> _value;
            return true;
        }
        bool read(const ELEMENT& e,const T& value) {
            return read(e) && _value == value;
        }
        void write(ELEMENT& e) const {
            ELEMENT attr;
            attr.put_value<T>(_value);
            e.add_child("<xmlattr>."+_name,attr);
        }
        std::string _name;
        T _value;
    };
    template <typename T,int comp>
    struct PVSMProperty {
        PVSMProperty(const std::string& name,const std::string& type="Element"):_name(name),_type(type) {}
        bool read(const ELEMENT& e,int id) {
            ostringstream oss;
            oss << id << "." << _name;
            for(ELEMENT::const_iterator beg=e.begin(),end=e.end(); beg!=end; beg++)
                if(PVSMAttribute<std::string>("id").read(beg->second,oss.str().c_str()))
                    return readElement(beg->second);
            return false;
        }
        bool readElement(const ELEMENT& e) {
            vector<bool> set(comp,false);
            for(ELEMENT::const_iterator beg=e.begin(),end=e.end(); beg!=end; beg++) {
                if(beg->first == _type) {
                    PVSMAttribute<int> index("index");
                    PVSMAttribute<T> value("value");
                    if(value.read(beg->second)) {
                        if(!index.read(beg->second))index._value=0;
                        _value[index._value]=value._value;
                        //already set, data corruption
                        if(set[index._value])return false;
                        set[index._value]=true;
                    }
                }
            }
            //check component consistency
            for(sizeType i=0; i<(sizeType)set.size(); i++)
                if(!set[i])return false;
            return true;
        }
        void write(ELEMENT& e,int id) const {
            ELEMENT prop;
            ostringstream oss;
            oss << id << "." << _name;
            PVSMAttribute<std::string> val("id");
            val._value=oss.str();
            val.write(prop);

            for(int i=0; i<comp; i++) {
                ELEMENT propE;
                PVSMAttribute<int> index("index");
                index._value=i;
                index.write(propE);

                PVSMAttribute<T> value("value");
                value._value=_value[i];
                value.write(propE);
                prop.add_child(_type,propE);
            }
            e.add_child("Property",prop);
        }
        template <typename T2>
        void set(const Eigen::Matrix<T2,comp,1>& m) {
            for(int i=0; i<comp; i++)_value[i]=m[i];
        }
        void set(const T& v) {
            BOOST_STATIC_ASSERT(comp == 1);
            _value[0]=v;
        }
        Eigen::Matrix<T,comp,1> get() const {
            return Eigen::Map<const Eigen::Matrix<T,comp,1> >(_value);
        }
        T operator[](int i) const {
            return _value[i];
        }
        T& operator[](int i) {
            return _value[i];
        }
        std::string _name;
        std::string _type;
        T _value[comp];
    };
    struct PVSMItem {
        PVSMItem():_id("id"),_name("name") {}
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return name == "Item" && _id.read(e) && _name.read(e);
        }
        virtual std::string write(ELEMENT& e) const {
            _id.write(e),_name.write(e);
            return "Item";
        }
        PVSMAttribute<int> _id;
        PVSMAttribute<std::string> _name;
    };
    struct PVSMObject {
        PVSMObject():_id("id"),_name("name") {}
        virtual ~PVSMObject() {}
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return _id.read(e);
        }
        virtual std::string write(ELEMENT& e) const {
            _id.write(e);
            return "";
        }
        virtual bool link(PVSMScene& scene) {
            _name._value=scene.findName(_id._value);
            return true;
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMObject);
        }
        PVSMAttribute<int> _id;
        PVSMAttribute<std::string> _name;
        boost::shared_ptr<PVSMObject> _parent;
    };
    template <typename T>
    struct PVSMRef : public PVSMProperty<int,1> {
        PVSMRef(const std::string& name):PVSMProperty(name,"Proxy") {
            _value[0]=-1;
        }
        virtual bool link(PVSMScene& scene,PVSMObject* parent) {
            _ref=scene.findObject<T>(_value[0]);
            if(!_ref)return false;
            _ref->_parent=scene.findPtr(parent->_id._value);
            return true;
        }
        virtual bool linkDef(PVSMScene& scene,PVSMObject* parent) {
            _ref=scene.findOrCreateObject<T>(_value[0]);
            _ref->_parent=scene.findPtr(parent->_id._value);
            return true;
        }
        boost::shared_ptr<T> _ref;
    };
    //filters
    struct PVSMFilter : public PVSMObject {
        PVSMFilter():_group("group"),_input("Input") {}
        virtual FILTER_TYPE getType() const =0;
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMObject::read(name,e) && name == "Proxy" && _group.read(e,"filters");
        }
        virtual std::string write(ELEMENT& e) const {
            _group.write(e);
            _input.write(e,_id._value);
            PVSMObject::write(e);
            return "Proxy";
        }
        virtual bool link(PVSMScene& scene) {
            return PVSMObject::link(scene) && _input.link(scene,this);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>((PVSMFilter*)NULL);
        }
        PVSMAttribute<std::string> _group;
        PVSMRef<PVSMObject> _input;
    };
    struct PVSMTransform3 : public PVSMObject {
        PVSMTransform3():_group("group"),_type("type"),_P("Position"),_R("Rotation"),_S("Scale") {}
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMObject::read(name,e) && name == "Proxy" &&
                   _group.read(e,"extended_sources") && _type.read(e,"Transform3") &&
                   _P.read(e,_id._value) && _R.read(e,_id._value) && _S.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _group.write(e);
            _type.write(e);
            _P.write(e,_id._value);
            _R.write(e,_id._value);
            _S.write(e,_id._value);
            PVSMObject::write(e);
            return "Proxy";
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMTransform3);
        }
        Mat4 get() const;
        Vec3 getScale() const;
        Mat3 getRotate() const;
        Vec3 getTranslate() const;
        static Mat3 getQuat(scalar rot,const Vec3& axis);
        PVSMAttribute<std::string> _group;
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _P,_R,_S;
    };
    struct PVSMTransformFilter : public PVSMFilter {
        PVSMTransformFilter():_type("type"),_transform("Transform") {}
        virtual FILTER_TYPE getType() const {
            return TRANSFORM;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMFilter::read(name,e) && _type.read(e,"TransformFilter") &&
                   _input.read(e,_id._value) && _transform.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _transform.write(e,_id._value);
            return PVSMFilter::write(e);
        }
        Mat4 get() const {
            return _transform._ref->get();
        }
        virtual bool link(PVSMScene& scene) {
            return PVSMFilter::link(scene) && _transform.link(scene,this);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMTransformFilter);
        }
        PVSMAttribute<std::string> _type;
        PVSMRef<PVSMTransform3> _transform;
    };
    //material
    template <int comp>
    struct PVSMMat : public PVSMProperty<GLfloat,comp> {
        PVSMMat(const std::string& name,GLenum ENUM):PVSMProperty<GLfloat, comp>(name),_ENUM(ENUM) {}
        void setOpenGL() const {
            glMaterialfv(GL_FRONT_AND_BACK,_ENUM,PVSMProperty<GLfloat, comp>::_value);
        }
        GLenum _ENUM;
    };
    struct PVSMPSize : public PVSMMat<1> {
        PVSMPSize();
        void setOpenGL() const;
    };
    struct PVSMLWidth : public PVSMMat<1> {
        PVSMLWidth();
        void setOpenGL() const;
    };
    template <int comp>
    struct PVSMLight : public PVSMProperty<GLfloat,comp> {
        PVSMLight(const std::string& name,GLenum ENUM):PVSMProperty<GLfloat, comp>(name),_ENUM(ENUM) {}
        void setOpenGL(GLuint i) const {
            glLightfv(GL_LIGHT0+i,_ENUM,PVSMProperty<GLfloat, comp>::_value);
        }
        GLenum _ENUM;
    };
    struct PVSMSource;
    struct PVSMMaterial : public PVSMObject {
        PVSMMaterial();
        void randomInit(const PVSMSource* src);
        virtual bool read(const std::string& name,const ELEMENT& e);
        virtual std::string write(ELEMENT& e) const;
        virtual boost::shared_ptr<PVSMObject> copy() const;
        void setOpenGLMaterial();
        void setOpenGLLight(PVSMScene& s);
        //data
        PVSMAttribute<std::string> _type;
        //material
        PVSMMat<4> _MAmb,_MDiff,_MSpec,_MEmi;
        PVSMMat<1> _MShi;
        PVSMLWidth _MLWidth;
        PVSMPSize _MPSize;
        //light
        PVSMLight<4> _LPos,_LDir;
        PVSMLight<4> _LAmb,_LDiff,_LSpec;
        PVSMLight<1> _LSpotE,_LCutOff,_LAC0,_LAC1,_LAC2;
        GLuint _lightId;
    };
    //sources
    class PVSMAction;
    struct PVSMSource : public PVSMObject {
        PVSMSource():_group("group"),_mat("Material") {}
        virtual SOURCE_TYPE getType() const =0;
        virtual bool read(const std::string& name,const ELEMENT& e) {
            bool ret=PVSMObject::read(name,e) && name == "Proxy" && _group.read(e,"sources");
            _mat.read(e,_id._value);
            return ret;
        }
        virtual std::string write(ELEMENT& e) const {
            _group.write(e);
            PVSMObject::write(e);
            _mat.write(e,_id._value);
            return "Proxy";
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>((PVSMObject*)NULL);
        }
        virtual bool link(PVSMScene& scene) {
            bool ret=PVSMObject::link(scene) && _mat.linkDef(scene,this);
            _mat._ref->randomInit(this);
            return ret;
        }
        PVSMAttribute<std::string> _group;
        PVSMRef<PVSMMaterial> _mat;
        std::vector<boost::shared_ptr<PVSMAction> > _actions;
    };
    struct PVSMObjMesh : public PVSMSource {
        PVSMObjMesh():_type("type"),_file("FileName") {}
        virtual SOURCE_TYPE getType() const {
            return OBJMESH;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"objreader") && _file.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _file.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMObjMesh);
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<std::string,1> _file;
    };
    struct PVSMBox : public PVSMSource {
        PVSMBox():_type("type"),_Ctr("Center"),_X("XLength"),_Y("YLength"),_Z("ZLength") {}
        virtual SOURCE_TYPE getType() const {
            return BOX;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"CubeSource") &&
                   _Ctr.read(e,_id._value) && _X.read(e,_id._value) &&
                   _Y.read(e,_id._value) && _Z.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _Ctr.write(e,_id._value);
            _X.write(e,_id._value);
            _Y.write(e,_id._value);
            _Z.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMBox);
        }
        virtual BBox<scalar> getBox() const {
            Vec3 ext(_X[0],_Y[0],_Z[0]);
            return BBox<scalar>(_Ctr.get()-ext*0.5f,_Ctr.get()+ext*0.5f);
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _Ctr;
        PVSMProperty<scalar,1> _X,_Y,_Z;
    };
    struct PVSMPlane : public PVSMSource {
        PVSMPlane():_type("type"),_Ctr("Origin"),_P1("Point1"),_P2("Point2") {}
        virtual SOURCE_TYPE getType() const {
            return PLANE;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"PlaneSource") &&
                   _Ctr.read(e,_id._value) && _P1.read(e,_id._value) && _P2.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _Ctr.write(e,_id._value);
            _P1.write(e,_id._value);
            _P2.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMPlane);
        }
        Plane getPlane() const {
            return Plane(_Ctr.get().cast<scalar>(),
                         _P1.get().cast<scalar>(),
                         _P2.get().cast<scalar>());
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _Ctr,_P1,_P2;
    };
    struct PVSMSphere : public PVSMSource {
        PVSMSphere():_type("type"),_Ctr("Center"),_Rad("Radius") {}
        virtual SOURCE_TYPE getType() const {
            return SPHERE;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"SphereSource") &&
                   _Ctr.read(e,_id._value) && _Rad.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _Ctr.write(e,_id._value);
            _Rad.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMSphere);
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _Ctr;
        PVSMProperty<scalar,1> _Rad;
    };
    struct PVSMLine : public PVSMSource {
        PVSMLine():_type("type"),_P1("Point1"),_P2("Point2") {}
        virtual SOURCE_TYPE getType() const {
            return LINE;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"LineSource") &&
                   _P1.read(e,_id._value) && _P2.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _P1.write(e,_id._value);
            _P2.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMLine);
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _P1,_P2;
    };
    struct PVSMPoint : public PVSMSource {
        PVSMPoint():_type("type"),_Ctr("Center") {}
        virtual SOURCE_TYPE getType() const {
            return POINT;
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return PVSMSource::read(name,e) && _type.read(e,"PointSource") && _Ctr.read(e,_id._value);
        }
        virtual std::string write(ELEMENT& e) const {
            _type.write(e);
            _Ctr.write(e,_id._value);
            return PVSMSource::write(e);
        }
        virtual boost::shared_ptr<PVSMObject> copy() const {
            return boost::shared_ptr<PVSMObject>(new PVSMPoint);
        }
        PVSMAttribute<std::string> _type;
        PVSMProperty<scalar,3> _Ctr;
    };
    //actions
    class PVSMAction
    {
    public:
        PVSMAction():_group(""),_stamp(-1),_active(true) {}
        PVSMAction(const std::string& group):_group(group),_stamp(-1),_active(true) {}
        virtual ~PVSMAction() {}
        virtual boost::shared_ptr<PVSMAction> copy(boost::shared_ptr<PVSMAction> self) {
            ASSERT(false);
            return boost::shared_ptr<PVSMAction>();
        }
        //display
        virtual bool preDraw(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool draw(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool postDraw(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool idle(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool frame(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool reshape(PVSMScene& s,boost::shared_ptr<PVSMSource> src,int w,int h) {
            return false;
        }
        //mouse
        virtual bool mouse(PVSMScene& s,boost::shared_ptr<PVSMSource> src,int button,int state,int x,int y) {
            return false;
        }
        virtual bool motion(PVSMScene& s,boost::shared_ptr<PVSMSource> src,int x,int y) {
            return false;
        }
        virtual bool passiveMotion(PVSMScene& s,boost::shared_ptr<PVSMSource> src,int x,int y) {
            return false;
        }
        virtual bool mouseWheel(PVSMScene& s,boost::shared_ptr<PVSMSource> src,int wheel,int dir,int x,int y) {
            return false;
        }
        //keyboard
        virtual bool keyboard(PVSMScene& s,boost::shared_ptr<PVSMSource> src,unsigned char key,int x,int y) {
            return false;
        }
        virtual bool keyboardUp(PVSMScene& s,boost::shared_ptr<PVSMSource> src,unsigned char key,int x,int y) {
            return false;
        }
        virtual void activate() {}
        virtual void deactivate() {}
        //initialize
        virtual bool initStage1(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool initStage2(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool initStage3(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        virtual bool initStage4(PVSMScene& s,boost::shared_ptr<PVSMSource> src) {
            return false;
        }
        //utility
        bool isInGroup(const PVSMSource& src,const std::string& forceGroup=std::string("")) const;
        void setStamp(sizeType stamp) {
            _stamp=stamp;
        }
        sizeType getStamp() const {
            return _stamp;
        }
        void setActive(bool active) {
            _active=active;
            if(_active)activate();
            else deactivate();
        }
        bool isActive() const {
            return _active;
        }
    protected:
        static void checkMode(int button,int state,int buttonMode,bool& mode);
        std::string _group;
        boost::unordered_map<std::string,SOURCE_TYPE> _names;
    private:
        sizeType _stamp;
        bool _active;
    };
    struct Info {
        int _w,_h;
        int _x,_y;
        int _wheel,_dir;
        int _button,_state;
        unsigned char _key;
    };
    //global setting
    struct GlobalSetting {
        GlobalSetting():_title("title"),_size("size"),_pos("pos"),_clearColor("clearColor"),_textColor("textColor") {
            _title._value="PVSMScene Application";
            _size.set(Eigen::Matrix<int,2,1>(800,600));
            _pos.set(Eigen::Matrix<int,2,1>(0,0));
            _clearColor.set(Vec4f(1.0f,1.0f,1.0f,1.0f));
            _textColor.set(Vec4f(1.0f,0.5f,0.5f,1.0f));
        }
        virtual bool read(const std::string& name,const ELEMENT& e) {
            return name == "GlobalSetting" &&
                   _title.read(e) && _size.read(e,-1) && _pos.read(e,-1),_clearColor.read(e,-1) && _textColor.read(e,-1);
        }
        virtual std::string write(ELEMENT& e) const {
            _title.write(e);
            _size.write(e,-1);
            _pos.write(e,-1);
            _clearColor.write(e,-1);
            _textColor.write(e,-1);
            return "GlobalSetting";
        }
        PVSMAttribute<std::string> _title;
        PVSMProperty<int,2> _size,_pos;
        PVSMProperty<GLfloat,4> _clearColor;
        PVSMProperty<GLfloat,4> _textColor;
    };
public:
    typedef boost::unordered_map<std::string,boost::shared_ptr<PVSMAction> > ACTION_MAP;
    //IO
    bool read(const std::string& path);
    bool write(const std::string& path) const;
    void clearUnused();
    GlobalSetting getGlobalSetting();
    //group management
    void delAction(boost::shared_ptr<PVSMAction> action);
    void addAction(boost::shared_ptr<PVSMAction> action);
    void delSceneAction(boost::shared_ptr<PVSMAction> action);
    void addSceneAction(boost::shared_ptr<PVSMAction> action);
    void updateScene();
    //events
    void render();
    void idle();
    void reshape(int w,int h);
    void mouse(int button,int state,int x,int y);
    void motion(int x,int y);
    void passiveMotion(int x,int y);
    void mouseWheel(int wheel,int dir,int x,int y);
    void keyboard(unsigned char key,int x,int y);
    void keyboardUp(unsigned char key,int x,int y);
    void init();
    //render event
    DefaultRenderer& getRenderer();
    const DefaultRenderer& getRenderer() const;
    void setRenderer(boost::shared_ptr<DefaultRenderer> rend);
    DefaultConsole& getConsole();
    const DefaultConsole& getConsole() const;
    int getFPS() const;
    scalar getFrameTime() const;
    void frame();
    sizeType getStamp() const;
    //utility
    template <typename T>
    static boost::shared_ptr<T> nextFilter(boost::shared_ptr<PVSMObject> obj) {
        while(obj->_parent) {
            boost::shared_ptr<T> ret=
                boost::dynamic_pointer_cast<T>(obj->_parent);
            if(ret)return ret;
            obj=obj->_parent;
        }
        return boost::shared_ptr<T>();
    }
    template <typename T2>
    static Vec3 transform(const Mat4& T,const Eigen::Matrix<T2,3,1>& v) {
        Vec4 homo;
        homo.block<3,1>(0,0)=v.template cast<scalar>();
        homo[3]=1.0f;
        homo=T*homo;
        return homo.block<3,1>(0,0)/homo[3];
    }
    template <typename T2>
    static Vec4 transformH(const Mat4& T,const Eigen::Matrix<T2,3,1>& v) {
        Vec4 homo;
        homo.block<3,1>(0,0)=v.template cast<scalar>();
        homo[3]=1.0f;
        return T*homo;
    }
    static Mat4 getTrans(boost::shared_ptr<PVSMScene::PVSMObject> obj) {
        Mat4 T=Mat4::Identity();
        boost::shared_ptr<PVSMScene::PVSMTransformFilter> filter=
            nextFilter<PVSMScene::PVSMTransformFilter>(obj);
        while(filter) {
            T=filter->get()*T;
            filter=nextFilter<PVSMScene::PVSMTransformFilter>(filter);
        }
        return T;
    }
protected:
    //IO
    void readElement(const std::string& name,const ELEMENT& e);
    void setCandidates();
    template <typename T>
    boost::shared_ptr<T> findObject(int id) const;
    template <typename T>
    boost::shared_ptr<T> findOrCreateObject(int& id);
    boost::shared_ptr<PVSMObject> findPtr(int id) const;
    int getNewId() const;
    std::string findName(int id) const;
    void registerName(int id,const std::string& name);
    //group management
    void visit(ACTION_TYPE type,const Info* f=NULL);
    void performAction(boost::shared_ptr<PVSMSource>& src,PVSMAction* action,ACTION_TYPE type,const Info* f);
    virtual void drawObject();
protected:
    //scene data
    std::vector<boost::shared_ptr<PVSMObject> > _candidateObjs;
    boost::unordered_map<int,boost::shared_ptr<PVSMObject> > _objs;
    boost::unordered_map<int,std::string> _names;
    //global setting
    GlobalSetting _gSetting;
    boost::shared_ptr<DefaultRenderer> _renderer;
    boost::shared_ptr<DefaultConsole> _console;
    sizeType _stamp;
    //scene action
    boost::unordered_set<boost::shared_ptr<PVSMAction> > _candidateActions;
    boost::unordered_set<boost::shared_ptr<PVSMAction> > _sceneActions;
};
class PVSMSwitch : public PVSMScene::PVSMAction
{
public:
    PVSMSwitch(boost::shared_ptr<PVSMAction> A,boost::shared_ptr<PVSMAction> B);
    virtual bool keyboard(PVSMScene& s,boost::shared_ptr<PVSMScene::PVSMSource> src,unsigned char key,int x,int y);
    void setKey(unsigned char key);
private:
    boost::shared_ptr<PVSMAction> _A,_B;
    unsigned char _key;
    bool _switchA;
};

PRJ_END

#endif
