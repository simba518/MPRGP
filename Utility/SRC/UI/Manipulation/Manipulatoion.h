#ifndef _MANIPULATOION_H_
#define _MANIPULATOION_H_

#include <QObject>
#include <QMouseEvent>
#include <QGLViewerExt.h>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <Log.h>
#include <ConTrackBall.h>
using namespace qglviewer;
using namespace Eigen;

namespace QGLVEXT{
  
  /**
   * @class LocalframeManipulatoion
   * 
   */
  class LocalframeManipulatoion:public QObject,public SelfRenderEle,public Selectable{

	Q_OBJECT
	
  public:
	LocalframeManipulatoion(pQGLViewerExt v);
	virtual void setShow(const bool show);
	bool isShown(){
	  return con_track_ball->isShown();
	}
	virtual void draw()const{}
	
	// select
	virtual int totalEleNum()const = 0;
	virtual void drawWithNames()const = 0;
	virtual void select(const vector<int> &sel_ids) = 0;

	// manipulate
	virtual void computeCurrentPosition(Vector3d& pos_xyz)const = 0;
	const Vector3d &currentPosition()const{
	  return current_pos;
	}
	virtual void applyTransform() = 0;
	void translateCamera(const double x,const double y, const double z){
	  con_track_ball->translateCamera(x,y,z);
	}

  public slots:
	void manipulate();
	void selectConTrackball(const int sel_axis);
	
  protected:
	bool enabled;
	bool shown;
	pQGLViewerExt viewer;
	ManipulatedFrame* frame;
	pConTrackBall con_track_ball;
	Vector3d current_pos;
  };
  typedef boost::shared_ptr<LocalframeManipulatoion> pLocalframeManipulatoion;

  class LocalframeManipulatoionExt:public LocalframeManipulatoion{

  public:
	LocalframeManipulatoionExt(pQGLViewerExt v):LocalframeManipulatoion(v){}
	virtual void setShow(const bool show);

	// select
	virtual int totalEleNum()const = 0;
	virtual void drawWithNames()const = 0;
	virtual void select(const vector<int> &sel_ids) = 0;
	
	// manipulate
	virtual void computeCurrentPosition(Vector3d &current_pos)const;
	virtual void applyTransform();

	virtual Matrix<double,3,-1> &getCurrentPositions() = 0;
	virtual const Matrix<double,3,-1> &getCurrentPositions()const = 0;

  protected:
	Matrix<double,3,-1> initial_pos;
  };

  class LocalframeManipulatoionCtrl:public QObject{
  
  	Q_OBJECT

  public:
  	LocalframeManipulatoionCtrl(pQGLViewerExt viewer,pLocalframeManipulatoion manipulate):
  	  viewer(viewer), manipulate(manipulate){

	  enabled = true;
  	  mouse_button = Qt::LeftButton;
  	  modify_key = Qt::ControlModifier;
  	  viewer->addSelfRenderEle(manipulate);
  	  connect(viewer,SIGNAL(mousePressSignal(QMouseEvent*)),this,SLOT(press(QMouseEvent*)));
  	  connect(viewer,SIGNAL(selectedIds(const vector<int>)),this,SLOT(select(const vector<int>)));
  	}
	void setEnable(const bool enable){
	  enabled = enable;
	}

  public slots:
  	void press (QMouseEvent *e){
	  if (enabled){
		if ((e->button()==mouse_button)&&(e->modifiers() == modify_key)){
		  manipulate->prepareSelection();
		  viewer->setSelector(manipulate);
		  viewer->select(QRect(e->pos().x(),e->pos().y(),10,10));
		}
	  }
  	}
  	void select(const vector<int> ids){
	  if (enabled){
		if (viewer->getSelector().get() == manipulate.get()){
		  manipulate->select(ids);
		  manipulate->setShow(ids.size() > 0);
		}
	  }
  	}
	void toggleEnable(){
	  enabled = enabled ? false:true;
	}

  private:
  	pQGLViewerExt viewer;
  	pLocalframeManipulatoion manipulate;
  	Qt::MouseButton mouse_button;
  	Qt::KeyboardModifiers modify_key;
	bool enabled;
  };
  typedef boost::shared_ptr<LocalframeManipulatoionCtrl> pLocalframeManipulatoionCtrl;  

}//end of namespace

#endif /*_MANIPULATOION_H_*/
