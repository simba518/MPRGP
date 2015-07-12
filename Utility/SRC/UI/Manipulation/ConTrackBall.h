#ifndef _CONTRACKBALL_H_
#define _CONTRACKBALL_H_

#include <boost/shared_ptr.hpp>
#include <QGLViewerExt.h>
#include <SelfRenderEle.h>
#include <AxisTorus.h>

namespace QGLVEXT{

  /**
   * @class ConTrackBall constrained track ball. the rotation will be
   * constrained to selected axis.
   * 
   */
  class ConTrackBall: public QObject, public MouseGrabber{

	Q_OBJECT

  public:
	ConTrackBall(pQGLViewerExt v, ManipulatedFrame* frame=NULL):viewer(v){
	  init();
	  addConstrainedFrame( NULL==frame ? viewer->camera()->frame():frame);
	  camera_translate[0] = 0.0f;
	  camera_translate[1] = 0.0f;
	  camera_translate[2] = 0.0f;
	}
	void addConstrainedFrame(ManipulatedFrame* frame){
	  assert(frame);
	  assert(constraint);
	  frame->setConstraint(constraint);
	}
  	void checkIfGrabsMouse(int x, int y, const Camera* const camera);

  	void setShow(const bool show){
	  m_show = show;
	  if (m_show){
		viewer->addSelfRenderEle(p_AxisTorus);
	  }else{
		viewer->removeSelfRenderEle(p_AxisTorus);
	  }
	}
  	void setEnable(const bool enable){
	  m_enabled = enable;
	}
	bool isShown()const{
	  return m_show;
	}
	bool isEnabled()const{
	  return m_enabled;
	}
	bool isHitted()const{
	  return m_hit;
	}

	void translate(const double x,const double y, const double z){
	  p_AxisTorus->translate(camera_translate[0]+x,
							 camera_translate[1]+y,
							 camera_translate[2]+z);
	}
	void scale(const double x,const double y, const double z){
	  p_AxisTorus->scale(x,y,z);
	}
	void translateCamera(const double x,const double y, const double z){
	  camera_translate[0] = x;
	  camera_translate[1] = y;
	  camera_translate[2] = z;
	}
	void draw()const{
	  p_AxisTorus->draw();
	}

  public slots:
  	void selectAxises(const vector<int> sel_group_ids);

  	void press(QMouseEvent* e);
  	void release(QMouseEvent* e);

  	void toggleShow(){
	  setShow((m_show ? false:true));
	}
	void toggleEnable(){
	  setEnable((m_enabled ? false:true));
	}

  signals:
	void selectObject(const int obj_id);

  protected:
	void init();

  private:
  	pQGLViewerExt viewer;
  	pAxisTorus p_AxisTorus;
  	AxisPlaneConstraint* constraint;

  	AxisPlaneConstraint::Type grabbed_type;
  	short constrained_axi; // 0: x, 1: y, 2: z, -1:none.
	bool m_enabled;
	bool m_show;
	bool m_hit;
	double camera_translate[3];
  };
  
  typedef boost::shared_ptr<ConTrackBall> pConTrackBall;
  
}//end of namespace

#endif /*_CONTRACKBALL_H_*/
