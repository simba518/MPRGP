#ifndef _DRAGCTRL_H_
#define _DRAGCTRL_H_

#include <QObject>
#include <QGLViewerExt.h>
using namespace qglviewer;

namespace QGLVEXT{

  /**
   * @class DragHook return the position of the point that is dragged, which is
   * used to set the z-deepth for the drag operation.
   * 
   */
  class DragHook{
	
  public:
	virtual void getDragedPoint(double point[3])const = 0;

  };
  typedef boost::shared_ptr<DragHook> pDragHook;

  class DragObserver{
	
  public: 
	virtual void selectDragEle(int id){}
	virtual void startDragScreen (int screen_x,int screen_y){}
	virtual void startDrag (double x,double y,double z){}
	virtual void dragTo (double x,double y,double z){}
	virtual void stopDrag (double x,double y,double z){}
  };
  typedef boost::shared_ptr<DragObserver> pDragObserver;

  /**
   * @class DragCtrl control the drag operation.
   */
  class DragCtrl:public QObject{

	Q_OBJECT
	
  public:
	DragCtrl(pQGLViewerExt viewer);
	DragCtrl(pQGLViewerExt viewer, pSelectable selector);
	void setKeyMouse(Qt::KeyboardModifiers key,Qt::MouseButton m_button){
	  this->modify_key = key;
	  this->mouse_button = m_button;
	}
	bool beginDrag()const{
	  return begin_drag;
	}
	void setDragHook(pDragHook drag_hook){
	  this->drag_hook = drag_hook;
	}
	void setObserver(pDragObserver observer){
	  _observer = observer;
	}
	void initZ_Deepth();

  public slots:
	bool press (QMouseEvent *e);
	bool release (QMouseEvent *e);
	bool move (QMouseEvent *e);
	void selectDragEle(const vector<int> ids);
	
  signals:
	void selectDragEle(int id);
	void startDrag (int screen_x,int screen_y);
	void startDrag (double x,double y,double z);
	void dragTo (double x,double y,double z);
	void stopDrag (double x,double y,double z);

  protected:
	void createConnections();
	
  private:
	pQGLViewerExt viewer;
	pDragHook drag_hook;
	pDragObserver _observer;
	pSelectable _selector;

	/**
	 * the mouse button and modify key that trags the draging operation, could
	 * be set by calling setKeyMouse(), and the default value is left+Ctrl.
	 */
	Qt::MouseButton mouse_button;
	Qt::KeyboardModifiers modify_key;

	double z_deepth;
	bool begin_drag;
	bool is_dragging;
	bool _selected;
  };

  typedef boost::shared_ptr< DragCtrl > pDragCtrl;
  
}//end of namespace

#endif /* _DRAGCTRL_H_ */
