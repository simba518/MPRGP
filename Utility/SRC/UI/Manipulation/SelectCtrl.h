#ifndef _SELECTCTRL_H_
#define _SELECTCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <QGLViewerExt.h>
#include <SelfRenderRect.h>
#include <Selectable.h>

namespace QGLVEXT{

  class SelectObserver{
	
  public: 
	virtual void addSelection(const vector<int> &sel_ids) = 0;
	virtual void removeSelection(const vector<int> &sel_ids) = 0;
  };
  typedef boost::shared_ptr<SelectObserver> pSelectObserver;

  enum SEL_STATUS{ADD_ELE,REMOVE_ELE,DO_NOTHING_ON_SEL};
  
  /**
   * @class SelectCtrl manager the basic select operation process.
   * 
   */
  class SelectCtrl:public QObject{
	
	Q_OBJECT
	
  public:
	SelectCtrl(pQGLViewerExt viewer,pSelectable selector,const int minSelRect=0);
	void setAddKeyMouse(Qt::KeyboardModifiers key,Qt::MouseButton button){
	  add_modify_key = key;
	  add_mouse_button = button;
	}
	void setRmKeyMouse(Qt::KeyboardModifiers key,Qt::MouseButton button){
	  rm_mouse_button = button;
	  rm_modify_key = key;
	}
	void setSelector( pSelectable selector){
	  this->selector = selector;
	}
	void setObserver( pSelectObserver observer){
	  _observer = observer;
	}
	bool isActive()const{ return enable_op;	}

  public slots:
	bool press (QMouseEvent *e);
	bool release (QMouseEvent *e);
	bool move (QMouseEvent *e);

	void endSelection(const vector<int> sel_ids);
	void togglePrintSelEle(){
	  print_selected_nodes = (print_selected_nodes ? false:true);
	}

	void enableSelOp(){ enable_op = true; }
	void disableSelOp(){ enable_op = false; }
	bool toggleSelOp(){
	  enable_op = enable_op ? false:true;
	  return enable_op;
	}

  signals:
	void addSelEleMsg(const vector<int> sel_ids); // add selected elements
	void removeSelEleMsg(const vector<int> sel_ids); // remove selected elements
	void rectChanged()const;
  	void select(QRect select_rect)const;

  protected:
	void createConnections();
	void printSelection(const vector<int> &selIds, bool print=true)const;
	
  protected:
	pQGLViewerExt viewer;
	pSelectable selector;
	pSelectObserver _observer;
	SEL_STATUS sel_status;
	bool enable_op;
	bool print_selected_nodes;

	pSelfRenderRect select_rect;
	bool begin_select;
	int minimalSelRect; // the minimal width and height for the select rect,
						// even the mouse is not move.

	/**
	 * the mouse button and modify key that trigs the selection operation
	 * (including add and remove operation),the default value is "shift+left
	 * button" for adding, and "shift+right button" for removing.
	 */
	Qt::MouseButton add_mouse_button;
	Qt::KeyboardModifiers add_modify_key;
	Qt::MouseButton rm_mouse_button;
	Qt::KeyboardModifiers rm_modify_key;
  };
  
  typedef boost::shared_ptr<SelectCtrl> pSelectCtrl;
  
}//end of namespace

#endif /* _SELECTCTRL_H_ */
