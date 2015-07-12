#ifndef _QGLVIEWEREXT_H_
#define _QGLVIEWEREXT_H_

#include <QGLViewerHeaders.h>
#include <stack>
#include <QList>
#include <SelfRenderEle.h>
#include <TextForRender.h>
#include <Selectable.h>
#include <fstream>
using namespace qglviewer;
using namespace std;

namespace QGLVEXT{

  class QGLViewerExt:public QGLViewer{

	Q_OBJECT

	public:
	QGLViewerExt (QWidget *parent);

	// render
	void addSelfRenderEle (pSelfRenderEle ele);
	void removeSelfRenderEle (pSelfRenderEle ele);
	void addTextForRender (pTextForRender ele);
	void removeTextForRender (pTextForRender ele);
	bool toggleRemoveAddSelfRenderEle(pSelfRenderEle ele);//true if added
	bool contains(pSelfRenderEle ele)const{return self_render_ele.contains(ele);}

	// selection
	void setSelector(pSelectable selector){ this->selector = selector;}
	pSelectable_const getSelector()const{ return selector;}

	// draging
	Vec getWorldCoords(int x,int y,double z)const;
	Vec getWorldCoords(const QPoint &screen_xy,double z)const;
	Vec getScreenCoords(double world_x,double world_y,double world_z)const;

	// other
	const QPoint &getMousePos()const{
	  return mouse_pos;
	}

  public slots:
	// rendering
	void update (){  
	  QGLViewer::updateGL();
	}
	void resetSceneBoundBox(double x0,double y0,double z0,
							double x1,double y1,double z1);
	void toggleDrawLights(){
	  // draw the lights's positions, for debuging.
	  draw_lights = draw_lights ? false : true;
	  this->update();
	}
	void show3DGrid();
	void drawBackground()const;

	// selection
	void select (QRect select_rect);
	void select	(const QMouseEvent *event);
 
	// animation
	void startAnimation (){ QGLViewer::startAnimation(); emit resetAnimation();	}
	void stopAnimation (){ QGLViewer::stopAnimation();  emit resetAnimation();	}
	void pauseAnimation (){
	  if (animationIsStarted()){
		QGLViewer::stopAnimation();
		update();
	  }else{
		QGLViewer::startAnimation();
	  }
	}

	void backward(){ QGLViewer::stopAnimation();  emit backwardAnimation(); }
	void forward(){  QGLViewer::stopAnimation(); emit forwardAnimation();   }
	void fastForward(){ QGLViewer::stopAnimation();  emit ffAnimation(); 	}
	void fastBackward(){QGLViewer::stopAnimation();  emit fbAnimation();  	}
	void stepByStep (){step_by_step = step_by_step ? false:true;}
	void animate (){
	  emit updateAnimation();
	  if(step_by_step){
		QGLViewer::stopAnimation();
	  }
	}

	// Mouse events
	virtual void mousePressEvent (QMouseEvent *e);
	virtual void mouseMoveEvent (QMouseEvent *e);
	virtual void mouseReleaseEvent (QMouseEvent *e);
	virtual void wheelEvent (QWheelEvent *e);

	// keypress events
	virtual void keyPressEvent (QKeyEvent *e);

	// push and pop status
	void pushStatus();
	void restoreStatus();
	void popStatus();

	// init file
	void loadStateFile();
	void saveStateFile();

  signals:
	void resetAnimation();
	void updateAnimation();
	void ffAnimation(); // fast forward
	void fbAnimation(); // fast backward
	void forwardAnimation();	
	void backwardAnimation();
	void selectedIds(const vector<int> sel_ids);

	void mousePressSignal (QMouseEvent *e);
	void mouseMoveSignal (QMouseEvent *e);
	void mouseReleaseSignal (QMouseEvent *e);
	void keyPressSignal (QKeyEvent *e);
	void wheelSignal (QWheelEvent *e);

  protected:
	void init();

	// rendering
	virtual void draw ();
	virtual void selfRender ();
	virtual void displayText();
	virtual void drawMouse();
	virtual void draw3DGrid()const;

	// other
	virtual QString helpString ()const;

	// selection
	void drawWithNames(){ if (selector != NULL)	selector->drawWithNames(); }
	void endSelection(const QPoint&p);

  protected:
	QList<pSelfRenderEle > self_render_ele;
	QList<pTextForRender> _textForRender;
	pSelectable selector;
	bool step_by_step; // control animation
	QPoint mouse_pos; // mouse position: (screen_x, screen_y).
	Qt::MouseButton buttonPressed;
	bool draw_lights; // draw the lights's positions, for debuging.
	std::stack<GLdouble *> MV_status; // record the modelview matrix.

	//add by qnn
	ofstream outfile;
	bool m_show3DGrid;
  };

  // typedef boost::shared_ptr< QGLViewerExt > pQGLViewerExt; 
  typedef QGLViewerExt * pQGLViewerExt;
  typedef const QGLViewerExt * pQGLViewerExt_const;

}//end of namespace

#endif /* _QGLVIEWEREXT_H_ */
