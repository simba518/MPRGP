#include <iostream>
using namespace std;

#include <QMouseEvent>
#include <boost/foreach.hpp>
#include <QFileDialog>
#include <Log.h>
#include "QGLViewerExt.h"

using namespace QGLVEXT;

QGLViewerExt::QGLViewerExt (QWidget *parent):QGLViewer(parent){

  step_by_step = false;
  draw_lights = false;
  buttonPressed = Qt::NoButton;
  m_show3DGrid = true;
}

void QGLViewerExt::addSelfRenderEle(pSelfRenderEle ele){

  if(!self_render_ele.contains(ele)){
	if (ele){
	  self_render_ele.push_back(ele);
	}else{
	  ERROR_LOG("the input pSelfRenderEle is null!");
	}
  }
  update();
}

void QGLViewerExt::removeSelfRenderEle (pSelfRenderEle ele){

  self_render_ele.removeAll(ele);
  update();
}

void QGLViewerExt::addTextForRender (pTextForRender ele){
  if(!_textForRender.contains(ele)){
	if (ele){
	  _textForRender.push_back(ele);
	}else{
	  ERROR_LOG("the input pTextForRender is null!");
	}
  }
  update();
}

void QGLViewerExt::removeTextForRender (pTextForRender ele){
  _textForRender.removeAll(ele);
  update();
}

bool QGLViewerExt::toggleRemoveAddSelfRenderEle (pSelfRenderEle ele){

  bool added = true;
  if(self_render_ele.contains(ele)){
	removeSelfRenderEle (ele);
	added = false;
  }else if (ele != NULL){
	addSelfRenderEle (ele);
  }else{
	ERROR_LOG("the input pSelfRenderEle is null.");
  }
  return added;
}

void QGLViewerExt::selfRender(){

  // render the elements with first order
  BOOST_FOREACH(pSelfRenderEle &ele,self_render_ele){
	if(ele->getRenderPriority() == FIRST_RENDER){
	  ele->render();
	}
  }

  // render the elements with second order
  BOOST_FOREACH(pSelfRenderEle &ele,self_render_ele){
	if(ele->getRenderPriority() == SECOND_RENDER){
	  ele->render();
	}
  }

  // render the elements with third order
  BOOST_FOREACH(pSelfRenderEle &ele,self_render_ele){
	if(ele->getRenderPriority() == THIRD_RENDER){
	  ele->render();
	}
  }

  // render the elements with fourth order
  BOOST_FOREACH(pSelfRenderEle &ele,self_render_ele){
	if(ele->getRenderPriority() == FOURTH_RENDER){
	  ele->render();
	}
  }

}

void QGLViewerExt::draw(){

  drawMouse();
	
  // drawBackground();
  glPushMatrix();
  selfRender();
  displayText();
  glPopMatrix();	
	
  if (m_show3DGrid){
	draw3DGrid();	
  }

  // debug lights
  if (draw_lights){
  	drawLight(GL_LIGHT0);
  	drawLight(GL_LIGHT1);
  	drawLight(GL_LIGHT2);
  	drawLight(GL_LIGHT3);
  }
}

void QGLViewerExt::displayText(){

  qglColor(QColor(0,0,0));
  glDisable(GL_LIGHTING);
  BOOST_FOREACH(pTextForRender textSet, _textForRender){
	BOOST_FOREACH(const TextWithPosition &t, *textSet){
	  drawText(t.getX(),t.getY(),t.getContent().c_str(),t.getFont());
	}
  }
  glEnable(GL_LIGHTING);
}

void QGLViewerExt::mousePressEvent (QMouseEvent *e){

  buttonPressed = e->button();
  this->mouse_pos = e->pos();
  QGLViewer::mousePressEvent(e);
  emit mousePressSignal(e);
}

void QGLViewerExt::mouseMoveEvent (QMouseEvent *e){

  this->mouse_pos = e->pos();
  QGLViewer::mouseMoveEvent(e);
  emit mouseMoveSignal(e);
  if (buttonPressed == Qt::NoButton){
	updateGL();
  }
}

void QGLViewerExt::mouseReleaseEvent (QMouseEvent *e){

  buttonPressed = Qt::NoButton;
  this->mouse_pos = e->pos();
  QGLViewer::mouseReleaseEvent(e);
  emit mouseReleaseSignal(e);
}

void QGLViewerExt::wheelEvent (QWheelEvent *e){

  QGLViewer::wheelEvent(e);
  emit wheelSignal(e);
}

void QGLViewerExt::keyPressEvent (QKeyEvent *e){

  QGLViewer::keyPressEvent(e);
  emit keyPressSignal (e);
}

QString QGLViewerExt::helpString() const{

  QString text("<h2>S i m p l e V i e w e r</h2>");
  return text;
}

void QGLViewerExt::resetSceneBoundBox(double x0,double y0,double z0,
									  double x1,double y1,double z1){
	
  Vec v_min(x0,y0,z0);
  Vec v_max(x1,y1,z1);

  const double scence_radius = (v_max-v_min).norm()/2.0f;
  if( scence_radius > 0.0f ){

	// load identity
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	GLdouble mvm[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
	QGLViewer::camera()->setFromModelViewMatrix(mvm);

	GLdouble prj[16];
	glGetDoublev(GL_PROJECTION_MATRIX, prj);
	QGLViewer::camera()->setFromModelViewMatrix(prj);

	// scale scene
	QGLViewer::setSceneCenter((v_min+v_max)/2.0f);
	QGLViewer::setSceneRadius(scence_radius);
	QGLViewer::showEntireScene();

	// set lights
	const Vec vd = camera()->viewDirection();

	float pos[4] = {0.0, 0.0, 0.0, 0.0};
	// pos[0] = -vd[0]*scence_radius; 	
	// pos[1] = -vd[1]*scence_radius; 
	// pos[2] = -vd[2]*scence_radius;
	pos[0] = -2; 	
	pos[1] = 2; 
	pos[2] = 2;
	glLightfv(GL_LIGHT0, GL_POSITION, pos);

	pos[0] = 2; 	
	pos[1] = 2; 
	pos[2] = 2;
	glLightfv(GL_LIGHT1, GL_POSITION, pos);
  }
}

void QGLViewerExt::init(){

  // // set lighting
  const double scalor = 0.5f;
  const GLfloat light_ambient[4]  = {0.8f*scalor, 0.8f*scalor, 0.8f*scalor, 1.0f*scalor};
  const GLfloat light_diffuse[4]  = {0.8f*scalor, 0.8f*scalor, 0.8f*scalor, 1.0f*scalor};
  const GLfloat light_specular[4] = {1.0f*scalor, 1.0f*scalor, 1.0f*scalor, 1.0f*scalor};

  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glEnable(GL_LIGHT0);

  glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);
  glEnable(GL_LIGHT1);

  glEnable(GL_LIGHTING);

  // other
  restoreStateFromFile();
  // setBackgroundColor(QColor(255,255,255));
  glBlendFunc(GL_ONE, GL_ONE);

  resetSceneBoundBox(-3,-3,-3,3,3,3);

  // Make camera the default manipulated frame.
  setManipulatedFrame( camera()->frame() );

  setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltModifier);
  setHandlerKeyboardModifiers(QGLViewer::FRAME, Qt::NoModifier);
  // setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlModifier);

  setMouseTracking(true);

  // // set mouse binding
  // setMouseBinding(Qt::MiddleButton, CAMERA, NO_MOUSE_ACTION);
  // setMouseBinding(Qt::LeftButton, CAMERA, NO_MOUSE_ACTION);
  // setMouseBinding(Qt::RightButton, CAMERA, NO_MOUSE_ACTION);
  // setWheelBinding(Qt::NoModifier, CAMERA, NO_MOUSE_ACTION);
   
  // setMouseBinding(Qt::MiddleButton+Qt::CTRL, CAMERA, ZOOM);
  // setMouseBinding(Qt::LeftButton+Qt::CTRL, CAMERA, ROTATE);
  // setMouseBinding(Qt::RightButton+Qt::CTRL, CAMERA, TRANSLATE);
  // setWheelBinding(Qt::ControlModifier, CAMERA, MOVE_FORWARD);
}

void QGLViewerExt::select(const QMouseEvent *event){
  /// do nothing, disable the default behavior.
}

void QGLViewerExt::select (QRect select_rect){

  if(selector != NULL){

	// set buffer size
	const int sel_num = selector->totalEleNum();
	setSelectBufferSize (sel_num*4+4);

	// perform selection
	select_rect.normalized();
	setSelectRegionWidth(select_rect.width());
	setSelectRegionHeight(select_rect.height());
	QGLViewer::select(select_rect.center());
  }
}

void QGLViewerExt::endSelection(const QPoint&p){

  const GLuint* select_buffer = selectBuffer();
  if (select_buffer != NULL && selector != NULL){

	// get the selected ids
	glFlush();
	const GLint nbHits = glRenderMode(GL_RENDER);
	vector<int> sel_ids;
	for (int i=0; i<nbHits; ++i){
	  /// @todo
	  // ASSERTIN ( i,0,(selector->totalEleNum())*4 );
	  // ASSERTIN ( (int)(select_buffer)[4*i+3],0, selector->totalEleNum() );
	  sel_ids.push_back((select_buffer)[4*i+3]);
	}

	// send the seleted ids
	emit selectedIds(sel_ids);

	// update the scence
	update();
  }
}

Vec QGLViewerExt::getWorldCoords(int x,int y,double z)const{

  const Vec src(x,y,z);
  const Vec mouse_pos = this->camera()->unprojectedCoordinatesOf(src);
  return mouse_pos;
}

Vec QGLViewerExt::getWorldCoords(const QPoint &screen_xy,double z)const{

  return getWorldCoords(screen_xy.x(),screen_xy.y(),z);
}

Vec QGLViewerExt::getScreenCoords(double w_x,double w_y,double w_z)const{

  const Vec src(w_x,w_y,w_z);
  const Vec screen_pos = this->camera()->projectedCoordinatesOf(src);
  return screen_pos;
}

void QGLViewerExt::pushStatus(){

  GLdouble *m = new GLdouble[16];
  this->camera()->getModelViewMatrix(m);
  MV_status.push(m);
}

void QGLViewerExt::restoreStatus(){

  if ( !MV_status.empty() ) {

	GLdouble *m = MV_status.top();
	QGLViewer::camera()->setFromModelViewMatrix(m);
	this->update();
  }
}

void QGLViewerExt::popStatus(){

  if ( !MV_status.empty() ) {

	GLdouble *m = MV_status.top();
	QGLViewer::camera()->setFromModelViewMatrix(m);
	this->update();
	MV_status.pop();
  }
}

void QGLViewerExt::loadStateFile(){

  QFileDialog  fileDialog;
  const QString fname = fileDialog.getOpenFileName(this);
  if (fname.size() > 0){
	this->setStateFileName( fname );
	this->restoreStateFromFile();
  }
}

void QGLViewerExt::saveStateFile(){

  QFileDialog  fileDialog;
  const QString fname = fileDialog.getSaveFileName(this);
  if (fname.size() > 0){
	this->setStateFileName( fname );
	this->saveStateToFile();
  }
}

void QGLViewerExt::drawMouse(){
  
  ///@todo the position of the mouse is not correct.
  int x = mouse_pos.x(), y = mouse_pos.y();
  startScreenCoordinatesSystem();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glColor3f(1.0f,1.0f,1.0f);
  glLineWidth(1.0);
  glBegin(GL_POLYGON);
  glVertex2f(x, y);
  glVertex2f(x, y+ 15);
  glVertex2f(x+3, y+11);
  glVertex2f(x+6, y+18);
  glVertex2f(x+8, y+18);
  glVertex2f(x+8, y+17);
  glVertex2f(x+5, y+10);
  glVertex2f(x+10, y+10);
  glEnd();
  glColor3f(0.0f,0.0f,0.0f);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glBegin(GL_LINE_LOOP);
  glVertex2f(x, y);
  glVertex2f(x, y+ 15);
  glVertex2f(x+3, y+11);
  glVertex2f(x+6, y+18);
  glVertex2f(x+8, y+18);
  glVertex2f(x+8, y+17);
  glVertex2f(x+5, y+10);
  glVertex2f(x+10, y+10);
  glEnd();
  glPopAttrib();
  stopScreenCoordinatesSystem();
}

void QGLViewerExt::show3DGrid(){
  
  m_show3DGrid = !m_show3DGrid;
}

void QGLViewerExt::draw3DGrid()const{

  glLineWidth(1.0f);
  glColor3f(1.0f,1.0f,1.0f);
  glRotated(90.0f,1.0f,0.0f,0.0f);
  drawGrid(QGLViewer::sceneRadius());
  glRotated(-90.0f,1.0f,0.0f,0.0f);

  glRotated(90.0f,0.0f,1.0f,0.0f);
  drawGrid(QGLViewer::sceneRadius());
  glRotated(-90.0f,0.0f,1.0f,0.0f);
	
}

void QGLViewerExt::drawBackground()const{

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  startScreenCoordinatesSystem();
  GLint viewport[4];
  camera()->getViewport(viewport);
  const GLint vx=viewport[0];
  const GLint vy=viewport[1];
  const GLint vw=viewport[2];
  const GLint vh=-viewport[3];

  glBegin(GL_QUADS);
  const double c = 0.22;
  glColor3f(c,c,c);
  glVertex2f(0, 0);
  glColor3f(c,c,c);
  glVertex2f(vx+vw,0);
  glColor3f(1.0,1.0,1.0); 
  glVertex2f(vx+vw,vy+vh);
  glColor3f(1.0,1.0,1.0); 
  glVertex2f(0,vy+vh);
  glEnd();
  stopScreenCoordinatesSystem();

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glPopAttrib();

}
