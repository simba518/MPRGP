#include <QMouseEvent>
#include <iostream>
using namespace std;

#include "DragCtrl.h"
using namespace QGLVEXT;

DragCtrl::DragCtrl(pQGLViewerExt _viewer):viewer(_viewer){

  begin_drag = false;
  is_dragging = false;
  _selected = false;
  mouse_button = Qt::LeftButton;
  modify_key = Qt::NoModifier;
  z_deepth = 0.5f;
  createConnections();
}

DragCtrl::DragCtrl(pQGLViewerExt _viewer, pSelectable selector):
  viewer(_viewer),_selector(selector){
  
  begin_drag = false;
  is_dragging = false;
  _selected = false;
  mouse_button = Qt::LeftButton;
  modify_key = Qt::NoModifier;
  z_deepth = 0.5f;
  createConnections();
}

void DragCtrl::createConnections(){

  connect( viewer, SIGNAL(mousePressSignal(QMouseEvent *)),this, SLOT(press(QMouseEvent *)) );
  connect( viewer, SIGNAL(mouseReleaseSignal(QMouseEvent *)),this, SLOT(release(QMouseEvent *)) );
  connect( viewer, SIGNAL(mouseMoveSignal(QMouseEvent *)),this, SLOT(move(QMouseEvent *)) );
  connect( viewer, SIGNAL(selectedIds(const vector<int> )),this, SLOT(selectDragEle(const vector<int> )) );
}

bool DragCtrl::press (QMouseEvent *e){

  begin_drag = false;
  is_dragging = false;
  if(viewer != NULL){
	if ((e->button()==mouse_button)&&(e->modifiers() == modify_key)){

	  begin_drag = true;
	  if(_observer)
		_observer->startDragScreen(e->pos().x(),e->pos().y());
	  emit startDrag (e->pos().x(),e->pos().y());

	  if(_selector){
		_selected = true;
		_selector->prepareSelection();
		viewer->setSelector(_selector);
		viewer->select(QRect(e->pos().x(),e->pos().y(),10,10));
	  }
	}	
  }
  return begin_drag;
}

bool DragCtrl::release (QMouseEvent *e){

  const bool perfromed = begin_drag;
  if (begin_drag || is_dragging){
	Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	if(_observer)
	  _observer->stopDrag(world_pos[0],world_pos[1],world_pos[2]);
	emit stopDrag(world_pos[0],world_pos[1],world_pos[2]);
  }
  begin_drag = false;
  is_dragging = false;
  return perfromed;
}

bool DragCtrl::move (QMouseEvent *e){

  if (begin_drag) {

	initZ_Deepth();
	is_dragging = true;
	begin_drag = false;
	const Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	if(_observer)
	  _observer->startDrag(world_pos[0],world_pos[1],world_pos[2]);
	emit startDrag(world_pos[0],world_pos[1],world_pos[2]);
  }else if (is_dragging){

	const Vec world_pos = viewer->getWorldCoords(e->pos(),z_deepth);
	if(_observer)
	  _observer->dragTo(world_pos[0],world_pos[1],world_pos[2]);
	emit dragTo(world_pos[0],world_pos[1],world_pos[2]);
  }
  return is_dragging;
}

void DragCtrl::initZ_Deepth(){

  if (drag_hook && viewer){
	double p[3];
	drag_hook->getDragedPoint(p);
	const Vec screen_pos = viewer->getScreenCoords(p[0],p[1],p[2]);
	z_deepth = screen_pos[2];
  }else{
	z_deepth = 0.5f;
  }
}

void DragCtrl::selectDragEle(const vector<int> ids){
  
  if(_selected && ids.size() > 0){
	if(_observer)
	  _observer->selectDragEle(ids[0]);
	emit selectDragEle(ids[0]);
  }else{
	begin_drag = false;
	is_dragging = false;
  }
  _selected = false;
}
