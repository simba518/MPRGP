#include <WindowsHeaders.h>
#include <GL/glew.h>
#include <math.h>
#include "ConTrackBall.h"
#include <assertext.h>
using namespace qglviewer;
using namespace QGLVEXT;

void ConTrackBall::selectAxises(const vector<int> sel_group_ids){

  if (m_enabled && viewer->getSelector().get() == p_AxisTorus.get()){

	if (sel_group_ids.size() <= 0 ){
	  constrained_axi = -1;
	  p_AxisTorus->selectAxis(constrained_axi);
	}else{
	  constrained_axi = sel_group_ids[0];
	  p_AxisTorus->selectAxis(constrained_axi);
	  viewer->update();
	}
  }
  m_hit = (p_AxisTorus->selectedAxis() >= 0);
  emit selectObject(p_AxisTorus->selectedAxis());
}

void ConTrackBall::checkIfGrabsMouse(int x,int y,const Camera*const camera){

  if (m_enabled){
	viewer->setSelector(p_AxisTorus);
	QRect select_rect(x,y,5,5);
	viewer->select (select_rect);
  }
}

void ConTrackBall::press(QMouseEvent* e){

  if (m_enabled && m_hit && constrained_axi >= 0){

	assert_in(constrained_axi,0,6);
	Vec dir(0.0,0.0,0.0);
	dir[constrained_axi%3] = 1.0;
	if (constrained_axi < 3){
	  constraint->setTranslationConstraintType(AxisPlaneConstraint::FORBIDDEN);
	  constraint->setRotationConstraintType(AxisPlaneConstraint::AXIS);
	  constraint->setRotationConstraintDirection(dir);
	}else if(constrained_axi < 6){
	  constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
	  constraint->setTranslationConstraintType(AxisPlaneConstraint::AXIS);
	  constraint->setTranslationConstraintDirection(dir);
	}else{
	  constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
	  constraint->setTranslationConstraintType(AxisPlaneConstraint::FREE);
	}
  }
}

void ConTrackBall::release(QMouseEvent* e){

  if (m_enabled){
	constraint->setRotationConstraintType(AxisPlaneConstraint::FREE);
	constraint->setTranslationConstraintType(AxisPlaneConstraint::FREE);
	constrained_axi = -1;
	p_AxisTorus->selectAxis(constrained_axi);
  }
}

void ConTrackBall::init(){
  
  assert(viewer);

  constrained_axi = 0;
  m_show = false;
  m_hit = false;
  m_enabled = false;
  p_AxisTorus = pAxisTorus(new AxisTorus());

  constraint = new WorldConstraint();
  constraint->setRotationConstraintType(AxisPlaneConstraint::FREE);

  connect( viewer,SIGNAL(mousePressSignal(QMouseEvent *)),this, SLOT(press(QMouseEvent *)));
  connect( viewer,SIGNAL(mouseReleaseSignal(QMouseEvent*)),this,SLOT(release(QMouseEvent*)));
  connect( viewer,SIGNAL(selectedIds(vector<int>)),this, SLOT(selectAxises(vector<int>)));
}
