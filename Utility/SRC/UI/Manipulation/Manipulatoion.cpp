#include "Manipulatoion.h"
using namespace QGLVEXT;

LocalframeManipulatoion::LocalframeManipulatoion(pQGLViewerExt v):viewer(v){

  frame = new ManipulatedFrame();
  current_pos.setZero();
  con_track_ball = pConTrackBall(new ConTrackBall(viewer,frame));
  con_track_ball->setShow(false);
  con_track_ball->setEnable(false);
  connect( frame, SIGNAL(manipulated()), this, SLOT(manipulate()) );
  connect( con_track_ball.get(), SIGNAL(selectObject(const int)), 
		   this, SLOT(selectConTrackball(const int )) );
}

void LocalframeManipulatoion::setShow(const bool show){

  con_track_ball->setShow(show);
  con_track_ball->setEnable(show);
  viewer->setManipulatedFrame(viewer->camera()->frame());

  if (isShown()){

	computeCurrentPosition(current_pos);
	frame->setTranslation(Vec(current_pos[0],current_pos[1],current_pos[2]));
	frame->setOrientation(qglviewer::Quaternion(Vec(1.0,0.0,0.0), 0.0));

	con_track_ball->translate(current_pos[0],current_pos[1],current_pos[2]);
	const double s = viewer->sceneRadius()/15.0f;
	con_track_ball->scale(s,s,s);
  }
}

void LocalframeManipulatoion::manipulate(){
  computeCurrentPosition(current_pos);
  con_track_ball->translate(current_pos[0],current_pos[1],current_pos[2]);
  this->applyTransform();
}

void LocalframeManipulatoion::selectConTrackball(const int sel_axis){

  if (sel_axis >= 0 && isShown()){
	viewer->setManipulatedFrame(frame);
  }else{
	viewer->setManipulatedFrame(viewer->camera()->frame());
  }
}

void LocalframeManipulatoionExt::setShow(const bool show){

  LocalframeManipulatoion::setShow(show);
  if (show){
	computeCurrentPosition(current_pos);
	const Vector3d xc = current_pos;
	initial_pos = this->getCurrentPositions();
	for (int i = 0; i < initial_pos.cols(); ++i){
	  initial_pos.col(i) -= xc;
	}
  }
}

void LocalframeManipulatoionExt::computeCurrentPosition(Vector3d &current_pos)const{

  const Matrix<double,3,-1> &xc = this->getCurrentPositions();
  current_pos.setZero();
  for (int i = 0; i < xc.cols(); ++i)
	current_pos += xc.col(i);
  if (xc.cols() > 0)
	current_pos *= 1.0f/xc.cols();
}

void LocalframeManipulatoionExt::applyTransform(){

  Matrix<double,3,-1> &x = this->getCurrentPositions();
  assert_eq(x.cols(), initial_pos.cols());
  assert_eq(x.rows(),3);
  assert_eq(initial_pos.rows(),3);
  for (int i = 0; i < initial_pos.cols(); ++i){
	const Vec p0(initial_pos(0,i), initial_pos(1,i), initial_pos(2,i));
	const Vec p = frame->inverseCoordinatesOf(p0);
	x(0,i) = p[0];
	x(1,i) = p[1];
	x(2,i) = p[2];
  }
}
