#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <QtGui/QFileDialog>
#include <AniCtrl.h>
#include <assertext.h>
using namespace boost;
using namespace QGLVEXT;

AniCtrl::AniCtrl
(pAniDataModel data_model,pQGLViewerExt viewer)
 :viewer(viewer){

  play_circle = false;
  is_pause = false;

  assert (viewer != NULL);
  assert (data_model != NULL);

  addAniDataModel(data_model);

  connect(viewer,SIGNAL(resetAnimation()), this,SLOT(reset()));

  connect(viewer,SIGNAL(updateAnimation()), this,SLOT(update()));
  connect(viewer,SIGNAL(backwardAnimation()), this,SLOT(backward()));
  connect(viewer,SIGNAL(forwardAnimation()), this,SLOT(forward()));

  connect(viewer,SIGNAL(ffAnimation()), this,SLOT(fastForward()));
  connect(viewer,SIGNAL(fbAnimation()), this,SLOT(fastBackward()));

  viewer->stopAnimation();

}

int AniCtrl::currentFrameNum()const{
  
  int current_frame_num = -1;
  if(data_models.size() > 0){
	current_frame_num = data_models.front()->currentFrameNum();
  }
  return current_frame_num;
}

bool AniCtrl::allReachEnd()const{
  
  bool reach_end = true;
  BOOST_FOREACH(pAniDataModel data_model, data_models){
	if(!data_model->reachEnd()){
	  reach_end = false;
	  break;
	}
  }
  return reach_end;
}

void AniCtrl::removeAllAniDataModel(){
  
  if(viewer != NULL){
	viewer->stopAnimation();
  }
  updateObserveViewers();
  data_models.clear();
}

void AniCtrl::reset(){

  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->reset();
  }
  updateObserveViewers();
  emit currentFrame(this->currentFrameNum());
}

void AniCtrl::update(){

  if (isPause()){
	return ;
  }

  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->forward();
  }
  emit currentFrame(this->currentFrameNum());

  if(allReachEnd()){
	if (playCircle()){
	  reset();
	}else if (viewer != NULL){
	  viewer->pauseAnimation();
	}
  }
  updateObserveViewers();
}

void AniCtrl::forward(){

  if (isPause()){
	return ;
  }

  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->forward();
  }
  emit currentFrame(this->currentFrameNum());
  if (viewer != NULL){
	viewer->update();
	updateObserveViewers();
	if(allReachEnd()){
	  if (playCircle()){
		reset();
	  }else{
		viewer->pauseAnimation();
	  }
	}
  }
}

void AniCtrl::fastForward(){

  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->fastForward();
  }
  emit currentFrame(this->currentFrameNum());
  if (viewer != NULL){
	viewer->update();
  }
  updateObserveViewers();
}

void AniCtrl::backward(){
	  
  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->backward();
  }
  if(viewer != NULL){
	viewer->update();
  }
  emit currentFrame(this->currentFrameNum());
}

void AniCtrl::fastBackward(){
	  
  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->fastBackward();
  }
  if(viewer != NULL){
	viewer->update();
  }
  emit currentFrame(this->currentFrameNum());
}

void AniCtrl::setCurrentFrame(const int frame_num){
  
  BOOST_FOREACH(pAniDataModel data_model, data_models){
	data_model->setCurrentFrame(frame_num);
  }
  if (viewer != NULL){
	viewer->update();
  }
}

void AniCtrl::setAnimationPeriod(const int period){
  
  assert_gt (period,0);
  if (viewer != NULL)
	viewer->setAnimationPeriod(period);
}

void AniCtrl::setFPS(const int fps){
  
  assert_gt (fps , 0);
  const int period = 1000/fps;
  this->setAnimationPeriod(period);
}

int AniCtrl::totalFrameNum()const{

  int T = 0;
  BOOST_FOREACH(pAniDataModel data_model, data_models){
	if ( data_model->totalFrameNum() > T){
	  T = data_model->totalFrameNum();
	}
  }  
  return T;
}

void AniCtrl::updateObserveViewers(){
  
  BOOST_FOREACH(pQGLViewerExt ele, observe_viewers){
	if(ele)
	  ele->update();
  }
}
