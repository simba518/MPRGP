#ifndef _ANISLIDERBAR_H_
#define _ANISLIDERBAR_H_

#include <boost/shared_ptr.hpp>
#include <AniCtrl.h>
#include <QSlider>

namespace QGLVEXT{
  
  /**
   * @class AniSliderBar
   * 
   */
  class AniSliderBar: public QSlider{

	Q_OBJECT
	
  public:
	AniSliderBar(QWidget *parent = 0):QSlider(parent){
  
	  QSlider::setMinimum(0);  
	  QSlider::setValue(0);
	  connect(this, SIGNAL(valueChanged(int)), this, SLOT(updateCtrlFrame(int)));
	}
	void setAniCtrl(pAniCtrl ctrl){

	  this->ctrl = ctrl;
	  QSlider::setValue(0);
	  if (ctrl != NULL){
		connect(ctrl.get(), SIGNAL(currentFrame(int)), this, SLOT(setCurrentFrame(int)));
	  }
	}
	
  public slots:
	void setCurrentFrame(const int frame_num){
  
	  if (ctrl != NULL){
		QSlider::setMaximum(ctrl->totalFrameNum());
		QSlider::setValue(frame_num);
	  }
	}
	void updateCtrlFrame(const int frame_num){
  
	  if (ctrl != NULL){
		QSlider::setMaximum(ctrl->totalFrameNum());
		ctrl->setCurrentFrame(frame_num);
	  }
	}

  private:
	pAniCtrl ctrl;

  };
  
  typedef boost::shared_ptr<AniSliderBar> pAniSliderBar;
  
}//end of namespace

#endif /*_ANISLIDERBAR_H_*/
