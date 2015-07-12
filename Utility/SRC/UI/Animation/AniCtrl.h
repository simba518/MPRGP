#ifndef _ANICTRL_H_
#define _ANICTRL_H_

#include <list>
using namespace std;

#include <QObject>
#include <boost/shared_ptr.hpp>
#include <AniDataModel.h>
#include <QGLViewerExt.h>

namespace QGLVEXT{
  
  /**
   * @class AniCtrl base animation controller
   * 
   */
  class AniCtrl:public QObject{

	Q_OBJECT
	
  public:
	AniCtrl(pAniDataModel data_model,pQGLViewerExt viewer);
	bool init(const string init_filename);
	void addObserveViewers(pQGLViewerExt v){observe_viewers.push_back(v);}
	void addAniDataModel(pAniDataModel dm){data_models.push_back(dm);}
	void removeAllAniDataModel();
	int currentFrameNum()const;
	int totalFrameNum()const;
	bool allReachEnd()const;
	void setPlayCircle(bool play_circle){this->play_circle = play_circle;}
	bool playCircle()const{ return play_circle;}
	bool isPause()const{ return is_pause;}

  public slots: 
	void reset();
	void update();

	void forward();
	void fastForward();

	void backward();
	void fastBackward();

	void setCurrentFrame(const int frame_num);

	void setAnimationPeriod(const int period); //in milliseconds
	void setFPS(const int fps); //frames per second

	void play(){ is_pause = false;}
	void pause(){ is_pause = true;}
	void togglePlay(){ is_pause = is_pause ? false:true;}
	void togglePlayCircle(){ play_circle = play_circle ? false:true;}
	void updateObserveViewers();

  signals:
	void currentFrame(int frame_num);

  protected:
	pQGLViewerExt viewer;
	list<pAniDataModel> data_models;
	list<pQGLViewerExt> observe_viewers;
	bool play_circle;
	bool is_pause;
  };
  
  typedef boost::shared_ptr<AniCtrl> pAniCtrl;
  
}//end of namespace

#endif /* _ANICTRL_H_ */
