#ifndef _ANIDATAMODEL_H_
#define _ANIDATAMODEL_H_

#include <string>
#include <boost/shared_ptr.hpp>
using namespace std;

namespace QGLVEXT{

/*default step size for fast forward*/
#define DEFAULT_FF_STEP_SIZE 5 
  
  /**
   * @class AniDataModel data model for animatoin displaying,loading and
   * saving.
   */
  class AniDataModel{
	
  public:
	// initialize
	AniDataModel(){
	  current_frame = 0;
	  sub_ani_step = 1;
	  sub_ani_start = 0;
	  sub_ani_end = -1;
	  ff_step_size = DEFAULT_FF_STEP_SIZE;
	  total_frame_num = 0;
	}
	void setTotalFrameNum(const int num){
	  total_frame_num = num > 0 ? num:0;
	}
	void setFFSpeed(const int step_size){
	  ff_step_size = step_size>0? ff_step_size:DEFAULT_FF_STEP_SIZE;
	}
	void setSubAniStep(const int num){
	  sub_ani_step = num>0? num:1;
	}
	void setSubAniEnd(const int num){
	  sub_ani_end = num>0 ? num:totalFrameNum();
	}
	void setSubAniStart(const int num){
	  sub_ani_start = num>0? num:0;
	}	

	// get data
	int currentFrameNum()const{
	  return warpFrameNum( current_frame );
	}
	int totalFrameNum()const{
	  return total_frame_num;
	}
	int endFrame()const{
	  if (sub_ani_end > 0){
		return (sub_ani_end >= totalFrameNum() ? totalFrameNum()-1: sub_ani_end);
	  }else{
		return totalFrameNum()-1;
	  }
	}
	int getSubAniStep()const{
	  return sub_ani_step;
	}
	int getSubAniEnd()const{
	  return sub_ani_end;
	}
	int getSubAniStart()const{
	  return sub_ani_start;
	}

	// animation control
	void setCurrentFrame(const int frame_num){
	  this->current_frame = warpFrameNum(frame_num);
	}
	void forward(){
	  current_frame += sub_ani_step;
	  this->setCurrentFrame(current_frame);
	}
	void backward(){
	  current_frame -= sub_ani_step;
	  this->setCurrentFrame(current_frame);
	}
	void fastForward(){
	  current_frame += ff_step_size;
	  this->setCurrentFrame(current_frame);	  
	}
	void fastBackward(){
	  current_frame -= ff_step_size;
	  this->setCurrentFrame(current_frame);	  
	}
	void reset(){
	    this->setCurrentFrame(sub_ani_start);
	}
	bool reachEnd()const{
	  return (currentFrameNum() >= endFrame());
	}

  protected:
	// return the valid frame number in [0,totalFrameNum-1].
	int warpFrameNum(const int frame_num)const{

	  int c_frame = frame_num;
	  if( frame_num >= endFrame() ){
		c_frame = endFrame();
	  }
	  if(frame_num < sub_ani_start){
		c_frame = sub_ani_start;
	  }
	  return c_frame;
	}

  protected:
	int total_frame_num;
	int current_frame;
	int ff_step_size; // step size of fast foward

	// sub animation
	int sub_ani_step; // step size for sub animation sequence.
	int sub_ani_start;// first frame for sub animation sequence.
	int sub_ani_end;  // last frame for sub animation sequence.
  };
  
  typedef boost::shared_ptr<AniDataModel> pAniDataModel;
  typedef boost::shared_ptr<const AniDataModel> pAniDataModel_const;
  
}//end of namespace

#endif /* _ANIDATAMODEL_H_ */
