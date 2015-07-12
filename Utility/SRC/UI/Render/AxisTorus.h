#ifndef _AXISTORUS_H_
#define _AXISTORUS_H_

#include <boost/shared_ptr.hpp>
#include <SelfRenderEle.h>
#include <Selectable.h>

namespace QGLVEXT{
  
  /**
   * @class AxisTorus
   * 
   */
  class AxisTorus:public SelfRenderEle, public Selectable{
	
  public: 
	AxisTorus(){

	  selected_axis = -1;

	  translation[0] = 0.0;
	  translation[1] = 0.0;
	  translation[2] = 0.0;

	  scalor[0] = 1.0;
	  scalor[1] = 1.0;
	  scalor[2] = 1.0;
	}
	void translate(const double x,const double y, const double z){
	  this->translation[0] = x;
	  this->translation[1] = y;
	  this->translation[2] = z;
	}
	void scale(const double x,const double y, const double z){
	  this->scalor[0] = x;
	  this->scalor[1] = y;
	  this->scalor[2] = z;
	}

	void draw()const;
	int totalEleNum ()const{
	  return 8;
	}
	void drawWithNames ()const;
	void selectAxis(int a){
	  selected_axis = a;
	}
	int selectedAxis()const{
	  return selected_axis;
	}
	
  protected:
	void drawTorus(float x, float y, float z,
				   float r, float g, float b)const;
	void drawAxis(float x, float y, float z,
				  float r, float g, float b)const;
	void drawOriginal(float r, float g, float b)const;
	
  private:
	int selected_axis;
	double translation[3], scalor[3];
  };
  
  typedef boost::shared_ptr<AxisTorus> pAxisTorus;
  
}//end of namespace

#endif /* _AXISTORUS_H_ */
