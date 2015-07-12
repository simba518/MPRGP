#ifndef _BBOX_H_
#define _BBOX_H_

#include <vector>
#include <eigen3/Eigen/Dense>

namespace UTILITY{

  template<typename real,typename VECTOR3D>
  class BBox{

  public:
	BBox (){ clear();}
	BBox (const std::vector<real> &vertexs){
	  clear();
	  reset(vertexs);
	}
	BBox (const real *vertexs,const int vertex_num){
	  clear();
	  reset (vertexs,vertex_num);
	}
	void clear(){
	  setVec3d(_leftDownBack,0,0,0);
	  setVec3d(_rightUpFront,0,0,0);
	  setVec3d(_center,0,0,0);
	}
	void reset(const std::vector<real> &vertexs){
	  if(vertexs.size() > 0){
		int vertexs_num = (int)vertexs.size()/3;
		const real *p_vertexs = &vertexs[0];
		reset(p_vertexs,vertexs_num);
	  }
	}
	void reset(const real *vertexs,const int vertex_num){
	  if (vertexs == NULL || vertex_num <= 0){
		clear();
		return ;
	  }
	  setVec3d(_leftDownBack,vertexs[0],vertexs[1],vertexs[2]);
	  setVec3d(_rightUpFront,vertexs[0],vertexs[1],vertexs[2]);
	  real v[3];
	  for (int i=1; i < vertex_num; i++){
		v[0] = vertexs[i*3+0];
		v[1] = vertexs[i*3+1];
		v[2] = vertexs[i*3+2];
		_leftDownBack[0] = minValue(_leftDownBack[0],v[0]);
		_leftDownBack[1] = minValue(_leftDownBack[1],v[1]);
		_leftDownBack[2] = minValue(_leftDownBack[2],v[2]);
		_rightUpFront[0] = maxValue(_rightUpFront[0],v[0]);
		_rightUpFront[1] = maxValue(_rightUpFront[1],v[1]);
		_rightUpFront[2] = maxValue(_rightUpFront[2],v[2]);
	  }
	  _center = _rightUpFront;
	  _center +=_leftDownBack;
	  _center *= (1.0/2.0);
	}
	const VECTOR3D& getCenter()const{ return _center;}
	real getMaxAxis()const{
	  real m = maxValue(getWidth(),getHeight());
	  m = maxValue(m,getDeepth());
	  return m;
	}

	real getWidth()const{return _rightUpFront[0] - _leftDownBack[0];}
	real getDeepth()const{return _rightUpFront[2] - _leftDownBack[2];}
	real getHeight()const{return _rightUpFront[1] - _leftDownBack[1];}
	const VECTOR3D& getMaxConner()const{ return _rightUpFront;}
	const VECTOR3D& getMinConner()const{ return _leftDownBack;}
	void getMaxConner(real v[3])const{
	  v[0] = _rightUpFront[0];
	  v[1] = _rightUpFront[1];
	  v[2] = _rightUpFront[2];
	}
	void getMinConner(real v[3])const{
	  v[0] = _leftDownBack[0];
	  v[1] = _leftDownBack[1];
	  v[2] = _leftDownBack[2];
	}
	BBox& add(const BBox &box2){
	  const VECTOR3D &b2Left = box2.getMinConner();
	  const VECTOR3D &b2right = box2.getMaxConner();
	  _leftDownBack[0] = minValue(_leftDownBack[0],b2Left[0]);
	  _leftDownBack[1] = minValue(_leftDownBack[1],b2Left[1]);
	  _leftDownBack[2] = minValue(_leftDownBack[2],b2Left[2]);
	  _rightUpFront[0] = maxValue(_rightUpFront[0],b2right[0]);
	  _rightUpFront[1] = maxValue(_rightUpFront[1],b2right[1]);
	  _rightUpFront[2] = maxValue(_rightUpFront[2],b2right[2]);
	  return *this;
	}

  protected:
	void setVec3d(VECTOR3D &v,real v1,real v2,real v3){
	  v[0] = v1; 	  v[1] = v2;      v[2] = v3;
	}
	real maxValue(real a,real b)const{ return a > b ? a:b;}
	real minValue(real a,real b)const{ return a < b ? a:b;}

  private:
	VECTOR3D _leftDownBack;//min conner
	VECTOR3D _rightUpFront;//max conner
	VECTOR3D _center;
  };

  typedef BBox<double,Eigen::Vector3d> BBoxD;

}

#endif /* _BBOX_H_ */
