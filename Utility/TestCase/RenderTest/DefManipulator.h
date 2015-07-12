#ifndef _DEFMANIPULATOR_H_
#define _DEFMANIPULATOR_H_

#include <Manipulatoion.h>
using namespace QGLVEXT;

// class Coord{
  
// public:
//   Coord(){}
//   Coord(Vec o, Vec x, Vec y, Vec z):o(o),x(x),y(y),z(z){}
//   void draw()const{
	
// 	glDisable(GL_LIGHTING);
// 	glPointSize(10);
// 	glColor3f(1.0, 0.0, 0.0);

// 	glBegin(GL_POINTS);
// 	glVertex3f(new_o[0], new_o[1], new_o[2]);
// 	glVertex3f(new_x[0], new_x[1], new_x[2]);
// 	glVertex3f(new_y[0], new_y[1], new_y[2]);
// 	glVertex3f(new_z[0], new_z[1], new_z[2]);
// 	glEnd();

// 	glBegin(GL_LINES);
// 	glVertex3f(new_o[0], new_o[1], new_o[2]);
// 	glVertex3f(new_x[0], new_x[1], new_x[2]);

// 	glVertex3f(new_o[0], new_o[1], new_o[2]);
// 	glVertex3f(new_y[0], new_y[1], new_y[2]);

// 	glVertex3f(new_o[0], new_o[1], new_o[2]);
// 	glVertex3f(new_z[0], new_z[1], new_z[2]);
// 	glEnd();
//   }

// public:
//   Vec o,x,y,z;
//   Vec new_o,new_x,new_y,new_z;
// };

// class DefManipulator:public LocalframeManipulatoion{
	
// public:
//   DefManipulator(pQGLViewerExt v):LocalframeManipulatoion(v){

// 	coord[0].new_o = coord[0].o = Vec(0,0,0);
// 	coord[0].new_x = coord[0].x = Vec(1,0,0);
// 	coord[0].new_y = coord[0].y = Vec(0,1,0);
// 	coord[0].new_z = coord[0].z = Vec(0,0,1);
// 	coord[1] = coord[0];
// 	coord[0].new_x = coord[1].x = Vec(-1,0,0);
// 	sel_coord = -1;
//   }

//   void setEnable(const bool enable){

// 	LocalframeManipulatoion::setEnable(enable);
// 	if(enable){
// 	  assert_in(sel_coord,0,1);
// 	  double pos_xyz[3];
// 	  currentPosition(pos_xyz);
// 	  const Vec origin(pos_xyz[0], pos_xyz[1], pos_xyz[2]);
// 	  coord[sel_coord].o = coord[sel_coord].new_o - origin;
// 	  coord[sel_coord].x = coord[sel_coord].new_x - origin;
// 	  coord[sel_coord].y = coord[sel_coord].new_y - origin;
// 	  coord[sel_coord].z = coord[sel_coord].new_z - origin;
// 	}
//   }
//   void draw()const{
	
// 	LocalframeManipulatoion::draw();
// 	coord[0].draw();
// 	coord[1].draw();
//   }

//   // select
//   int totalEleNum()const{
// 	return 2;
//   }
//   void drawWithNames()const{
	
// 	glPushName(0);
// 	coord[0].draw();
// 	glPopName();

// 	glPushName(1);
// 	coord[1].draw();
// 	glPopName();
//   }
//   void select(const vector<int> &sel_ids){
// 	if (sel_ids.size() <= 0){
// 	  sel_coord = -1;
// 	}else{
// 	  sel_coord = sel_ids[0];
// 	}
//   }

//   // manipulate
//   void currentPosition(double pos_xyz[3])const{
// 	assert_in(sel_coord,0,1);
// 	pos_xyz[0] = coord[sel_coord].new_o[0];
// 	pos_xyz[1] = coord[sel_coord].new_o[1];
// 	pos_xyz[2] = coord[sel_coord].new_o[2];
//   }
//   void applyTransform(){
	
// 	assert_in(sel_coord,0,1);
// 	coord[sel_coord].new_o = frame->inverseCoordinatesOf(coord[sel_coord].o);
// 	coord[sel_coord].new_x = frame->inverseCoordinatesOf(coord[sel_coord].x);
// 	coord[sel_coord].new_y = frame->inverseCoordinatesOf(coord[sel_coord].y);
// 	coord[sel_coord].new_z = frame->inverseCoordinatesOf(coord[sel_coord].z);
//   }

// private:
//   int sel_coord;
//   Coord coord[2];
// };

class DefManipulatorExt:public LocalframeManipulatoionExt{
	
public:
  DefManipulatorExt(pQGLViewerExt v):LocalframeManipulatoionExt(v){
	sel_coord = -1;
	pos[0] = Matrix<double,3,-1>(3,4);
	pos[0] << 
	  0,1,0,0,
	  1,0,1,0,
	  0,0,0,1;
	pos[1] = pos[0];
	pos[1].col(3) *= -1.0f;
  }

  void draw()const{
	LocalframeManipulatoion::draw();
	drawCoords(pos[0]);
	drawCoords(pos[1]);
  }
  void drawCoords(const Matrix<double,3,-1> &coord)const{

	const Vector3d o = coord.col(0);
	const Vector3d x = coord.col(1);
	const Vector3d y = coord.col(2);
	const Vector3d z = coord.col(3);

	glDisable(GL_LIGHTING);
	glPointSize(10);
	glColor3f(1.0, 0.0, 0.0);

	glBegin(GL_POINTS);
	glVertex3f(o[0], o[1], o[2]);
	glVertex3f(x[0], x[1], x[2]);
	glVertex3f(y[0], y[1], y[2]);
	glVertex3f(z[0], z[1], z[2]);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(o[0], o[1], o[2]);
	glVertex3f(x[0], x[1], x[2]);

	glVertex3f(o[0], o[1], o[2]);
	glVertex3f(y[0], y[1], y[2]);

	glVertex3f(o[0], o[1], o[2]);
	glVertex3f(z[0], z[1], z[2]);
	glEnd();
  }

  // select
  int totalEleNum()const{
	return 2;
  }
  void drawWithNames()const{

	glPushName(0);
	drawCoords(pos[0]);
	glPopName();

	glPushName(1);
	drawCoords(pos[1]);
	glPopName();
  }
  void select(const vector<int> &sel_ids){
	if (sel_ids.size() <= 0){
	  sel_coord = -1;
	}else{
	  sel_coord = sel_ids[0];
	}
  }

  Matrix<double,3,-1> &getCurrentPositions(){
	assert_in(sel_coord,0,1);
	return pos[sel_coord];
  }
  const Matrix<double,3,-1> &getCurrentPositions()const{
	assert_in(sel_coord,0,1);
	return pos[sel_coord];
  }

private:
  int sel_coord;
  Matrix<double,3,-1> pos[2];
};

#endif /* _DEFMANIPULATOR_H_ */
