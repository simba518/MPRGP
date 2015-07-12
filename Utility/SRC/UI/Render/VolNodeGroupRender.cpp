#include <WindowsHeaders.h>
#include <boost/foreach.hpp>
#include "VolNodeGroupRender.h"
#include <PhongShader.h>
#include <GL/gl.h>
#include <GL/glew.h>
using namespace UTILITY;

VolNodeGroupRender::VolNodeGroupRender(){
  
  //set defualt color
  color[0] = 0.8f;
  color[1] = 0.1f;
  color[2] = 0.2f;
  color[3] = 1.0f;
  size = 8.0f;
  radius = 0.05f;
}

VolNodeGroupRender::VolNodeGroupRender(const double color [4],double size){

  setColor(color);
  setPointSize(size);
  setSphereRadius(0.002f);
}

void VolNodeGroupRender::setColor(const double c[4]){

  assert (c[0] >=0.0f && c[0] <= 1.0f);
  assert (c[1] >=0.0f && c[1] <= 1.0f);
  assert (c[2] >=0.0f && c[2] <= 1.0f);
  assert (c[3] >=0.0f && c[3] <= 1.0f);

  color[0] = c[0];
  color[1] = c[1];
  color[2] = c[2];
  color[3] = c[3];
}

void VolNodeGroupRender::drawPoints(pTetMesh_const vol_mesh, 
									int v_id, const double *u)const{

  glDisable(GL_LIGHTING);
  glPointSize(size);
  glColor4d(color[0], color[1], color[2],color[3]);
  glBegin(GL_POINTS);

  const Vector3d &v = vol_mesh->nodes()[v_id];
  if(u != NULL){

	double x = v[0]+u[v_id*3+0];
	double y = v[1]+u[v_id*3+1];
	double z = v[2]+u[v_id*3+2];
	glVertex3d(x,y,z);
  }else{
	glVertex3d(v[0],v[1],v[2]);
  }
  glEnd();

  glEnable(GL_LIGHTING);
}

void VolNodeGroupRender::drawShpere(pTetMesh_const vol_mesh, 
									int v_id, const double *u)const{

  glEnable(GL_LIGHTING);
  const Vector3d &v = vol_mesh->nodes()[v_id];
  if(u != NULL){
	const double x = v[0]+u[v_id*3+0];
	const double y = v[1]+u[v_id*3+1];
	const double z = v[2]+u[v_id*3+2];
	DrawSphere(x,y,z,radius,10,10);
  }else{
	DrawSphere(v[0],v[1],v[2],radius,10,10);
  }
}

void VolNodeGroupRender::drawPoints(pTetMesh_const vol_mesh, 
									const vector<int> &group, const double *u)const{

  glDisable(GL_LIGHTING);
  glPointSize(size);
  glColor4d(color[0], color[1], color[2],color[3]);
  glBegin(GL_POINTS);
  BOOST_FOREACH(int v_id, group){

	const Vector3d &v = vol_mesh->nodes()[v_id];
	if(u != NULL){
	  double x = v[0]+u[v_id*3+0];
	  double y = v[1]+u[v_id*3+1];
	  double z = v[2]+u[v_id*3+2];
	  glVertex3d(x,y,z);
	}else{
	  glVertex3d(v[0],v[1],v[2]);
	}
  }
  glEnd();
  glEnable(GL_LIGHTING); 
}

void VolNodeGroupRender::drawShpere(pTetMesh_const vol_mesh, 
									const vector<int> &group, const double *u)const{

  if (group.size() <= 0){
	return ;
  }

  double x = 0;
  double y = 0;
  double z = 0;
  glEnable(GL_LIGHTING);
  BOOST_FOREACH(int v_id, group){

	const Vector3d &v = vol_mesh->nodes()[v_id];
	if(u != NULL){
	  x += v[0]+u[v_id*3+0];
	  y += v[1]+u[v_id*3+1];
	  z += v[2]+u[v_id*3+2];
	}else{
	  x += v[0];
	  y += v[0];
	  z += v[0];
	}
  }

  x /= (double)(group.size());
  y /= (double)(group.size()); 
  z /= (double)(group.size());
  DrawSphere(x,y,z,radius,10,10);
}

void VolNodeGroupRender::drawPoints(pTetMesh_const vol_mesh, 
									const set<int> &group, const double *u)const{
  
  glDisable(GL_LIGHTING);
  glPointSize(size);
  glColor4d(color[0], color[1], color[2],color[3]);
  glBegin(GL_POINTS);
  BOOST_FOREACH(int v_id, group){

	const Vector3d &v = vol_mesh->nodes()[v_id];
	if(u != NULL){
	  double x = v[0]+u[v_id*3+0];
	  double y = v[1]+u[v_id*3+1];
	  double z = v[2]+u[v_id*3+2];
	  glVertex3d(x,y,z);
	}else{
	  glVertex3d(v[0],v[1],v[2]);
	}
  }
  glEnd();
  glEnable(GL_LIGHTING); 
}

void VolNodeGroupRender::drawShpere(pTetMesh_const vol_mesh, 
									const set<int> &group, const double *u)const{
  
  if (group.size() <= 0){
	return ;
  }

  double x = 0;
  double y = 0;
  double z = 0;
  BOOST_FOREACH(int v_id, group){

	const Vector3d &v = vol_mesh->nodes()[v_id];
	if(u != NULL){
	  x += v[0]+u[v_id*3+0];
	  y += v[1]+u[v_id*3+1];
	  z += v[2]+u[v_id*3+2];
	}else{
	  x += v[0];
	  y += v[0];
	  z += v[0];
	}
  }
  x /= (double)(group.size());
  y /= (double)(group.size()); 
  z /= (double)(group.size());

  UTILITY::PhongShader::getInstance()->bind();

  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_LIGHTING_BIT);

  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glColor4d(color[0], color[1], color[2],color[3]);
  DrawSphere(x,y,z,radius,40,40);
  glDisable(GL_COLOR_MATERIAL);

  glPopAttrib();

  UTILITY::PhongShader::getInstance()->unbind();
}

void VolNodeGroupRender::DrawSphere(float x, float y, float z, float fRadius, int M, int N)const{
  
  glPushMatrix();
  glTranslatef(x, y, z);
  gluSphere(gluNewQuadric(), fRadius, M, N);
  glPopMatrix();
}
