#include <boost/foreach.hpp>
#include "DragGroupSelDraw.h"

#include <GL/glew.h>
#include <GL/gl.h>
using namespace UTILITY;

void DragGroupSelDraw::drawAllGroupsWithPoints(pTetMesh_const vol_mesh,
											   const vector<set<int> >&drag_groups, 
											   const double*vol_u){

  glFlush();
  int i = 0;
  BOOST_FOREACH(const set<int> &group, drag_groups){
	glPushName(i);
	drawPoints(vol_mesh, group, vol_u);
	glPopName();
	++i;
  }
}

void DragGroupSelDraw::drawAllGroupsWithShpere(pTetMesh_const vol_mesh,
											   const vector<set<int> >&drag_groups, 
											   const double *vol_u, 
											   const double radius){

  glFlush();
  int i = 0;
  BOOST_FOREACH(const set<int> &group, drag_groups){
	glPushName(i);
	drawSphere(vol_mesh, group, vol_u, radius);
	glPopName();
	++i;
  }
}

void DragGroupSelDraw::drawPoints(pTetMesh_const vol_mesh,
								  const set<int> &group,
								  const double *vol_u){

  glBegin(GL_POINTS);
  BOOST_FOREACH(int vertex_id, group ){
	const Vector3d &v = vol_mesh->nodes()[vertex_id];
	const double x = v[0] + vol_u[vertex_id*3+0];
	const double y = v[1] + vol_u[vertex_id*3+1];
	const double z = v[2] + vol_u[vertex_id*3+2];		
	glVertex3d(x,y,z);
  }
  glEnd();  
}

void DragGroupSelDraw::drawSphere(pTetMesh_const vol_mesh,
								  const set<int> &group,
								  const double *vol_u, const double radius){

  assert_gt (radius, 0.0f);
  double x = 0;
  double y = 0;
  double z = 0;
  BOOST_FOREACH(int vertex_id, group ){
	const Vector3d &v = vol_mesh->nodes()[vertex_id];
	x += v[0] + vol_u[vertex_id*3+0];
	y += v[1] + vol_u[vertex_id*3+1];
	z += v[2] + vol_u[vertex_id*3+2];
  }
  if (group.size() > 0){
	x /= (double)group.size();
	y /= (double)group.size();
	z /= (double)group.size();
	assert_gt (radius, 0.0f);
	DrawSphere(x,y,z,radius,10,10);
  }
}

void DragGroupSelDraw::DrawSphere(float x, float y, float z, float fRadius, int M, int N){

  glPushMatrix();
  glTranslatef(x, y, z);
  gluSphere(gluNewQuadric(), fRadius, M, N);
  glPopMatrix();
}
