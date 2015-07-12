#include <WindowsHeaders.h>
#include <GL/glew.h>
#include <math.h>
#include "AxisTorus.h"
#include <QGLViewer/config.h>
#include <iostream>
#include <Log.h>
using namespace std;
using namespace QGLVEXT;

void AxisTorus::draw()const{

  if (selected_axis == 0) {
	drawTorus(0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f); // x
  }else{
	drawTorus(0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f); // x
  }if (selected_axis == 1){
	drawTorus(0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f); // y
  }else{
	drawTorus(0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f); // y
  }if (selected_axis == 2){
	drawTorus(1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f); // z
  }else{
	drawTorus(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f); // z
  }if (selected_axis == 3){
	drawAxis(1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f);  // x
  }else{
	drawAxis(1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f);  // x
  }if (selected_axis == 4){
	drawAxis(0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f);  // y
  }else{
	drawAxis(0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f);  // y
  }if (selected_axis == 5){
	drawAxis(0.0f, -1.0f, 0.0f, 1.0f, 1.0f, 0.0f); // z
  }else{
	drawAxis(0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f); // z
  }if (selected_axis == 6){
	drawOriginal(1.0f,1.0f,0.0f);  // o
  }else{
	drawOriginal(1.0f,0.0f,0.0f);  // o
  }
}

void AxisTorus::drawWithNames ()const{

  glFlush();

  glPushName(0);
  drawTorus(0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f);
  glPopName();

  glPushName(1);
  drawTorus(0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f);
  glPopName();

  glPushName(2);
  drawTorus(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
  glPopName();

  glPushName(3);
  drawAxis(1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f);
  glPopName();

  glPushName(4);
  drawAxis(0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f);
  glPopName();

  glPushName(5);
  drawAxis(0.0f, -1.0f, 0.0f, 1.0f, 1.0f, 0.0f);
  glPopName();

  glPushName(6);
  drawOriginal(1.0f, 0.0f, 0.0f);
  glPopName();
}

void AxisTorus::drawTorus(float x, float y, float z,
						  float r, float g, float b)const{

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glDisable(GL_LIGHTING);
  glColor4f(r,g,b,1.0f);

  glPushMatrix();
  glTranslated(translation[0], translation[1], translation[2]);
  glScaled(scalor[0], scalor[1], scalor[2]);

  glRotatef(90.0f, x, y, z);
  glLineWidth(1);
  const double TWOPI = 2 * M_PI;
  const int numc = 100;
  glBegin(GL_LINE_LOOP);
  for (int i = 0; i < numc; i++){
	glNormal3d(5*cos(1.0*i*TWOPI/numc), 0.0, 5*sin(1.0*i*TWOPI/numc));
	glVertex3d(5*cos(1.0*i*TWOPI/numc), 0.0, 5*sin(1.0*i*TWOPI/numc));
  }
  glEnd();

  glPopMatrix();
  glPopAttrib();
}

void AxisTorus::drawAxis(float x, float y, float z,
						 float r, float g, float b)const{
  
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glColor4f(r,g,b,1.0f);

  glPushMatrix();
  glTranslated(translation[0], translation[1], translation[2]);
  glScaled(scalor[0], scalor[1], scalor[2]);

  glRotatef(90.0f, x, y, z);

  glLineWidth(4);
  glBegin(GL_LINES);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(5.0f, 0.0f, 0.0f);
  glEnd();

  glPopMatrix();
  glPopAttrib();
}

void AxisTorus::drawOriginal(float r, float g, float b)const{

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glColor4f(r,g,b,1.0f);

  glPushMatrix();
  glTranslated(translation[0], translation[1], translation[2]);
  glScaled(scalor[0], scalor[1], scalor[2]);

  gluSphere(gluNewQuadric(), 0.8f, 30,30);

  glPopMatrix();
  glPopAttrib();
}
