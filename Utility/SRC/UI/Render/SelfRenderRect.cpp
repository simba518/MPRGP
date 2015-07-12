#include <iostream>
using namespace std;
#include <WindowsHeaders.h>
#include <GL/gl.h>
#include "SelfRenderRect.h"

using namespace QGLVEXT;

void SelfRenderRect::draw()const{

  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT);
  if(p_qglviewer != NULL && !(this->isNull()) ){

	p_qglviewer->startScreenCoordinatesSystem();

	glEnable(GL_BLEND);
	glDisable(GL_LIGHTING);

	glColor4f(0.0, 0.0, 0.8f, 0.0f);
	glBegin(GL_QUADS);
	glVertex2i(this->left(),  this->top());
	glVertex2i(this->right(), this->top());
	glVertex2i(this->right(), this->bottom());
	glVertex2i(this->left(),  this->bottom());
	glEnd();

	glDisable(GL_BLEND);
	glLineWidth(6.0);
	glColor4f(0.8f, 0.0f, 0.8f, 0.8f);
	glBegin(GL_LINE_LOOP);
	glVertex2i(this->left(),  this->top());
	glVertex2i(this->right(), this->top());
	glVertex2i(this->right(), this->bottom());
	glVertex2i(this->left(),  this->bottom());
	glEnd();

	glPointSize(6.0);
	glBegin(GL_POINTS);
	glVertex2i(this->left(),  this->top());
	glVertex2i(this->right(), this->top());
	glVertex2i(this->right(), this->bottom());
	glVertex2i(this->left(),  this->bottom());
	glEnd();

	glEnable(GL_LIGHTING);

	p_qglviewer->stopScreenCoordinatesSystem();	
  }
  glPopAttrib();
}
