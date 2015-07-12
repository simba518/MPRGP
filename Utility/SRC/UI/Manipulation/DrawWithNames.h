#ifndef _DRAWWITHNMAES_H_
#define _DRAWWITHNMAES_H_

#include <iostream>
using namespace std;
#include <vector>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
#ifdef WIN32
#include "windows.h"
#endif
#include <GL/gl.h>

namespace QGLVEXT{

  typedef std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > VectorV3d;
  typedef std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f> > VectorV3f;
  
  /**
   * @class DrawWithNames providing drawwithnmaes functions for performing
   * select operation.
   * 
   */
  class DrawWithNames{
	
  public:
	static void draw(const VectorV3f &verts){

	  for (int i = 0; i < (int)verts.size(); ++i){

		glPushName(i);
		glBegin(GL_POINTS);
		glVertex3f(verts[i].x(),verts[i].y(),verts[i].z());
		glEnd();
		glPopName();
	  }
	}

	static void draw(const VectorV3d &verts){

	  for (int i = 0; i < (int)verts.size(); ++i){

		glPushName(i);
		glBegin(GL_POINTS);
		glVertex3d(verts[i].x(),verts[i].y(),verts[i].z());
		glEnd();
		glPopName();
	  }
	}
	
  };
  
}//end of namespace

#endif /*_DRAWWITHNMAES_H_*/
