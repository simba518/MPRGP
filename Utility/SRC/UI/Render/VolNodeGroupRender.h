#ifndef _VOLNODEGROUPRENDER_H_
#define _VOLNODEGROUPRENDER_H_

#include <set>
#include <vector>
using namespace std;

#include <assertext.h>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <TetMesh.h>

namespace UTILITY{
  
  enum NODE_SHAPE {DRAW_POINT, DRAW_SHPERE};

  /**
   * @class VolNodeGroupRender draw verterxes in a group of a volume meshes.
   * 
   */
  class VolNodeGroupRender{
	
  public:
	VolNodeGroupRender();
	VolNodeGroupRender(const double color [4],double size);
	void setColor(const double color [4]);
	void setPointSize(const double size){
	  assert_gt (size , 0);
	  this->size = size;
	}
	void setSphereRadius(const double r){
	  assert_gt (r , 0);
	  this->radius = r;
	}

	// draw several groups of vertice, T could be vector<int>, set<int>.
	template<typename T>
	void draw(pTetMesh_const vol_mesh, const vector<T >&groups, const double *vol_u, const NODE_SHAPE s=DRAW_POINT)const{
	  BOOST_FOREACH(const T &ele, groups){
		draw(vol_mesh,ele,vol_u,s);
	  }
	}

	// draw a group of vertice, T could be vector<int>, set<int> or int.
	template<typename T> 
	void draw(pTetMesh_const vol_mesh, const T &group, const double *vol_u, const NODE_SHAPE s=DRAW_POINT)const{

	  if(vol_mesh == NULL){
		return;
	  }

	  switch (s){
	
	  case DRAW_POINT: drawPoints(vol_mesh, group, vol_u);
		break;
	
	  case DRAW_SHPERE: drawShpere(vol_mesh, group, vol_u);
		break;
	
	  default: drawPoints(vol_mesh, group, vol_u);
		break;
	  }
	}

  protected:
	void drawPoints(pTetMesh_const vol_mesh, int vertex_id, const double *vol_u)const;
	void drawShpere(pTetMesh_const vol_mesh, int vertex_id, const double *vol_u)const;

	void drawPoints(pTetMesh_const vol_mesh, const set<int> &group, const double *vol_u)const;
	void drawShpere(pTetMesh_const vol_mesh, const set<int> &group, const double *vol_u)const;

	void drawPoints(pTetMesh_const vol_mesh, const vector<int> &group, const double *vol_u)const;
	void drawShpere(pTetMesh_const vol_mesh, const vector<int> &group, const double *vol_u)const;

  protected:
	void DrawSphere(float x, float y, float z, float fRadius, int M, int N)const;

  private:
	double color[4]; //color of the nodes.
	double size; //size of the points or radius of the sphere.
	double radius;
  };
  
  typedef boost::shared_ptr<VolNodeGroupRender> pVolNodeGroupRender;
  typedef boost::shared_ptr<const VolNodeGroupRender> pVolNodeGroupRender_const;
  
}//end of namespace

#endif /*_VOLNODEGROUPRENDER_H_*/
