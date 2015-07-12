#ifndef _DRAGGROUPSELDRAW_H_
#define _DRAGGROUPSELDRAW_H_

#include <assertext.h>
#include <TetMesh.h>
#include <vector>
#include <set>
using namespace std;

namespace UTILITY{

  // @class DragGroupSelDraw draw selected nodes of volume mesh.
  class DragGroupSelDraw{
	
  public:
	static void drawAllGroupsWithPoints(pTetMesh_const vol_mesh,
										const vector<set<int> >&drag_groups, 
										const double*vol_u);

	static void drawAllGroupsWithShpere(pTetMesh_const vol_mesh,
										const vector<set<int> >&drag_groups, 
										const double *vol_u, 
										const double radius);

  protected:
	static void drawPoints(pTetMesh_const vol_mesh,
						   const set<int> &group,
						   const double *vol_u);

	static void drawSphere(pTetMesh_const vol_mesh,
						   const set<int> &group,
						   const double *vol_u, const double radius);

	static void DrawSphere(float x, float y, float z, float fRadius, int M, int N);

  };
  
}//end of namespace

#endif /*_DRAGGROUPSELDRAW_H_*/
