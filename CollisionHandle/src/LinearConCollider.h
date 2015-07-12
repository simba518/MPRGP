#ifndef _LINEARCONCOLLIDER_H_
#define _LINEARCONCOLLIDER_H_

#include <vector>
#include <FEMCollider.h>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <MPRGPSolver.h>
#include <FEMMesh.h>
#include <FEMGeom.h>
#include "ClothMesh.h"
#include "ClothCollision.h"
using namespace Eigen;
using namespace std;
USE_PRJ_NAMESPACE

typedef vector<Vector4d,Eigen::aligned_allocator<Vector4d> > VVec4d;
typedef vector<VVec4d > VVVec4d;
typedef vector<Triplet<double,int> > TRIPS;

// n^t*x+p>= 0
struct GeomConCache{
  
  GeomConCache(){
	vert_id = plane_id = -1;
	normal.setZero();
  }

  GeomConCache(const int vi, const int pi, const Vector3d &n):
	vert_id(vi), plane_id(pi), normal(n){
	assert_ge(vi,0);
	assert_ge(pi,0);
	assert_in(n.norm(),(1-1e-8),(1+1e-8));
  }

  bool handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n, VVVec4d &linear_con);

  void getConstraints(TRIPS &trips, vector<double> &rhs, const VVVec4d &linear_con)const;
  
  static int addConPlane(VVec4d &con_planes, const Vector4d &p);
  
  int vert_id;
  int plane_id;
  Vector3d normal;
};

// self collision constraints between a triangle and a vertex.
// c = n^t*(xi-a*xj-b*xk-c*xl), where a,b,c is the weight coordinates.
struct SurfaceSelfConCache{

  SurfaceSelfConCache(){
	i = j = k = l = -1;
	a = b = c = -1.0f;
	n.setZero();
	x0.setZero();
  }

  bool handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV,VVVec4d &linear_con);

  void getConstraints(TRIPS &trips, vector<double> &rhs, const VVVec4d &linear_con)const;

  int i,j,k,l;
  int pi, pj, pk, pl;
  double a,b,c;
  Vector3d n, x0;

protected:
  bool convertToLinearCon(VVVec4d &linear_con);
};

/**
 * Handle both self and geometric collisions, and return them as linear constraints.
 * 1. self collision constraints.
 * 2. geometric collision constraints.
 */
class LinearConCollider:public FEMCollider{
  
public:
  LinearConCollider(const bool decouple_con = false):decouple_constraints(decouple_con){
	reset(0);
  }

  void reset(const size_t num_verts){

	linear_con.clear();
	linear_con.resize(num_verts);
	geom_con.clear();
	surface_self_con.clear();

	coll_as_vert.resize(num_verts);
	coll_as_vert.assign(num_verts, false);

	coll_as_face.resize(num_verts);
	coll_as_face.assign(num_verts, false);
  }

  void setDecoupleConstraints(const bool decouple_con){
	this->decouple_constraints = decouple_con;
  }

  void handle(boost::shared_ptr<FEMBody> b, const Vector4d &plane, const double delta);

  void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV);

  void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);

  const VVVec4d &getLinearCon()const{return linear_con;}

  // get A and c for: A*x >= c
  void getConstraints(SparseMatrix<double> &A, VectorXd &c)const;

  void print()const;

  bool collided(const int vert_id)const{
	assert_in(vert_id, 0, coll_as_vert.size()-1);
	assert_in(vert_id, 0, coll_as_face.size()-1);
	return coll_as_vert[vert_id] || coll_as_face[vert_id];
  }

  void project(Vec3d &v, const int vert_id)const;

protected:
  VVVec4d linear_con;
  vector<GeomConCache> geom_con;
  vector<SurfaceSelfConCache> surface_self_con;
  bool decouple_constraints; // if true, one constraint for each triangle and each vertex.
  vector<bool> coll_as_vert, coll_as_face;
};

class ContinueCollider:public LinearConCollider, public ClothCollision::CollisionHandler{

public:
  ContinueCollider(FEMMesh&mesh,boost::shared_ptr<FEMGeom> geom,
				   const bool decouple_con = false):
	LinearConCollider(decouple_con),femmesh(mesh),geom(geom){}
  void init();
  void collide(const VectorXd &last_pos, const VectorXd &cur_pos);

protected:
  void handle(boost::shared_ptr<ClothMesh::ClothVertex> V1,
			  boost::shared_ptr<ClothMesh::ClothTriangle> T2,
			  const Vec3d n,const Vec4d& omg,scalarD t);
  void handle(boost::shared_ptr<ClothMesh::ClothEdge> E1,
			  boost::shared_ptr<ClothMesh::ClothEdge> E2,
			  const Vec3d n,const Vec4d& omg,scalarD t);

private:
  FEMMesh &femmesh;
  boost::shared_ptr<FEMGeom> geom;
  boost::shared_ptr<ClothMesh> surface, scene;
  ClothCollision CCD;
  vector<pair<int,int> > vol_cloth_map;
};

#endif /* _LINEARCONCOLLIDER_H_ */
