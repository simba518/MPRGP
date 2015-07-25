#ifndef _TETMESH_H_
#define _TETMESH_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <algorithm>
#include <HashedId.h>
#include <assertext.h>
#include <BBox.h>
#include <VTKWriter.h>
#include <AuxTools.h>
using namespace std;
using namespace Eigen;

namespace UTILITY{
 
  template <typename T>
  struct ElasticMaterial{

	void reset(const int elements_num,const T rho,
			   const T E/*Young's*/,
			   const T v/*Poisson*/){
	  assert_ge(elements_num,0);
	  assert_gt(rho,0.0f);
	  assert_gt(E,0.0f);
	  assert_gt(v,0.0f);
	  assert_in(v,0.0f,0.5f);
	  _rho.assign(elements_num,rho);
	  const Matrix<T,2,1> gl = toLameConstant(E,v);
	  _G.assign(elements_num,gl[0]);
	  _lambda.assign(elements_num,gl[1]);
	}
	void reset(const int elements_num){
	  reset(elements_num,1000.0f,2e6,0.45);
	}
	void reset(const T rho,const T E,const T v){
	  const size_t elements_num = _rho.size();
	  reset(elements_num,rho,E,v);
	}
  
	/// yp: [ Young's modulus(E), Poisson's ratio(v)]
	/// gl: [ Shear modulus(G),  Lamé's first parameter(Lambda)]
	/// http://en.wikipedia.org/wiki/Lam%C3%A9%27s_first_parameter
	static Matrix<T,2,1> toLameConstant(const Matrix<T,2,1>& yp){
	  Matrix<T,2,1> gl;
	  gl.x()=yp.x()/(2.0f*(1.0f+yp.y()));
	  gl.y()=yp.x()*yp.y()/((1.0f+yp.y())*(1.0f-2.0f*yp.y()));
	  return gl;
	}
	static Matrix<T,2,1> fromLameConstant(const Matrix<T,2,1>& gl){
	  Matrix<T,2,1> yp;
	  yp.x()=gl.x()*(3.0f*gl.y()+2.0f*gl.x())/(gl.x()+gl.y());
	  yp.y()=gl.y()/(2.0f*(gl.x()+gl.y()));
	  return yp;
	}
	static Matrix<T,2,1> toLameConstant(const T y,const T p){
	  Matrix<T,2,1> yp;  yp[0] = y; yp[1] = p;
	  return toLameConstant(yp);
	}
	static Matrix<T,2,1> fromLameConstant(const T g, const T L){
	  Matrix<T,2,1> gl;  gl[0]=g; gl[1]=L;
	  return fromLameConstant(gl);
	}

	static T getPoissonV(const T&G, const T&L){return L/(2.0*(L+G));}
	static T getYoungE(const T&G, const T&L){return G*(3.0*L+2.0*G)/(L+G);}
	static vector<T> getPoissonV(const vector<T>&G, const vector<T>&L){
	  assert_eq(G.size(),L.size());
	  vector<T> v(G.size());
	  for (size_t i = 0; i < G.size(); ++i)
		v[i] = getPoissonV(G[i],L[i]);
	  return v;
	}
	static vector<T> getYoungE(const vector<T>&G, const vector<T>&L){
	  assert_eq(G.size(),L.size());
	  vector<T> E(L.size());
	  for (size_t i = 0; i < L.size(); ++i)
		E[i] = getYoungE(G[i],L[i]);
	  return E;
	}

	vector<T> getPoissonV()const{
	  return getPoissonV(_G, _lambda);
	}
	vector<T> getYoungE()const{
	  return getYoungE(_G, _lambda);
	}

	std::vector<T> _rho; // Density
	std::vector<T> _G; // Shear modulus.
	std::vector<T> _lambda; // Lamé's first parameter.
  };

  template <typename T>
  class tetrahedronTpl{

  public:
	typedef Matrix<T,2,1> PT2;
	typedef Matrix<T,3,1> PT;
	typedef Matrix<T,4,1> PT4;
	typedef Matrix<T,6,1> PT6;
	typedef Matrix<T,3,3> MAT3;
  public:
	tetrahedronTpl(){}
	tetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d)
	  :_a(a),_b(b),_c(c),_d(d){
	  _swap=false;
	  if(volume() < 0.0f){
		swap(_c,_d);
		_swap=true;
	  }
	}
	PT4 bary(const PT& pt) const{
	  MAT3 A;
	  A.col(0)=_a-_d;
	  A.col(1)=_b-_d;
	  A.col(2)=_c-_d;
	  PT abc=A.inverse()*(pt-_d);
	  return PT4(abc.x(),abc.y(),abc.z(),1.0f-abc.sum());
	}
	bool isInside(const PT& pt) const{
	  PT4 bc=bary(pt);
	  return bc.x() >= 0 && bc.y() >= 0 && bc.z() >= 0 && bc.w() >= 0;
	}
	T volume() const{
	  return (_b-_a).cross(_c-_a).dot(_d-_a)/6.0f;
	}
	const PT& getNode(const int& i) const{return (&_a)[i];}
	const Vector3d center()const{
	  const Vector3d pos = (_a+_b+_c+_d)*0.25;
	  return pos;
	}

  public:
	//data
	PT _a;
	PT _b;
	PT _c;
	PT _d;
	bool _swap;
  };
  typedef tetrahedronTpl<double> tetrahedron;

  typedef std::vector<Vector4i,aligned_allocator<Vector4i> > VVec4i;
  typedef std::vector<Vector3i,aligned_allocator<Vector3i> > VVec3i;
  typedef std::vector<Vector3d,aligned_allocator<Vector3d> > VVec3d;
  typedef boost::unordered_map<HashedId,std::pair<int,int>,HashedIdHash> FaceId;
  typedef std::vector<boost::unordered_set<int> > VectorUseti;
 
  class TetMesh{

  public:
	// init
	TetMesh(){}
	TetMesh(const VVec3d& nodes, const VVec4i& tets){
	  reset(nodes,tets);
	}
	void reset(const VVec3d& nodes, const VVec4i& tets);
	template<class VECTOR_DOUBLE,class VECTOR_INT>
	void reset(const VECTOR_DOUBLE& nodes, const VECTOR_INT& tets){
	  assert_eq(nodes.size()%3,0);
	  assert_eq(tets.size()%4,0);
	  VVec3d nodesV(nodes.size()/3);
	  VVec4i tetsV(tets.size()/4);
	  for (size_t i = 0; i < nodesV.size(); ++i){
		nodesV[i][0] = nodes[i*3+0];
		nodesV[i][1] = nodes[i*3+1];
		nodesV[i][2] = nodes[i*3+2];
	  }
	  for (size_t i = 0; i < tetsV.size(); ++i){
		tetsV[i][0] = tets[i*4+0];
		tetsV[i][1] = tets[i*4+1];
		tetsV[i][2] = tets[i*4+2];
		tetsV[i][3] = tets[i*4+3];
	  }
	  reset(nodesV,tetsV);
	}
	template<class VECTOR>
	void applyDeformation(const VECTOR &u){
	  assert_eq(u.size(),_nodes.size()*3);
	  for (size_t i = 0; i < _nodes.size(); ++i)
		_nodes[i] += u.segment(i*3,3);
	}

	// set 
	void setSingleMaterial(const double&dens,const double&E,const double&v){
	  _mtl.reset(_tets.size(),dens,E,v);
	}
	void setRestPos(const VVec3d& pos){
	  _nodes = pos;
	}

	// get 
	const VVec3d& nodes() const{return _nodes;}
	template<class VECTOR>
	void nodes(VECTOR &x) const{
	  x.resize(_nodes.size()*3);
	  for (size_t i = 0; i < _nodes.size(); ++i){
		x[i*3+0] = _nodes[i][0];
		x[i*3+1] = _nodes[i][1];
		x[i*3+2] = _nodes[i][2];
	  }
	}
	const VVec4i& tets() const{return _tets;}
	template<class VECTOR>
	void tets(VECTOR &x) const{
	  x.resize(_tets.size()*4);
	  for (size_t i = 0; i < _tets.size(); ++i){
		x[i*4+0] = _tets[i][0];
		x[i*4+1] = _tets[i][1];
		x[i*4+2] = _tets[i][2];
		x[i*4+3] = _tets[i][3];
	  }
	}
	const FaceId& faceId() const{return _faceId;}
	const VVec3i& surface() const{return _surface;}
	template<class VECTOR_I>
	void surface(VECTOR_I &f) const{
	  f.resize(_surface.size()*3);
	  for (size_t i = 0; i < _surface.size(); ++i){
		f[i*3+0] = _surface[i][0];
		f[i*3+1] = _surface[i][1];
		f[i*3+2] = _surface[i][2];
	  }
	}
	const VVec3d& normal() const{return _normal;}
	const ElasticMaterial<double>& material() const{return _mtl;}
	ElasticMaterial<double>& material() {return _mtl;}
	double volume(const int i) const{
	  return tetrahedron(_nodes[_tets[i][0]],
	  					 _nodes[_tets[i][1]],
	  					 _nodes[_tets[i][2]],
	  					 _nodes[_tets[i][3]]).volume();
	}
	const VVec4i &faceNeighTet()const{
	  return _faceNeighTet;
	}
	const VectorUseti &nodeNeighNode()const{
	  return _nodeNeighNode;
	}
	const Vector3d &node(const int ele, const int i)const{
	  assert_in(i,0,3);
	  assert_in(ele,0,(int)_tets.size()-1);
	  assert_in(_tets[ele][i],0,(int)_nodes.size()-1);
	  return _nodes[_tets[ele][i]];
	}
	const tetrahedron getTet(const int ele)const{
	  assert_in(ele,0,(int)_tets.size()-1);
	  return tetrahedron(nodes()[_tets[ele][0]],nodes()[_tets[ele][1]],
						 nodes()[_tets[ele][2]],nodes()[_tets[ele][3]]);
	}
	int getContainingElement(const Vector3d &pos)const;
	int getClosestElement(const Vector3d &pos)const;

	int buildInterpWeights(const VectorXd &vertices,vector<int> &nodes,
						   VectorXd&weights,const double zeroThreshold=0.0f)const;
	static void interpolate(const vector<int> &tetNodes,const VectorXd &weights,
							const VectorXd& u,VectorXd& uTarget);
	static void buildInterpMatrix(const vector<int> &tetNodes,const VectorXd &weights,
								  const int tetNodesNum, SparseMatrix<double>& A);

	BBoxD getBBox()const;

	template<class VECTOR3D, class VECTOR>
	void computeGravity(const VECTOR3D &g, VECTOR &gravity)const{

	  gravity.resize(nodes().size()*3);
	  if (gravity.size() > 0){
		memset(&gravity[0], 0, sizeof(gravity[0])*nodes().size()*3);
	  }
	  const double invNumEleVert = 1.0/4.0;
	  for(int el=0; el < _tets.size(); el++){
		const double volume = tetrahedron(nodes()[_tets[el][0]],nodes()[_tets[el][1]],
										  nodes()[_tets[el][2]],nodes()[_tets[el][3]]).volume();
		const double mass_x_invN = _mtl._rho[el]*volume*invNumEleVert;
		for(int j=0; j<4; j++){
		  const int vj = _tets[el][j];
		  gravity[vj*3+0] += mass_x_invN*g[0];
		  gravity[vj*3+1] += mass_x_invN*g[1];
		  gravity[vj*3+2] += mass_x_invN*g[2];
		}
	  }
	}

	// io
	bool load(const std::string& filename);
	bool loadElasticMtl(const std::string& filename);
	bool writeElasticMtl(const std::string& filename)const;

	bool write(const std::string& filename)const;
	template<class VECTOR>
	bool write(const std::string& filename,const VECTOR &u)const{
	  const VVec3d nodes = _nodes;
	  const_cast<TetMesh*>(this)->applyDeformation(u);
	  bool succ = write(filename);
	  const_cast<TetMesh*>(this)->_nodes = nodes;
	  return succ;
	}
	bool write(const std::string& filename,const MatrixXd &U)const;
	template<class VECTOR>
	bool write(const std::string& filename,const vector<VECTOR> &U)const{
	  bool succ = true;
	  for (int i = 0; i < U.size() && succ; ++i)
		succ = write(filename+TOSTR(i)+".vtk",U[i]);
	  return succ;
	}
	bool writeVTK(const std::string& filename,const bool binary=true)const;
	template<class VECTOR>
	bool writeVTK(const std::string& filename,const VECTOR &u,
				  const bool binary=true)const{
	  const VVec3d nodes = _nodes;
	  const_cast<TetMesh*>(this)->applyDeformation(u);
	  bool succ = writeVTK(filename,binary);
	  const_cast<TetMesh*>(this)->_nodes = nodes;
	  return succ;
	}
	bool writeVTK(const std::string& filename,const MatrixXd &U,
				  const bool binary=true)const;
	template<class VECTOR>
	bool writeVTK(const std::string& filename,const vector<VECTOR> &U,
				  const bool binary=true)const{
	  bool succ = true;
	  for (int i = 0; i < U.size() && succ; ++i)
		succ = writeVTK(filename+TOSTR(i)+".vtk",U[i],binary);
	  return succ;
	}
	bool writeElasticMtlVTK(const string filename)const;

  protected:
	void computeSurfaceNormal();

  private:
	VVec3d _nodes;
	VVec4i _tets;
	ElasticMaterial<double> _mtl;

	FaceId _faceId; // 3 vertices of each face, and 2 tets it belongs to.
	VVec4i _faceNeighTet; // 4 neighbor tets(same face) of one tet.
	VectorUseti _nodeNeighNode; // the neighbor nodes of each node.
	VVec3i _surface; // surface's faces
	VVec3d _normal; // surface's faces' normal
  };

  typedef boost::shared_ptr<TetMesh> pTetMesh;
  typedef boost::shared_ptr<const TetMesh> pTetMesh_const;
}

#endif /* _TETMESH_H_ */
