#include <MeshRender.h>
using namespace QGLVEXT;
using namespace UTILITY;
using namespace Eigen;

class ObjRender:public SelfRenderEle{
public:
  ObjRender(const Objmesh &obj):_obj(obj){}
  void draw()const{UTILITY::draw(_obj);}
private:
  const Objmesh _obj;
};

class TetRender:public SelfRenderEle{
public:
  TetRender(const TetMesh &tet,const VectorXd &u):_tet(tet),_u(u){}
  TetRender(const TetMesh &tet):_tet(tet){}
  void move(const VectorXd &u){
	assert_eq(u.size(),_tet.nodes().size()*3);
	_u = u;
  }
  void draw()const{
	if (_u.size() > 0){
	  UTILITY::draw(_tet,&_u[0]);
	}else{
	  UTILITY::draw(_tet);
	}
  }
private:
  const TetMesh _tet;
  VectorXd _u; 
};
