#include "ComputeMassMat.h"

using namespace ELASTIC_OPT;

#define ADD_INFLUENCE_AD(a,b,M)							\
  M.elem((int)tet[(int)a],(int)tet[(int)b])+=influence;	\
  M.elem((int)tet[(int)b],(int)tet[(int)a])+=influence;

void ComputeMassMat::compute(SXMatrix &M,const TetMesh&mesh, const vector<SX> &density){
  
  assert_eq(density.size(),mesh.tets().size());
  this->_density = density;

  _M.setZero();
  computeCompactM(_M,mesh);

  const int n = mesh.nodes().size();
  M.resize(n*3,n*3);
  M.setZero();

  for (int i = 0; i < _M.size1(); ++i){
	for (int j = 0; j < _M.size2(); ++j){
	  if( _M.hasNZ(i,j) ){
		M.elem(i*3+0,j*3+0) = _M.elem(i,j);
		M.elem(i*3+1,j*3+1) = _M.elem(i,j);
		M.elem(i*3+2,j*3+2) = _M.elem(i,j);
	  }
	}
  }
}

void ComputeMassMat::computeDiag(SXMatrix&M,const TetMesh&mesh, const vector<SX>&density){
	  
  assert_eq(density.size(),mesh.tets().size());
  this->_density = density;

  _M.setZero();
  computeCompactM(_M,mesh);

  const int n = mesh.nodes().size();
  M.resize(n*3,n*3);
  M.setZero();

  for (int i = 0; i < _M.size1(); ++i){
	for (int j = 0; j < _M.size2(); ++j){
	  if( _M.hasNZ(i,j) ){
		M.elem(i*3+0,i*3+0) += _M.elem(i,j);
		M.elem(i*3+1,i*3+1) += _M.elem(i,j);
		M.elem(i*3+2,i*3+2) += _M.elem(i,j);
	  }
	}
  }
}

void ComputeMassMat::assembleMass(const TetMesh& mesh,const Vector4i& tet,const SX& density,SXMatrix &M)const{

  SX influence=tetrahedron(mesh.nodes()[tet[0]],mesh.nodes()[tet[1]],
						   mesh.nodes()[tet[2]],mesh.nodes()[tet[3]]).volume()/20.0f;
  influence*=density;

  ADD_INFLUENCE_AD(0,0,M);
  ADD_INFLUENCE_AD(1,1,M);
  ADD_INFLUENCE_AD(2,2,M);
  ADD_INFLUENCE_AD(3,3,M);

  ADD_INFLUENCE_AD(1,0,M);
  ADD_INFLUENCE_AD(2,0,M);
  ADD_INFLUENCE_AD(3,0,M);

  ADD_INFLUENCE_AD(1,2,M);
  ADD_INFLUENCE_AD(2,3,M);
  ADD_INFLUENCE_AD(3,1,M);
}

void ComputeMassMat::computeCompactM(SXMatrix &M,const TetMesh& mesh)const{

  const int n = (int)mesh.nodes().size();
  M.resize(n,n);
  for(int i=0;i<(int)mesh.tets().size();i++)
	assembleMass(mesh,mesh.tets()[i],_density[i],M);
}
