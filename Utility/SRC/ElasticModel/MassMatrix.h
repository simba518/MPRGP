#ifndef _MASSMATRIX_H_
#define _MASSMATRIX_H_

#include <boost/shared_ptr.hpp>
#include <TetMesh.h>
#include <eigen3/Eigen/Sparse>

namespace UTILITY{

#define ADD_INFLUENCE(a,b,M)								\
  M.coeffRef((int)tet[(int)a],(int)tet[(int)b])+=influence;	\
  M.coeffRef((int)tet[(int)b],(int)tet[(int)a])+=influence;
  
  /// @class MassMatrix calculate the lumped mass matrix for tet mesh.
  class MassMatrix{
	
  public:
	// compute lumped matrix
	void compute(Eigen::DiagonalMatrix<double,-1>&M,const TetMesh&mesh){
	  computeCompactM(_M,mesh);
	  //accumulate entries
	  M.resize(mesh.nodes().size()*3);
	  M.setZero();
	  for(int k=0;k<_M.outerSize();++k)
	  	for(Eigen::SparseMatrix<double>::InnerIterator it(_M,k);it;++it)
	  	  M.diagonal().block<3,1>(it.row()*3,0)+=Vector3d::Constant(it.value());
	}

	// compute unlumped matrix
	void compute(Eigen::SparseMatrix<double>&M,const TetMesh&mesh,const bool lower=false){

	  computeCompactM(_M,mesh);

	  //accumulate entries
	  typedef Eigen::Triplet<double> E_Triplet;
	  std::vector<E_Triplet> tri;
	  for(int k=0;k<_M.outerSize();++k){
		for(Eigen::SparseMatrix<double>::InnerIterator it(_M,k);it;++it){
		  const int r3 = it.row()*3;
		  const int c3 = it.col()*3;
		  const double v = it.value();
		  assert_gt(v,0.0);
		  if (!lower || r3 >= c3){
			tri.push_back(E_Triplet(r3+0,c3+0,v));
			tri.push_back(E_Triplet(r3+1,c3+1,v));
			tri.push_back(E_Triplet(r3+2,c3+2,v));
		  }
		}
	  }

	  const int n3 = mesh.nodes().size()*3;
	  M.resize(n3,n3);
	  M.setFromTriplets(tri.begin(), tri.end());
	}

	static void lump(const Eigen::SparseMatrix<double> &fullM,
					 Eigen::DiagonalMatrix<double,-1>&M/*lumped M*/){
	  const int n3 = fullM.cols();
	  M.resize(n3);
	  M.setZero();
	  for (int i = 0; i < n3; ++i)
		M.diagonal()[i] = fullM.col(i).sum();
	}

  protected:
	void assembleMass(const TetMesh& mesh,const Vector4i& tet,const double& density,
					  Eigen::SparseMatrix<double> &M)const{

	  double influence=tetrahedron(mesh.nodes()[tet[0]],mesh.nodes()[tet[1]],
								   mesh.nodes()[tet[2]],mesh.nodes()[tet[3]]).volume()/20.0f;
	  influence*=density;

	  ADD_INFLUENCE(0,0,M);
	  ADD_INFLUENCE(1,1,M);
	  ADD_INFLUENCE(2,2,M);
	  ADD_INFLUENCE(3,3,M);

	  ADD_INFLUENCE(1,0,M);
	  ADD_INFLUENCE(2,0,M);
	  ADD_INFLUENCE(3,0,M);

	  ADD_INFLUENCE(1,2,M);
	  ADD_INFLUENCE(2,3,M);
	  ADD_INFLUENCE(3,1,M);
	}
	void computeCompactM(Eigen::SparseMatrix<double> &M,const TetMesh& mesh)const{
	  M.resize((int)mesh.nodes().size(),(int)mesh.nodes().size());
	  for(int i=0;i<(int)mesh.tets().size();i++)
	  	assembleMass(mesh,mesh.tets()[i],mesh.material()._rho[i],M);
	}

  private:
	Eigen::SparseMatrix<double> _M; // compact mass matrix.
  };
  
  typedef boost::shared_ptr<MassMatrix> pMassMatrix;
  
}//end of namespace

#endif /* _MASSMATRIX_H_ */
