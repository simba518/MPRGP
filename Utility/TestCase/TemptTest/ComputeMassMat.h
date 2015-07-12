#ifndef _COMPUTEMASSMAT_H_
#define _COMPUTEMASSMAT_H_

#include <TetMesh.h>
#include <eigen3/Eigen/Sparse>
#include <CASADITools.h>
using namespace UTILITY;
using CasADi::SX;
using CasADi::SXMatrix;

namespace ELASTIC_OPT{
  
  /**
   * @class compute lumped mass matrix using casadi.
   * 
   */
  class ComputeMassMat{

  public:
	void compute(SXMatrix &M,const TetMesh&mesh, const bool diagonal=true){
	  vector<SX> density(mesh.material()._rho.size());
	  for (int i = 0;i<mesh.material()._rho.size(); ++i)
		density[i] = mesh.material()._rho[i];
	  if (diagonal){
		computeDiag(M,mesh, density);
	  }else{
		compute(M,mesh, density);
	  }
	}
	void computeDiag(SXMatrix &M,const TetMesh&mesh, const vector<SX> &density);
	void compute(SXMatrix &M,const TetMesh&mesh, const vector<SX> &density);

  protected:
	void assembleMass(const TetMesh&mesh,const Vector4i&tet,const SX&density,SXMatrix&M)const;
	void computeCompactM(SXMatrix &M,const TetMesh& mesh)const;

  private:
	SXMatrix _M; // compact mass matrix.
	vector<SX> _density;
  };

  typedef boost::shared_ptr<ComputeMassMat> pComputeMassMat;
  
}//end of namespace

#endif /*_COMPUTEMASSMAT_H_*/
