#ifndef _SPARSEGENEIGENSOLVER_H_
#define _SPARSEGENEIGENSOLVER_H_

/**
 * examples for solving eigenvalue problems using arapack can be found at\n
 * http://www.inf-cv.uni-jena.de/proj/linal/eigen_8h-source.html\n purpose:\n
 * solve the general eigen-problem: Ax = lambda*Bx with the n_eig largest
 * eigen-values, where A and B symetric and only the lower part of the matrix
 * are stored.
 */
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <assertext.h>
#include <string>

namespace EIGEN3EXT{

  /**
   *\fn bool eigensym( T *valA,int nnzA,  int *irowA,  int *pcolA,T *valB,int nnzB,  int *irowB,  int *pcolB,int n,int n_eig,T *vals,T *vecs)
   *
   *\brief  solve the general eigen-problem: Ax = lambda*Bx with the n_eig smallest eigen-values.
   *
   *\param[in] valA|valB pointer to an array that stores the nonzero elements of A and B.
   *\param[in] nnzA|nnzB Number of nonzero elements in A and B.
   *\param[in] irowA|irowB  pointer to an array that stores the row indices of the nonzeros in A and B.
   *\param[in] pcolA|pcolB pointer to an array of pointers to the beginning of each column of A (B) in valA (valB).
   *\param[in] n Dimension of the problem.
   *\param[in] n_eig number of eigenvectors to solve
   *
   *\param[out] vals  return the eigen values
   *\param[out] vecs  return the eigen vectors
   *\return success or not,and results is writed to vals and vecs
   *
   *\note
   *(1)the both the sparse matrix A and B should be symetric and only the lower part of the matrix are stored.
   *(2)the memory-space of vals and vecs should be allocated outside(before calling this function).\n
   *(3)onely float and double is supported, because of the limitation of arpack.\n
   *(4)number of the required eigen numbers should small than the dimension of the problem:\n
   *n_eig < n \n
   *\see eigensymTest.cpp
   * 
   * @see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html for mode and which
   * @note By using shift mode 'S', we need use "LM" instead of "SM".
   * some possible pair for (mode, which):
   * ('S', "LM"): smallest eigenvalues. (we can use ('R', "SM"), but it is slow.)
   * ('R', "LM"): largest eigenvalues.
   * ('S', "BE") or ('R','BE'): smallest and largest eigenvalues.
   */
  template <typename T>
  bool eigensym( const T *valA,int nnzA,  const int *irowA,  const int *pcolA,
				 const T *valB,int nnzB,  const int *irowB,  const int *pcolB,
				 int n,int n_eig,
				 T *vals,T *vecs,
				 const char mode = 'S',
				 const std::string which = "LM");
  
  /**
   * @class EigenSparseGenEigenSolver providing the interface for solving the
   * general eigen-value problem for the sparsematrix of the eigen3 library.
   * 
   * @see SparseGenEigenSolver
   * 
   * @note 
   * (1) only the lower part of K and M should be provided, otherwise there
   * will be error like "the coefffients in the matrix is not consistent" when
   * calling the arpack.
   * (2) the matrix K and M should be column compressed, which means that the 
   * function K.makeCompressed() and M.makeCompressed() shoulbe be called after 
   * initialization.
   * for mode and which:
   * @see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html 
   * some possible pair for (mode, which):
   * ('S', "LM"): smallest eigenvalues. (we can use ('R', "SM"), but it is slow.)
   * ('R', "LM"): largest eigenvalues.
   * ('S', "BE") or ('R','BE'): smallest and largest eigenvalues.
   */
  class EigenSparseGenEigenSolver{

  public:
	template <typename real>
	static bool solve (const Eigen::SparseMatrix<real> &K,
					   const Eigen::SparseMatrix<real> &M,
					   Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> &eig_vec,
					   Eigen::Matrix<real,Eigen::Dynamic,1> &eig_val,
					   const int max_eig_num,
					   const char mode = 'S',
					   const std::string which = "LM"){

	  //check the parameters
	  assert_gt(max_eig_num,0);
	  assert_le(max_eig_num,K.rows());
	  assert_gt(K.nonZeros(),0);
	  assert_gt(M.nonZeros(),0);
	  assert_eq(K.rows(),M.rows());
	  assert_eq(K.cols(),M.cols());
	  assert_eq(K.rows(),K.cols());

	  //allocate memory
	  eig_vec.resize(K.rows(),max_eig_num);
	  eig_val.resize(max_eig_num);
  
	  //caculate 
	  const int k_nonzero_num = K.nonZeros();
	  const int k_rows = K.rows();
	  const real *k_data = K.valuePtr();
	  const int *k_rowind = K.innerIndexPtr();
	  const int *k_colptr = K.outerIndexPtr();

	  const int m_nonzero_num = M.nonZeros();
	  const real *m_data = M.valuePtr();
	  const int *m_rowind = M.innerIndexPtr();
	  const int *m_colptr = M.outerIndexPtr();

	  real *p_eig_val = &(eig_val[0]);
	  real *p_eig_vec = &(eig_vec(0,0));

	  bool succ = eigensym(k_data,k_nonzero_num,k_rowind,k_colptr,
	  					   m_data,m_nonzero_num,m_rowind,m_colptr,
	  					   k_rows,max_eig_num,
	  					   p_eig_val,p_eig_vec,mode,which);

	  return  succ;
	}

	template <typename real>
	static bool solve (const Eigen::SparseMatrix<real> &K,
					   const Eigen::DiagonalMatrix<real,-1> &diagM,
					   Eigen::Matrix<real,Eigen::Dynamic,Eigen::Dynamic> &eig_vec,
					   Eigen::Matrix<real,Eigen::Dynamic,1> &eig_val,
					   const int max_eig_num,
					   const char mode = 'S',
					   const std::string which = "LM"){
	  Eigen::SparseMatrix<real> M(diagM.rows(),diagM.cols());
	  typedef Eigen::Triplet<real> T;
	  std::vector<T> triplets;
	  triplets.reserve(diagM.diagonal().rows());
	  for (size_t i = 0; i < diagM.diagonal().rows(); ++i)
		triplets.push_back(T(i,i,diagM.diagonal()[i]));
	  M.setFromTriplets(triplets.begin(),triplets.end());
	  return solve(K,M,eig_vec,eig_val,max_eig_num,mode,which);
	}
	
  };
    
}//end of namespace

#endif /* _SPARSEGENEIGENSOLVER_H_ */
