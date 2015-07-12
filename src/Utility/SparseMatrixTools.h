#ifndef _SPARSEMATRIXTOOLS_H_
#define _SPARSEMATRIXTOOLS_H_

#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <boost/foreach.hpp>
#include <assertext.h>

namespace EIGEN3EXT{

  /************************************ create ********************************/
  template <class T>  
  const Eigen::SparseMatrix<T> &createFromDense(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &M, 
												Eigen::SparseMatrix<T> &S, const T tol = 1e-16){

	typedef Eigen::Triplet<T> E_Triplet;
	std::vector<E_Triplet> striplet;
	striplet.reserve(M.size());
	for (int i = 0; i < M.rows(); ++i) {
	  for (int j = 0; j < M.cols(); ++j) {
		if ( fabs(M(i,j)) >= tol )
		  striplet.push_back( E_Triplet(i,j,M(i,j)) );
	  }
	}
	S.resize(M.rows(), M.cols());
	S.setFromTriplets(striplet.begin(), striplet.end());
	return S;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> createFromDense(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &M,
											   const T tol = 1e-16){
	
	Eigen::SparseMatrix<T> S;
	return createFromDense(M,S,tol);
  }

  template <class T>
  const Eigen::SparseMatrix<T> &eye(const int n, Eigen::SparseMatrix<T> &A, const T value){

	assert_ge(n,0);
	typedef Eigen::Triplet<T> Tri;
	std::vector<Tri> triplets;
	triplets.reserve(n);
	for (int i = 0; i < n; ++i){
	  triplets.push_back( Tri(i, i, value) );
	}
	A.resize(n,n);
	A.reserve(n);
	A.setFromTriplets( triplets.begin(),triplets.end() );
	A.makeCompressed();
	return A;
  }

  template <class T>
  const Eigen::SparseMatrix<T> eye(const int n, const T value){

	assert_ge(n,0);
	Eigen::SparseMatrix<T> S;
	return eye(n,S,value);
  }

  template <class T, class MatrixType> 
  const Eigen::SparseMatrix<T> &eye(Eigen::SparseMatrix<T> &A_block_diag, const std::vector<MatrixType> &block_mats){

	const size_t mat_num = block_mats.size();
	if (mat_num <=0 ){
	  A_block_diag.resize(0,0);
	  return A_block_diag;
	}

	const size_t n = block_mats[0].size()*mat_num;
	typedef Eigen::Triplet<T> Tri;
	std::vector<Tri> triplets;
	triplets.reserve(n);

	int rows = 0;
	int cols = 0;
	for (size_t mi = 0; mi < mat_num; ++mi){
	  const int r = block_mats[mi].rows();
	  const int c = block_mats[mi].cols();
	  for (int i = 0; i < r; ++i){
		for (int j = 0; j < c; ++j){
		  triplets.push_back( Tri(i+rows, j+cols, block_mats[mi](i,j)) );
		}
	  }
	  rows += r;
	  cols += c;
	}

	A_block_diag.resize(rows,cols);
	A_block_diag.reserve(triplets.size());
	A_block_diag.setFromTriplets( triplets.begin(),triplets.end() );
	A_block_diag.makeCompressed();
	return A_block_diag;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> random(const int r, const int c, const T scalor){

	Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> M = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>::Random(r,c);
	M *= scalor;
	Eigen::SparseMatrix<T> S;
	return createFromDense(M,S);
  }


  
  /************************************ block *********************************/
  template <class T>
  const std::vector<Eigen::Triplet<T> > &getTriplet(const Eigen::SparseMatrix<T>&S,std::vector<Eigen::Triplet<T> >&tri){
	
	tri.reserve(S.nonZeros());
	for(int k=0;k<S.outerSize();++k)
	  for(typename Eigen::SparseMatrix<T>::InnerIterator it(S,k);it;++it)
		tri.push_back(Eigen::Triplet<T>(it.row(),it.col(),it.value()));
	return tri;
  }

  template <class T>
  std::vector<Eigen::Triplet<T> > getTriplet(const Eigen::SparseMatrix<T> &S){
	std::vector<Eigen::Triplet<T> > tri;
	getTriplet(S,tri);
	return tri;
  }

  /** 
   * Generate the sparse matrix P which will remove the rows of Matrix A in
   * remove_rows_set by using:
   * A' = P*A.
   * 
   * @param total_rows rows of the original matrix A.
   * @param remove_rows_set the rows to be remove.
   * @param P output matrix.
   * @param remove if it is false, the rows in remove_rows_set will be
   * preserved, and others will be removed.
   * @note all indices in remove_rows_set should be in [0,total_rows-1].
   * @usage (see MatrixReshape.pdf):
   * remove rows: A = P*A
   * add rows: A = P.t*A
   * remove cols: A = A*P.t
   * add cols: A = A*P
   */	
  template <class T>
  const Eigen::SparseMatrix<T> &genReshapeMatrix(const int total_rows, 
												 const std::set<int> &remove_rows_set, 
												 Eigen::SparseMatrix<T> &P,
												 const bool remove = true){
	  
	typedef Eigen::Triplet<T> Tri;

	const int rows = remove ? total_rows - (int)remove_rows_set.size():(int)remove_rows_set.size();
	const int cols = total_rows;
	const int nonzeros = rows;

	std::vector<Tri> P_triplets;
	P_triplets.reserve(nonzeros);
	  
	if(remove){
	  for (int i = 0; i < total_rows; ++i){
		if ( remove_rows_set.find(i) == remove_rows_set.end() ){
		  P_triplets.push_back( Tri((int)P_triplets.size(), i, 1) );
		}
	  }
	}else{
	  for (int i = 0; i < total_rows; ++i){
		if ( remove_rows_set.find(i) != remove_rows_set.end() ){
		  P_triplets.push_back( Tri((int)P_triplets.size(), i, 1) );
		}
	  }
	}

	P.resize(rows,cols);
	P.reserve( nonzeros );
	P.setFromTriplets( P_triplets.begin(),P_triplets.end() );
	P.makeCompressed();
	return P;
  }

  template <class T>
  Eigen::SparseMatrix<T> genReshapeMatrix(const int total_rows, 
										  const std::set<int> &remove_rows_set,			 
										  const bool remove = true){
	Eigen::SparseMatrix<T> P;
	genReshapeMatrix(total_rows,remove,P,remove);
	return P;
  }

  /** 
   * Generate a P that will remove i-th r sub-rows from Matrix A,
   * A' = P*A.
   * example:
   *     |r1|                                          
   * A = |r2|, r = 2, remove_rows_set = {0}, then A' = |r3|.
   *     |r3|
   */
  template <class T>
  const Eigen::SparseMatrix<T> &genReshapeMatrix(const int total_rows, 
												 const int r,
												 const std::set<int> &remove_rows_set,
												 Eigen::SparseMatrix<T> &P, 
												 const bool remove = true){
	  
	std::set<int> rm_rows_set;
	BOOST_FOREACH(int ele, remove_rows_set){
	  if (ele*r >= 0 && ele*r + r <= total_rows){
		for (int i = 0; i < r; ++i){
		  rm_rows_set.insert(ele*r + i);
		}
	  }
	}
	genReshapeMatrix(total_rows,rm_rows_set,P, remove);
	return P;
  }

  template <class T>
  Eigen::SparseMatrix<T> genReshapeMatrix(const int total_rows, 
										  const int r,
										  const std::set<int> &remove_rows_set,		    
										  const bool remove = true){
	Eigen::SparseMatrix<T> P;
	genReshapeMatrix(total_rows,r,remove_rows_set,P, remove);
	return P;
  }

  template <class T> 
  const Eigen::SparseMatrix<T> block(const Eigen::SparseMatrix<T> &S, 
									 const int r0, const int c0, 
									 const int rows, const int cols){
	return S.block(r0,c0,rows,cols);
  }

  template <class T>
  std::vector<Eigen::Triplet<T> > &addToBlock(std::vector<Eigen::Triplet<T> > &M,
											  const std::vector<Eigen::Triplet<T> > &sub, 
											  const int r0, const int c0){
	assert_ge(r0,0);
	assert_ge(c0,0);
	M.reserve(M.size()+sub.size());
	for (size_t i = 0; i < sub.size(); ++i)
	  M.push_back(Eigen::Triplet<T>(sub[i].row()+r0,sub[i].col()+c0,sub[i].value()));
	return M;
  }

  template <class T>
  std::vector<Eigen::Triplet<T> > &addToBlock(std::vector<Eigen::Triplet<T> > &M,
											  const Eigen::SparseMatrix<T> &sub,
											  const int r0, const int c0){
	std::vector<Eigen::Triplet<T> > tripletList;
	getTriplet(sub,tripletList);
	return addToBlock(M,tripletList,r0,c0);
  }
 
  template <class T>
  Eigen::SparseMatrix<T> &addToBlock(Eigen::SparseMatrix<T> &M,
									 const Eigen::SparseMatrix<T> &sub,
									 const int r0, const int c0){
	assert_le(r0+sub.rows(),M.rows());
	assert_le(c0+sub.cols(),M.cols());
	std::vector<Eigen::Triplet<T> > tripletList;
	getTriplet(M,tripletList);
	addToBlock(tripletList,sub,r0,c0);
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	return M;
  }

  template <class T>
  Eigen::SparseMatrix<T> getLower(const Eigen::SparseMatrix<T> &M){
	
	std::vector<Eigen::Triplet<T> > tri;
	tri.reserve(M.nonZeros());
	for(int k=0;k<M.outerSize();++k)
	  for(typename Eigen::SparseMatrix<T>::InnerIterator it(M,k);it;++it)
		if (it.row() >= it.col())
		  tri.push_back(Eigen::Triplet<T>(it.row(),it.col(),it.value()));
	Eigen::SparseMatrix<T> L(M.rows(),M.cols());
	L.setFromTriplets(tri.begin(),tri.end());
	return L;
  }

  /************************************ inverse *********************************/
  template <class T>
  const Eigen::SparseMatrix<T> &inverse(const Eigen::SparseMatrix<T> &P, Eigen::SparseMatrix<T> &inv_P){

	/// @todo function inverse(), very slow.
	Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> M = P;
	M = M.inverse();
	createFromDense(M, inv_P);
	return inv_P;
  }

  template <class T>
  const Eigen::SparseMatrix<T> inverse(const Eigen::SparseMatrix<T> &P){

	Eigen::SparseMatrix<T> inv_P;
	return inverse(P,inv_P);
  }

}

#endif /* _SPARSEMATRIXTOOLS_H_ */
