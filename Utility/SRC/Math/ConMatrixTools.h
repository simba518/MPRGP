#ifndef _CONMATRIXTOOLS_H_
#define _CONMATRIXTOOLS_H_

#include <vector>
#include <set>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <boost/foreach.hpp>
#include <assertext.h>
using std::vector;
using std::set;
using namespace Eigen;

namespace UTILITY{

  typedef vector<Eigen::Triplet<double> > VecT;
  
  inline void insertI3x3(SparseMatrix<double> &m, const int node_id, const int c_index){
  
	assert_ge(node_id,0);
	assert_ge(m.cols(),3);
	assert_ge(m.rows(),3);
	assert_gt(m.cols(),node_id*3+2);
	assert_gt(m.rows(),c_index*3+2);

	m.insert(c_index*3+0,node_id*3+0) = 1;
	m.insert(c_index*3+1,node_id*3+1) = 1;
	m.insert(c_index*3+2,node_id*3+2) = 1;
  }

  // Compute con matrix fixed all of the provided nodes, not only barycenter.
  template<class VECTOR_I>
  inline void computeConM(const VECTOR_I &con_nodes, SparseMatrix<double> &C, const int total_node_num){
  
	assert_gt(total_node_num,0);
	assert_gt(con_nodes.size(),0);

	C.resize(con_nodes.size()*3, total_node_num*3);
	C.reserve(con_nodes.size()*3);

	vector<int> con_nodes_vec;
	con_nodes_vec.reserve((con_nodes.size()));
	BOOST_FOREACH(int ele, con_nodes)
	  con_nodes_vec.push_back(ele);
	
	for (int i = 0; i < (int)con_nodes_vec.size(); ++i)
	  insertI3x3(C, con_nodes_vec[i], i);
  }

  // Compute con matrix fixed all of the provided nodes, not only barycenter.
  template<class VECTOR_I>
  inline void computeConM(const VECTOR_I &con_nodes, VecT &trip_C, const int total_node_num){
	SparseMatrix<double> mat;
	computeConM(con_nodes, mat, total_node_num);
	trip_C.clear();
	trip_C.reserve(mat.nonZeros());
	for (int k=0; k <mat.outerSize(); ++k)
	  for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
		trip_C.push_back(Triplet<double>(it.row(),it.col(),it.value()));
  }

  // compute the constraint matrix of the barycenter constraints.
  void computeBaryCenterConM(const set<int> &con_nodes,int total_node_num,SparseMatrix<double> &con_m);

  // compute the constraint matrix of the barycenter constraints.
  void computeBaryCenterConM(const vector<set<int> > &con_nodes,int total_node_num,VecT &con_M_triplet);

  // compute the constraint matrix of the barycenter constraints.
  inline void computeBaryCenterConM(const vector<set<int> > &con_nodes,int total_node_num,SparseMatrix<double> &con_m){
	int nz = 0;
	VecT C_triplet;
	computeBaryCenterConM(con_nodes,total_node_num,C_triplet);
	con_m.resize (con_nodes.size()*3,total_node_num*3);
	con_m.reserve (nz);
	con_m.setFromTriplets(C_triplet.begin(), C_triplet.end());
  }

  // compute the reduced constraint matrix of the barycenter constraints.
  inline void computeRedBaryCenterConM(const set<int> &con_nodes,const MatrixXd &U,MatrixXd &con_m){
	SparseMatrix<double> sparse_con_m;
	int total_node_num = U.rows() / 3;
	computeBaryCenterConM(con_nodes,total_node_num,sparse_con_m);
	con_m = sparse_con_m * U;
  }

  // compute the reduced constraint matrix of the barycenter constraints.
  inline void computeRedBaryCenterConM(const vector<set<int> > &con_nodes,const MatrixXd &U, MatrixXd &con_m){
	const int total_node_num = U.rows()/3;
	const int r = U.cols();
	con_m.resize(3*con_nodes.size(),r);
	int i = 0;
	BOOST_FOREACH(const set<int> &group, con_nodes){
	  SparseMatrix<double> sparse_con_m;
	  computeBaryCenterConM(group,total_node_num,sparse_con_m);
	  con_m.block(i*3,0,3,r) = sparse_con_m*U;
	  ++i;
	}
  }

  // compute the reduced constraint matrix of the barycenter constraints.
  inline void computeRedBaryCenterConM(const vector<set<int> > &con_nodes,const MatrixXd &U,vector<MatrixXd> &C){
	C.resize(con_nodes.size());
	vector<set<int> >::const_iterator it = con_nodes.begin();
	for( int i = 0; it != con_nodes.end(); it ++ ,i++ ){
	  computeRedBaryCenterConM(*it,U,C[i]);
	}
  }
}

#endif /* _CONMATRIXTOOLS_H_ */
