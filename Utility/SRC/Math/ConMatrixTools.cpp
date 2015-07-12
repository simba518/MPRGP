#include <ConMatrixTools.h>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/ref.hpp>

void insert_I3x3(SparseMatrix<double> &m,int node_id){

  assert_ge(node_id,0);
  assert_gt(m.cols(),0);
  assert_gt(m.rows(),0);
  assert_gt(m.cols(),node_id*3+2);
  m.insert(0,node_id*3+0) = 1;
  m.insert(1,node_id*3+1) = 1;
  m.insert(2,node_id*3+2) = 1;
}

void UTILITY::computeBaryCenterConM(const set<int> &con_nodes,int total_node_num,SparseMatrix<double> &con_m){

  assert_gt(total_node_num,0);
  assert_gt(con_nodes.size(),0);
  con_m.resize(3,total_node_num*3);
  for_each(con_nodes.begin(),con_nodes.end(),
  		   boost::lambda::bind(&insert_I3x3,boost::ref(con_m),boost::lambda::_1));
  con_m *= 1.0f/(con_nodes.size());
}

void UTILITY::computeBaryCenterConM(const vector<set<int> > &con_nodes,int total_node_num,VecT &C_triplet){

  int nz = 0;
  for (size_t i = 0; i < con_nodes.size(); ++i){
	nz += con_nodes[i].size()*3;
  }
  C_triplet.reserve(nz);
  int row_0 = 0;
  vector<set<int> >::const_iterator node_it = con_nodes.begin();
  for( int i = 0; node_it != con_nodes.end(); node_it ++ ,i++ ){
	SparseMatrix<double> C;
	computeBaryCenterConM(*node_it,total_node_num,C);
	for (int k=0; k<C.outerSize(); ++k){
	  for (SparseMatrix<double>::InnerIterator it(C,k); it; ++it){
		const int r = it.row()+row_0;
		const int c = it.col();
		assert_in(r,0,(int)con_nodes.size()*3);
		assert_in(c,0,total_node_num*3);
		C_triplet.push_back(Eigen::Triplet<double>(r,c,it.value()));
	  }
	}
	row_0 += C.rows();
  }
  assert_eq((int)con_nodes.size()*3,row_0);
}
