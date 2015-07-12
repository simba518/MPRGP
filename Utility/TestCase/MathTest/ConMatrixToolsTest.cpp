#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ConMatrixTools.h>
using namespace Eigen;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(ConMatrixToolsTest)

BOOST_AUTO_TEST_CASE(insertI3x3Test){
  
  SparseMatrix<double> C(9,12);
  insertI3x3(C,0,0);
  insertI3x3(C,3,2);

  SparseMatrix<double> corect_C(9,12);
  corect_C.insert(0,0) = 1;
  corect_C.insert(1,1) = 1;
  corect_C.insert(2,2) = 1;

  corect_C.insert(6,9) = 1;
  corect_C.insert(7,10) = 1;
  corect_C.insert(8,11) = 1;

  const MatrixXd dC = C;
  const MatrixXd dcorect_C = corect_C;
  ASSERT_EQ_SMALL_MAT(dC,dcorect_C);
}

BOOST_AUTO_TEST_CASE(computeConMTest){
	  
  vector<int> con_nodes;
  con_nodes.push_back(0);
  con_nodes.push_back(1);
  con_nodes.push_back(3);

  SparseMatrix<double> C;
  int total_node_num = 4;
  computeConM(con_nodes,C,total_node_num);
  
  SparseMatrix<double> corect_C(9,12);
  corect_C.insert(0,0) = 1;
  corect_C.insert(1,1) = 1;
  corect_C.insert(2,2) = 1;

  corect_C.insert(3,3) = 1;
  corect_C.insert(4,4) = 1;
  corect_C.insert(5,5) = 1;

  corect_C.insert(6,9) = 1;
  corect_C.insert(7,10) = 1;
  corect_C.insert(8,11) = 1;

  const MatrixXd dC = C;
  const MatrixXd dcorect_C = corect_C;
  ASSERT_EQ_SMALL_MAT(dC,dcorect_C);
}

BOOST_AUTO_TEST_CASE(computeBaryCenterConMTest){
  
  int total_node_num = 4;
  set<int> con_nodes;
  SparseMatrix<double> con_m,correct_con_m;
  con_nodes.insert(0);
  con_nodes.insert(3);
  
  computeBaryCenterConM(con_nodes,total_node_num,con_m);

  correct_con_m.resize(3,total_node_num*3);
  correct_con_m.insert(0,0) = 1.0f;
  correct_con_m.insert(1,1) = 1.0f;
  correct_con_m.insert(2,2) = 1.0f;
  correct_con_m.insert(0,0+9) = 1.0f;
  correct_con_m.insert(1,1+9) = 1.0f;
  correct_con_m.insert(2,2+9) = 1.0f;
  correct_con_m *= 0.5f;

  ASSERT_EQ(correct_con_m.rows(),con_m.rows());
  ASSERT_EQ(correct_con_m.cols(),con_m.cols());
  ASSERT_LE((correct_con_m-con_m).norm(),1e-6);

  vector<set<int> > con_nodes_vec;
  con_nodes_vec.push_back(con_nodes);
  computeBaryCenterConM(con_nodes_vec,total_node_num,con_m);
  ASSERT_EQ(correct_con_m.rows(),con_m.rows());
  ASSERT_EQ(correct_con_m.cols(),con_m.cols());
  ASSERT_LE((correct_con_m-con_m).norm(),1e-6);

  con_nodes_vec.push_back(con_nodes);
  computeBaryCenterConM(con_nodes_vec,total_node_num,con_m);
  ASSERT_EQ(correct_con_m.rows()*2,con_m.rows());
  ASSERT_EQ(correct_con_m.cols(),con_m.cols());

  const MatrixXd m1 = con_m;
  const MatrixXd m2 = correct_con_m;
  const MatrixXd m3 = m1.topLeftCorner(m1.rows()/2,m1.cols());
  const MatrixXd m4 = m1.bottomLeftCorner(m1.rows()/2,m1.cols());
  ASSERT_EQ_SMALL_MAT(m2,m4);
  ASSERT_EQ_SMALL_MAT(m3,m4);
}

BOOST_AUTO_TEST_CASE(computeRedBaryCenterConMTest){
  
  vector<set<int> > con_nodes_group;

  int total_node_num = 4;
  set<int> con_nodes;
  con_nodes.insert(2);
  con_nodes_group.push_back(con_nodes);// (2)

  con_nodes.clear();
  con_nodes.insert(0);
  con_nodes.insert(3);
  con_nodes_group.push_back(con_nodes);// (0,3)
  
  MatrixXd U(12,4); // basis matrix.
  for(int i=0; i < U.rows(); i++){
	for(int j = 0; j < U.cols(); j++)
	  U(i,j) = i * U.cols() + j;
  }

  vector<MatrixXd> C;
  computeRedBaryCenterConM(con_nodes_group,U,C);

  SparseMatrix<double> sparse_con_m;
  sparse_con_m.resize(3,total_node_num*3);
  sparse_con_m.insert(0,0) = 1.0f;
  sparse_con_m.insert(1,1) = 1.0f;
  sparse_con_m.insert(2,2) = 1.0f;
  sparse_con_m.insert(0,0+9) = 1.0f;
  sparse_con_m.insert(1,1+9) = 1.0f;
  sparse_con_m.insert(2,2+9) = 1.0f;
  sparse_con_m *= 0.5f;
  const MatrixXd correct_con_m = sparse_con_m * U;
  ASSERT_EQ(C[1],correct_con_m);

  SparseMatrix<double> sparse_con_m_1;
  sparse_con_m_1.resize(3,total_node_num*3);
  sparse_con_m_1.insert(0,0+6) = 1.0f;
  sparse_con_m_1.insert(1,1+6) = 1.0f;
  sparse_con_m_1.insert(2,2+6) = 1.0f;
  const MatrixXd correct_con_m2 = sparse_con_m_1 * U;

  ASSERT_EQ(C[0],correct_con_m2);
}

BOOST_AUTO_TEST_SUITE_END()
