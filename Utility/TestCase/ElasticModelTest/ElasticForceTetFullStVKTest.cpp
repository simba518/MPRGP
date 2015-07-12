#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include "../CommonTest/TetMeshFactoryForTest.h"
#include <ElasticForceTetFullStVK.h>
#include <MatrixIO.h>
#include <DefGradTet.h>
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(ElasticForceTetFullStVKTest)

BOOST_AUTO_TEST_CASE(testKTet2){
  
  // initialize
  pTetMesh_const vol_mesh = TetMeshFactoryForTest::tet2();
  ElasticForceTetFullStVK stvk_fun(vol_mesh);
  TEST_ASSERT(stvk_fun.prepare());

  // forces
  VectorXd X0;
  vol_mesh->nodes(X0);
  VectorXd out(X0.size());
  out.setZero();
  stvk_fun.force(X0,out);
  ASSERT_EQ ( out.norm(), 0.0f);

  const VectorXd X = VectorXd::Random(X0.size());
  stvk_fun.force(X,out);
  ASSERT_GT ( out.norm(), 1e-3);

  const VectorXd X1 = X0+VectorXd::Ones(X0.size());
  stvk_fun.force(X1,out);
  ASSERT_EQ ( out.norm(), 0.0f);

  // K(x)
  const SparseMatrix<double> &K = stvk_fun.K(X);
  ASSERT_EQ (K.rows(),K.cols());
  ASSERT_EQ (K.rows(),(int)vol_mesh->nodes().size()*3);
  const SparseMatrix<double> Kt = K.transpose();
  ASSERT_GT ( K.norm(), 1e-3);
  ASSERT_LT ( (K - Kt).norm(), 1e-7*K.norm());

  // check eigen values of K(0)
  const Matrix<double,-1,-1> K0 = stvk_fun.K(X0);
  SelfAdjointEigenSolver<Matrix<double,-1,-1> > eigensolver(K0);
  TEST_ASSERT (eigensolver.info() == Eigen::Success);

  /// @bug the eigenvalues of K(x) are all minus.
  // cout << eigensolver.eigenvalues() << endl;

  // K(x)*dx
  out.setZero();
  const VectorXd dx = VectorXd::Random(vol_mesh->nodes().size()*3);
  stvk_fun.Kdx(dx, X, out);
  ASSERT_GT ( out.norm(), 1e-3);

  const SparseMatrix<double> &K2 = stvk_fun.K(X);
  cout << out.transpose() << endl;
  cout << (K2*dx).transpose() << endl;

  // stvk_fun.force(X,out);
  // cout << out.transpose() << endl;

  ASSERT_LE ((out-(K2*dx)).norm(), 1e-12*out.norm());

}

BOOST_AUTO_TEST_CASE(testK){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"beam.abq";
  pTetMesh mesh = pTetMesh(new TetMesh);
  TEST_ASSERT ( mesh->load(tet_fname) );
  mesh->material().reset(1000.0f,2E6,0.45);
  const int n = mesh->nodes().size();

  const string stiff_mat_f = std::string(TEST_DATA_DIR)+"beam_sparse_K.b";
  SparseMatrix<double> ref_K;
  TEST_ASSERT ( load(ref_K,stiff_mat_f) );
  ASSERT_EQ(ref_K.rows(),n*3);
  ASSERT_EQ(ref_K.cols(),n*3);

  ElasticForceTetFullStVK elas(mesh);
  TEST_ASSERT(elas.prepare());
  VectorXd x0(n*3);
  mesh->nodes(x0);
  const SparseMatrix<double> &K = elas.K(x0);
  
  ASSERT_EQ(K.rows(),n*3);
  ASSERT_EQ(K.cols(),n*3);
  ASSERT_LT((K+ref_K).norm(),1e-5); /// @bug here we use -K.
}

BOOST_AUTO_TEST_CASE(testEnergy){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"beam.abq";
  pTetMesh mesh = pTetMesh(new TetMesh);
  TEST_ASSERT ( mesh->load(tet_fname) );
  mesh->material().reset(1000.0f,2E6,0.45);
  const int n = mesh->nodes().size();

  const string u_str = std::string(TEST_DATA_DIR)+"disp_for_test_elastic_energy.b";
  VectorXd u;
  TEST_ASSERT( load(u_str,u) );
  VectorXd x;
  mesh->nodes(x);
  ASSERT_EQ(x.size(),u.size());
  x += u;

  ElasticForceTetFullStVK elas(mesh);
  TEST_ASSERT(elas.prepare());
  const double energy = elas.energy(x);
  ASSERT_EQ_TOL(energy,11720539881.8642311,1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
