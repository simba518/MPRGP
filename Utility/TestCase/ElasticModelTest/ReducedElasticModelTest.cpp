#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <MatrixIO.h>
#include <eigen3/Eigen/Dense>
#include <ReducedElasticModel.h>
using namespace Eigen;
using namespace EIGEN3EXT;
using namespace SIMULATOR;

BOOST_AUTO_TEST_SUITE(ReducedElasticModelTest)

BOOST_AUTO_TEST_CASE(testCubaturedElasticModel_all_tets){

  // load data: tetmesh, elastic mtl, and B.
  const string tet_file = std::string(TEST_DATA_DIR)+"dino_scaled.abq";
  pTetMesh tet_mesh = pTetMesh(new TetMesh());
  TEST_ASSERT(tet_mesh->load(tet_file));
  const string elas_file = std::string(TEST_DATA_DIR)+"dino_scaled.elastic";
  TEST_ASSERT(tet_mesh->loadElasticMtl(elas_file));

  const string basis_file = std::string(TEST_DATA_DIR)+"dino_scaled_extended_basis.b";
  MatrixXd B;
  TEST_ASSERT(EIGEN3EXT::load(basis_file, B));
  ASSERT_GT(B.rows(),0);
  ASSERT_GT(B.cols(),0);
  ASSERT_EQ(B.rows(), tet_mesh->nodes().size()*3);

  // init elastic model
  CubaturedElasticModel cub_model(tet_mesh);
  DirectReductionElasticModel direct_model(tet_mesh);
  cub_model.setModalBasis(B);
  direct_model.setModalBasis(B);
  TEST_ASSERT( cub_model.prepare() );
  TEST_ASSERT( direct_model.prepare() );

  // check f(q), K(q)
  const int r = cub_model.reducedDim();
  const VectorXd q = VectorXd::Random(r);

  VectorXd cub_f, direct_f;
  cub_model.evaluateF(q, cub_f);
  direct_model.evaluateF(q, direct_f);
  ASSERT_GT(cub_f.norm(), 1e3);
  ASSERT_EQ_TOL( (cub_f - direct_f).norm(), 0.0f, 1e-12*cub_f.norm());

  MatrixXd cub_K, direct_K;
  cub_model.evaluateK(q, cub_K);
  direct_model.evaluateK(q, direct_K);
  ASSERT_GT(cub_K.norm(), 1e3);
  ASSERT_EQ_TOL( (cub_K - direct_K).norm(), 0.0f, 1e-12*cub_K.norm());

  // check f(0)
  const VectorXd q0 = VectorXd::Random(r)*0.0f;
  cub_model.evaluateF(q0, cub_f);
  direct_model.evaluateF(q0, direct_f);
  ASSERT_EQ_TOL( cub_f.norm(), 0.0f, 1e-11);
  ASSERT_EQ_TOL( direct_f.norm(), 0.0f, 1e-11);
}

BOOST_AUTO_TEST_CASE(testCubaturedElasticModel_cub_tets){

  // load data: tetmesh, elastic mtl, and B.
  const string tet_file = std::string(TEST_DATA_DIR)+"dino_scaled.abq";
  pTetMesh tet_mesh = pTetMesh(new TetMesh());
  TEST_ASSERT(tet_mesh->load(tet_file));
  const string elas_file = std::string(TEST_DATA_DIR)+"dino_scaled.elastic";
  TEST_ASSERT(tet_mesh->loadElasticMtl(elas_file));

  const string basis_file = std::string(TEST_DATA_DIR)+"dino_scaled_extended_basis.b";
  MatrixXd B;
  TEST_ASSERT(EIGEN3EXT::load(basis_file, B));
  ASSERT_GT(B.rows(),0);
  ASSERT_GT(B.cols(),0);
  ASSERT_EQ(B.rows(), tet_mesh->nodes().size()*3);

  // init elastic model
  CubaturedElasticModel cub_model(tet_mesh);
  vector<int> cub_tets;
  vector<double> weights;
  cub_tets.push_back(1);
  cub_tets.push_back(10);
  cub_tets.push_back(100);
  cub_tets.push_back(1000);
  weights.push_back(2.0);
  weights.push_back(3.0);
  weights.push_back(1.0);
  weights.push_back(4.0);
  cub_model.setCubature(weights, cub_tets);

  cub_model.setModalBasis(B);
  TEST_ASSERT( cub_model.prepare() );

  // check f(0)
  VectorXd cub_f, direct_f;
  const int r = cub_model.reducedDim();
  const VectorXd q0 = VectorXd::Random(r)*0.0f;
  cub_model.evaluateF(q0, cub_f);
  ASSERT_EQ_TOL( cub_f.norm(), 0.0f, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
