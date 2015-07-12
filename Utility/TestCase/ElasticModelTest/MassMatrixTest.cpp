#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MatrixIO.h>
#include <MassMatrix.h>
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

BOOST_AUTO_TEST_SUITE(MassMatrixTest)

BOOST_AUTO_TEST_CASE(TestCompMassMat){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"beam.abq";
  TetMesh mesh;
  TEST_ASSERT ( mesh.load(tet_fname) );
  mesh.material().reset(1000.0f,2E6,0.45);
  const int n = mesh.nodes().size();

  const string mass_mat_f = std::string(TEST_DATA_DIR)+"beam_sparse_M.b";
  SparseMatrix<double> ref_M;
  TEST_ASSERT ( load(ref_M,mass_mat_f) );
  ASSERT_EQ(ref_M.rows(),n*3);

  // test lumped matrix
  MassMatrix mass;
  DiagonalMatrix<double,-1> M;
  mass.compute(M,mesh);

  DiagonalMatrix<double,-1> lumped_ref_M;
  MassMatrix::lump(ref_M,lumped_ref_M);

  ASSERT_EQ(M.rows(),n*3);
  ASSERT_EQ(M.rows(),lumped_ref_M.rows());
  const MatrixXd diffM = lumped_ref_M.diagonal()-(M.diagonal());
  ASSERT_LE(diffM.norm(),1e-10);

  // test unlumped matrix
  SparseMatrix<double> fullM;
  mass.compute(fullM,mesh);
  ASSERT_EQ(fullM.rows(),n*3);
  ASSERT_EQ(fullM.cols(),n*3);
  SparseMatrix<double> fullM_t = fullM.transpose();
  ASSERT_LE((fullM_t-fullM).norm(),1e-10);
  ASSERT_LE((fullM-ref_M).norm(),1e-10);

}

BOOST_AUTO_TEST_SUITE_END()

// *ELASTIC
// 2E6, 0.45
// *DENSITY
// 1000.000000000000000
