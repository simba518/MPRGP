
void testK(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"beam.abq";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  mesh->material().reset(1000.0f,2E6,0.45);
  const int n = mesh->nodes().size();

  const string stiff_mat_f = std::string(TEST_DATA_DIR)+"beam_sparse_K.b";
  SparseMatrix<double> ref_K;
  assert ( load(ref_K,stiff_mat_f) );
  assert_eq(ref_K.rows(),n*3);
  assert_eq(ref_K.cols(),n*3);

  ComputeStiffnessMat elas(mesh);
  assert(elas.prepare());
  VectorXd x0(n*3);
  mesh->nodes(x0);
  const SXMatrix &K = elas.K(x0);
  
  MatrixXd Km;
  CASADI::convert(K,Km);
  const MatrixXd ref_Km = ref_K;
  
  assert_eq(Km.rows(),n*3);
  assert_eq(Km.cols(),n*3);
  assert_lt((Km+ref_Km).norm(),1e-5);
  
  cout<< "(K-K_ref).norm(): " << (Km+ref_Km).norm() << endl;
  cout << "test end.." << endl;
}

void test_diagM(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"two_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"two_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  MassMatrix mass_c;
  Eigen::DiagonalMatrix<double,-1> diagM;
  mass_c.compute(diagM, *mesh);

  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M, *mesh);
  MatrixXd Me;
  CASADI::convert(M,Me);

  assert_eq(Me.rows(), Me.cols());
  assert_eq(Me.rows(), diagM.rows());
  assert_eq(Me.rows(), diagM.cols());
  const MatrixXd mm = Me-MatrixXd(diagM);
  assert_lt(mm.norm(),1e-8);
  cout<< "(M-Mref).norm() = " << mm.norm() << endl;
}

void testM(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"two_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"two_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  MassMatrix mass_c;
  SparseMatrix<double> Me;
  mass_c.compute(Me, *mesh);

  ComputeMassMat mass;
  SXMatrix M;
  mass.compute(M, *mesh, false);
  MatrixXd Ms;
  CASADI::convert(M,Ms);

  assert_eq(Me.rows(), Me.cols());
  assert_eq(Me.rows(), Ms.rows());
  assert_eq(Me.rows(), Ms.cols());
  const MatrixXd mm = Ms-MatrixXd(Me);
  assert_lt(mm.norm(),1e-8);
  cout<< "(M-Mref).norm() = " << mm.norm() << endl;
}

void testKSMall(){
  
  const string tet_fname = std::string(TEST_DATA_DIR)+"one_tet.abq";
  const string tet_mtl_fname = std::string(TEST_DATA_DIR)+"one_tet.elastic";
  pTetMesh mesh = pTetMesh(new TetMesh);
  assert ( mesh->load(tet_fname) );
  assert ( mesh->loadElasticMtl(tet_mtl_fname) );
  const int n = mesh->nodes().size();

  VectorXd x0(n*3);
  mesh->nodes(x0);

  ElasticForceTetFullStVK elas_c(mesh);
  assert(elas_c.prepare());
  const MatrixXd ref_Km = elas_c.K(x0);

  ComputeStiffnessMat elas(mesh);
  assert(elas.prepare());
  elas.K(x0);
  const SXMatrix &K = elas.K(x0);
  MatrixXd Km;
  CASADI::convert(K,Km);
  
  assert_eq(Km.rows(),n*3);
  assert_eq(Km.cols(),n*3);
  assert_lt((Km-ref_Km).norm(),1e-5);
  
  cout<< "(K-K_ref).norm(): " << (Km-ref_Km).norm() << endl;
  cout << "test end.." << endl;
}


#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdlib.h>
#include <float.h>
#include "QPSolver.h"
using namespace std;

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

int main(int argc, char *argv[]){

  FixedSparseMatrix<double> A;
  MPRGPQPSolver<double>::Vec B(1),L(1),H(1);
  A.resize(1,1);
  A.insert(0,0) = 1;
  B << -1;
  L << 1.0f;
  H << 10;

  // solve for min 1/2*x^t*A*x-x^t*B, s.t.  L<=x<=H.
  MPRGPQPSolver<double> f(A,B,L,H);
  f.setCallback(boost::shared_ptr<Callback<double> >(new Callback<double>()));
  MPRGPQPSolver<double>::Vec x(1);
  x << 3.0;
  f.solve(x);
  printf("\n");
  cout << x.transpose() << endl;

  return 0;
}


// void recoverFineSmooth(){

//   TRACE_FUN();

//   const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";
//   MatrixXd eig_W;
//   VectorXd eig_lambda;
//   bool succ = true;
//   { // load opt value
// 	succ = load(data_root+"tempt/opt_W_smooth.b",eig_W); assert(succ);
// 	succ = load(data_root+"tempt/opt_lambda_smooth.b",eig_lambda); assert(succ);
//   }

//   { // load real value
// 	// succ = load(data_root+"tempt/eigenvectors_smooth.b",eig_W); assert(succ);
// 	// succ = load(data_root+"tempt/eigenvalues_smooth.b",eig_lambda); assert(succ);
//   }

//   MatrixXd W(eig_W.rows(),3);
//   W.col(0) = eig_W.col(0);
//   W.col(1) = eig_W.col(1);
//   W.col(2) = eig_W.col(6);
//   VectorXd lambda(3);
//   lambda << eig_lambda[0], eig_lambda[1], eig_lambda[6];
//   cout<< "\n\nLambda: " << lambda.transpose() << "\n\n";

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth2",
// 			W,lambda,1e-3,1e-3,1,1e12, true);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth2",
// 			W,lambda,1e-3,1e-3,1,1e12, false);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth3",
// 			W,lambda,1e-6,1e-6,1,1e12, true);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth3",
// 			W,lambda,1e-6,1e-6,1,1e12, false);
//   }

//   {
//   	MaterialFitting_MA_K mtlfit_k;
//   	mtlfit_k.useHessian(false);
//   	recover(&mtlfit_k, data_root+"model/",
//   			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth",
//   			W,lambda,100,100,1,1e12, true);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth",
// 			W,lambda,100,100,1,1e12, false);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth1",
// 			W,lambda,1,1,1,1e12, true);
//   }

//   {
// 	MaterialFitting_MA_K mtlfit_k;
// 	recover(&mtlfit_k, data_root+"model/",
// 			data_root+"model/tempt_mesh.elastic","./tempt/test_fine_smooth1",
// 			W,lambda,1,1,1,1e12, false);
//   }

// }

// void recoverFineNonSmooth(){

//   TRACE_FUN();

//   const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";

//   MatrixXd eig_W, W;
//   VectorXd eig_lambda, lambda;

//   bool succ = true;
//   { // load opt value
// 	succ = load(data_root+"tempt/opt_W.b",eig_W); assert(succ);
// 	succ = load(data_root+"tempt/opt_lambda.b",eig_lambda); assert(succ);
// 	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
//   }

//   lambda = eig_lambda;
//   W = eig_W;

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth32",
// 			W,lambda,1e-3,1e-3,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth33",
// 			W,lambda,1e-1,1e-1,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth34",
// 			W,lambda,100,100,100,1e12,true);
//   }

//   exit(0);

//   lambda.resize(2);
//   W.resize(eig_W.rows(),2);

//   lambda << eig_lambda[0], eig_lambda[1];
//   W.col(0) = eig_W.col(0);
//   W.col(1) = eig_W.col(1);

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth16",
// 			W,lambda,1e-5,1e-5,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth16",
// 			W,lambda,1e-5,1e-5,100,1e12,false);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth15",
// 			W,lambda,1e-7,1e-7,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth15",
// 			W,lambda,1e-7,1e-7,100,1e12,false);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth17",
// 			W,lambda,1.0,1.0,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth17",
// 			W,lambda,1.0,1.0,100,1e12,false);
//   }

//   lambda.resize(3);
//   W.resize(eig_W.rows(),3);

//   lambda << eig_lambda[0], eig_lambda[1], eig_lambda[6];
//   W.col(0) = eig_W.col(0);
//   W.col(1) = eig_W.col(1);
//   W.col(2) = eig_W.col(6);

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth18",
// 			W,lambda,1e-7,1e-7,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth18",
// 			W,lambda,1e-7,1e-7,100,1e12,false);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth19",
// 			W,lambda,1e-5,1e-5,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth19",
// 			W,lambda,1e-5,1e-5,100,1e12,false);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth20",
// 			W,lambda,1.0,1.0,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth20",
// 			W,lambda,1.0,1.0,100,1e12,false);
//   }

//   { // load real value
// 	succ = load(data_root+"tempt/eigenvectors.b",eig_W); assert(succ);
// 	succ = load(data_root+"tempt/eigenvalues.b",eig_lambda); assert(succ);
// 	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
//   }

//   lambda.resize(2);
//   W.resize(eig_W.rows(),2);

//   lambda << eig_lambda[0], eig_lambda[1];
//   W.col(0) = eig_W.col(0);
//   W.col(1) = eig_W.col(1);

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth152",
// 			W,lambda,1e-7,1e-7,100,1e12,true);
//   }

//   {
// 	MaterialFitting_EV_MA_K mtlfit_k;
// 	mtlfit_k.useHessian(false);
// 	recover(&mtlfit_k, data_root+"model/",data_root+"model/mesh.elastic",
// 			"./tempt/test_fine_non_smooth152",
// 			W,lambda,1e-7,1e-7,100,1e12,false);
//   }

// }
