#include <MatrixTools.h>
#include <MatrixIO.h>
#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <Timer.h>
#include "MaterialFitting.h"
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace ELASTIC_OPT;

void computeEigenValues(const string data_root, const string mtl_file,
						const int eigenNum,const set<int> &fixednodes,
						const string save_to){
  
  // load data.
  INFO_LOG("load data");
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(mtl_file); assert(succ);
  const int num_tet = (int)tetmesh->tets().size();
  const int n = tetmesh->nodes().size();

  // save material
  INFO_LOG("save material");
  succ = tetmesh->writeElasticMtlVTK("./tempt/material_correct"); assert(succ);

  // compute K
  INFO_LOG("compute K");
  ElasticForceTetFullStVK elas(tetmesh);
  assert(elas.prepare());
  VectorXd x0(n*3);
  tetmesh->nodes(x0);
  SparseMatrix<double> K = elas.K(x0);
  K *= -1.0f;

  // compute mass matrix
  INFO_LOG("compute mass matrix");
  MassMatrix mass;
  DiagonalMatrix<double,-1> M;
  mass.compute(M,*tetmesh);

  // remove fixed nodes
  SparseMatrix<double> P;
  EIGEN3EXT::genReshapeMatrix(K.rows(),3,fixednodes,P);
  K = P*(K*P.transpose());
  const SparseMatrix<double> Msub = P*(M*P.transpose());
  const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(Msub);

  /// compute W, lambda
  INFO_LOG("compute W, lambda");
  MatrixXd W_full;
  VectorXd lambda_full;
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W_full,lambda_full,eigenNum); 
  assert(succ);

  // const MatrixXd W = W_full.rightCols(lambda_full.size()-6);
  // const VectorXd lambda = lambda_full.tail(lambda_full.size()-6);

  // test Wt*M*W
  INFO_LOG("check results");
  const MatrixXd WtMWI = W_full.transpose()*Msub*W_full-MatrixXd::Identity(W_full.cols(),W_full.cols());
  const MatrixXd KW_MWLambda = K*W_full-(Msub*W_full)*lambda_full.asDiagonal();
  cout << "(WtMW-I).norm(): " << WtMWI.norm() << endl;
  cout << "(KW_MWLambda).norm(): " << KW_MWLambda.norm() << endl;
  cout << "eigenvalues: " << lambda_full.transpose() << endl;
  cout << "norm(Klower): " << Klower.norm() << endl;
  cout << "norm(M): " << Msub.norm() << endl;

  // add fixed nodes to W, and lambda, then save
  const MatrixXd W = P.transpose()*W_full;
  const VectorXd lambda = lambda_full;
  succ = write(save_to+"lambda.b", lambda); assert(succ);
  succ = write(save_to+"W.b", W); assert(succ);
  
  const MatrixXd K_dense = K;
  const MatrixXd M_dense = Msub;
  VectorXd M_diag(Msub.rows());
  for (int i = 0; i < M_dense.rows(); ++i){
	M_diag[i] = M_dense(i,i);
  }
  succ = write(save_to+"K.b", K_dense); assert(succ);
  succ = write(save_to+"M.b", M_diag); assert(succ);
  succ = tetmesh->write(save_to+"mesh.vol"); assert(succ);
}

void recover(MaterialFitting *mtlfit_k,
			 const string data_root, const string mtl_file, const string save_to,
			 const MatrixXd &W, const VectorXd &lambda, const MatrixXd &Scale,
			 const double Gs=1e-3, const double Gl=1e-3,
			 const double lower=1, const double upper=1e10, const bool use_ipopt=true){

  { // fit G, l
  	INFO_LOG("fit G, l");
  	mtlfit_k->loadTetMesh(data_root+"mesh.abq");
  	mtlfit_k->loadMtl(mtl_file);
  	mtlfit_k->loadFixednodes(data_root+"/con_nodes.bou");
  	mtlfit_k->setWLambda(W, lambda);
	mtlfit_k->setScaleMatrix(Scale);
  	mtlfit_k->computeK();
  	mtlfit_k->computeM();

	// save input W, lambda, M , fixed nodes, and mesh.vol
	bool succ = mtlfit_k->saveAllInputs(save_to);  assert(succ);
	// save K
	SparseMatrix<double> K0;
  	mtlfit_k->computeK(K0);
	const MatrixXd dense_K = K0;

  	mtlfit_k->removeFixedDOFs();
  	mtlfit_k->setBounds(lower,upper);
  	mtlfit_k->setMuSmoothGL(Gs, Gl);
  	mtlfit_k->setMuSmoothEv(0, 0);
  	mtlfit_k->setMuSmoothDensity(0.0f);
  	mtlfit_k->setMuAverageDensity(0.0f);
  	mtlfit_k->assembleObjfun();

	if (use_ipopt){
	  mtlfit_k->solveByIpopt();
	  mtlfit_k->saveResults(save_to);
	}else{
	  mtlfit_k->solveByMPRGP(save_to+"_x.b");
	  mtlfit_k->saveResults(save_to+"_MPRGP_");
	}
  	mtlfit_k->printResult();
  }

  { // check resutls
	SparseMatrix<double> K, K0;
  	mtlfit_k->computeK(K,mtlfit_k->getShearGResult(),mtlfit_k->getLameResult());
  	mtlfit_k->computeK(K0);

	DiagonalMatrix<double,-1> M0_diag, M_diag;
  	mtlfit_k->computeM(M_diag, mtlfit_k->getDensityResult());
  	mtlfit_k->computeM(M0_diag);

	const SparseMatrix<double> P = mtlfit_k->getMatrixForRemovingFixedDOFs();
	K = P*K*P.transpose();
	K0 = P*K0*P.transpose();
	const MatrixXd W_sub = P*W;
	const SparseMatrix<double> M0 = P*M0_diag*P.transpose();
	const SparseMatrix<double> M = P*M_diag*P.transpose();

  	const MatrixXd la = lambda.asDiagonal();
	const MatrixXd m2 = (W_sub.transpose()*K0*W_sub-(W_sub.transpose()*M0*W_sub)*la);
	const MatrixXd m3 = (W_sub.transpose()*K*W_sub-(W_sub.transpose()*M0*W_sub)*la);
	const MatrixXd m4 = Scale.cwiseProduct(m2);
	const MatrixXd m5 = Scale.cwiseProduct(m3);

  	cout << "\n\nW^t*K0*W-W^t*M0*W*Lambda:  "<<m2.norm()<<"\n\n"<<m2<<"\n\n";
  	cout << "\n\nW^t*K*W-W^t*M0*W*Lambda:  "<<m3.norm()<<"\n\n"<<m3<<"\n\n";
  	cout << "\n\nS:W^t*K0*W-W^t*M0*W*Lambda:  "<<m4.norm()<<"\n\n"<<m4<<"\n\n";
  	cout << "\n\nS:W^t*K*W-W^t*M0*W*Lambda:  "<<m5.norm()<<"\n\n"<<m5<<"\n\n";

	MatrixXd W;
	VectorXd lambda;
	const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
	const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(M);
	bool succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W,lambda,20);
	assert(succ);
	cout << cout.precision(6) << "\neigenvalues: " << lambda.transpose() << "\n\n";
  }

}

void loadAndScale(const string data_root_opt, MatrixXd &W, VectorXd &lambda,MatrixXd &C,const double h, const double rho, const int r = -1){
  
  MatrixXd S, eig_W;
  VectorXd eig_lambda;
  // load W and lambda
  bool succ = load(data_root_opt+"opt_scaled_W.b",eig_W); assert(succ);
  succ = load(data_root_opt+"opt_S.b",S); assert(succ);
  const MatrixXd t = eig_W;
  eig_W = t.leftCols(S.rows())*S;
  eig_W *= 1.0f/rho;
  
  succ = load(data_root_opt+"opt_lambda.b",eig_lambda); assert(succ);
  eig_lambda *= rho/(h*h);

  // use the first r modes.
  if ( r > 0 && r < eig_W.cols()){
	W = eig_W;
	eig_W = W.leftCols(r);
	lambda = eig_lambda;
	eig_lambda = lambda.head(r);
  }

  // compute C
  MatrixXd Z,F;
  succ = load(data_root_opt+"opt_Z.b",Z); assert(succ);
  succ = load(data_root_opt+"opt_F.b",F); assert(succ);
  if (Z.rows() > eig_lambda.size()){
	const MatrixXd tz = Z;
	const MatrixXd tf = F;
	Z = tz.topRows(eig_lambda.size());
	F = F.topRows(eig_lambda.size());
  }

  const int T = 150;
  VectorXd scale(eig_lambda.size());
  double s_max = -1.0f;
  for (int i = 0; i < Z.rows(); ++i){
	const double dz = Z.row(i).head(T).norm();
	const double df = F.row(i).head(T).norm();
	double s = dz/(df*eig_lambda[i]);
	scale[i] = s*s;
	if (scale[i] > s_max)
	  s_max = scale[i];
	cout<<i<<",\t" << dz<<",\t" << df << ",\t" << scale[i]<<",\t"<<eig_lambda[i]<<"\n";
  }
  scale *= 1.0f/s_max;

  C.resize(scale.size(), scale.size());
  for (int i = 0; i < C.rows(); ++i){
    for (int j = i; j < C.cols(); ++j)
	  C(j,i) = C(i,j) = scale[i]>scale[j] ? scale[j]:scale[i];
  }

  W = eig_W;
  lambda = eig_lambda;
  cout<< "\n\nLambda: " << lambda.transpose() << "\n\n";
  cout << "Scale matrix: \n\n" << C << "\n\n";
}

void recoverRealCoarseSmooth(){
  
  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";

  MatrixXd W, eig_W, Scale;
  VectorXd lambda, eig_lambda, mu_e;

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/mesh_smooth_small.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/coarse/coarse_real_smooth"); assert(succ);
  }

  { // load real value
  	bool succ = load(data_root+"/eigenvectors_smooth_mtl_80.b",eig_W); assert(succ);
  	succ = load(data_root+"/eigenvalues_smooth_mtl_80.b",eig_lambda); assert(succ);
  	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
  }

  const int r = 2;
  Scale = MatrixXd::Ones(r,r);

  W = eig_W.leftCols(r);
  lambda = eig_lambda.head(r);
  mu_e.resize(6);
  mu_e << 0, 1e-8, 1e-2, 1.0, 10, 100;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	mtlfit_k.setFullSpaceConPen(mu_e[i]);
	recover(&mtlfit_k,data_root, data_root_opt+"/mesh.elastic", 
			"./tempt/coarse/s/test_coarse_smooth_null_p"+TOSTR(mu_e[i]),
			W, lambda, Scale, 1e-4, 1e-4 , 1e-2, 1e16, true);
  }

}

void recoverRealCoarseNonSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";

  computeEigenValues(data_root, data_root+"/mesh.elastic", 16, 
					 set<int>(), "./tempt/data/no_con_coarse_non_smooth_ma_");

  MatrixXd W, eig_W, Scale;
  VectorXd lambda, eig_lambda, mu_e;

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/mesh.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/coarse/coarse_real"); assert(succ);
  }

  { // load real value
  	bool succ = load(data_root+"/eigenvectors_2mtl_80.b",eig_W); assert(succ);
  	succ = load(data_root+"/eigenvalues_2mtl_80.b",eig_lambda); assert(succ);
  	cout << "lambda: " << eig_lambda.head(20).transpose() << "\n\n";
  }

  const int r = 10;
  Scale = MatrixXd::Ones(r,r);

  W = eig_W.leftCols(r);
  lambda = eig_lambda.head(r);

  mu_e.resize(2);
  mu_e << 1e-4, 1e-2;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	mtlfit_k.setFullSpaceConPen(mu_e[i]);
	recover(&mtlfit_k,data_root, data_root+"/mesh.elastic", 
			// "./tempt/coarse/ns/test_coarse_non_smooth_ev_null_p"+TOSTR(mu_e[i]),
			"./tempt/data/coarse_non_smooth_ma_",
			W, lambda, Scale, 0, 0, 1e-3, 1e12, true);
  }
  
}

void recoverCoarseNonSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam-coarse/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam-coarse/model/";

  MatrixXd W, Scale;
  VectorXd lambda, mu_e;
  loadAndScale(data_root_opt, W, lambda, Scale, 0.03f,1.0f, 10);

  mu_e.resize(3);
  mu_e << 1e-8, 1e-4, 1e-2;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	recover(&mtlfit_k,data_root, data_root_opt+"/mesh.elastic", 
			"./tempt/data/coarse_non_smooth_opt_",
			W, lambda,Scale,mu_e[i],mu_e[i],1e-2,1e10,true);
  }

}

void recoverRealFineSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";

  MatrixXd eig_W, W, Scale;
  VectorXd eig_lambda, lambda, mu_e;

  computeEigenValues(data_root+"model/", data_root+"model/mtl_smooth.elastic", 16,
					 set<int>(), "./tempt/data/no_con_fine_smooth_ma_");
  { // load real value
  	bool succ = load(data_root+"tempt/eigenvectors_smooth.b",eig_W); assert(succ);
  	succ = load(data_root+"tempt/eigenvalues_smooth.b",eig_lambda); assert(succ);
  	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
  }

  const int r = 10;
  Scale = MatrixXd::Ones(r,r);
  W = eig_W.leftCols(r);
  lambda = eig_lambda.head(r);

  mu_e.resize(3);
  mu_e << 0, 1e-6, 1e-10;

  for (int i = 0; i < mu_e.size(); ++i){
	
	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	recover(&mtlfit_k,data_root+"model/", data_root+"model/mtl_smooth.elastic", 
			"./tempt/fine_smooth/test_fine_smooth_real_full"+TOSTR(mu_e[i]), 
			W, lambda, Scale, mu_e[i],mu_e[i],1,1e12,true);
  }
}

void recoverFineSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";

  MatrixXd eig_W, W, Scale;
  VectorXd eig_lambda, lambda;
  loadAndScale(data_root_opt+"smooth/", eig_W, eig_lambda, Scale, 0.04f, 100.0f);

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/mtl_smooth.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/fine_smooth/fine_smooth_real"); assert(succ);
  }

  W = eig_W;
  lambda = eig_lambda;
  VectorXd mu_e(1);
  mu_e << 1e-6;
  for (int i = 0; i < mu_e.size(); ++i){
	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	mtlfit_k.setFullSpaceConPen(0);
	recover(&mtlfit_k,data_root, data_root+"/tempt_mesh.elastic", 
			"./tempt/fine_smooth/check_test_fine_smooth_ev"+TOSTR(mu_e[i]),
			// "./tempt/data/fine_smooth_opt_",
			W, lambda, Scale, mu_e[i],mu_e[i],1.0,1e8,true);
  }

}

void recoverRealFineNonSmooth(){
  
  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";

  MatrixXd W, eig_W, Scale;
  VectorXd lambda, eig_lambda, mu_e;

  computeEigenValues(data_root+"model/", data_root+"model/mesh.elastic", 16, 
					 set<int>(), "./tempt/data/no_con_fine_non_smooth_ma_");
  { // load real value
  	bool succ = load(data_root+"tempt/eigenvectors.b",eig_W); assert(succ);
  	succ = load(data_root+"tempt/eigenvalues.b",eig_lambda); assert(succ);
  	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
  }

  const int r = 10;
  Scale = MatrixXd::Ones(r,r);
  W = eig_W.leftCols(r);
  lambda = eig_lambda.head(r);

  mu_e.resize(3);
  mu_e << 0, 1e-6, 1e-10;
  for (int i = 0; i < mu_e.size(); ++i){
	
	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	recover(&mtlfit_k,data_root+"/model/", data_root+"/model/mesh.elastic", 
			"./tempt/fine/test_fine_non_smooth_real_full"+TOSTR(mu_e[i]), 
			W, lambda, Scale, mu_e[i],mu_e[i],1.0,1e12,true);
  }
}

void recoverFineNonSmooth(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/beam/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/beam_fine/model/";

  MatrixXd W, Scale;
  VectorXd lambda;
  loadAndScale(data_root_opt, W, lambda, Scale, 0.04f,100.0f);

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"/model/mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/model/mesh.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/fine/fine_real"); assert(succ);
  }

  VectorXd mu_e(4);
  mu_e << 1, 50, 1e-6, 1e-4;
  for (int i = 0; i < mu_e.size(); ++i){
	
	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(false);
	recover(&mtlfit_k,data_root+"/model/", data_root+"/model/tempt_mesh.elastic", 
			"./tempt/fine/check_fine_non_smooth_ev"+TOSTR(mu_e[i]),
			// "./tempt/data/fine_non_smooth_opt_",
			W, lambda, Scale, mu_e[i],mu_e[i], 1.0, 1e8,true);
  }

}

void recoverDino(){

  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/dino-mtlopt/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/dino_scaled/model/";

  MatrixXd W, Scale;
  VectorXd lambda, mu_e;
  loadAndScale(data_root_opt, W, lambda, Scale, 0.03f,1.0f, -1);

  mu_e.resize(1);
  mu_e << 1e-6;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	recover(&mtlfit_k,data_root, data_root_opt+"/mesh.elastic", 
			"./tempt/data/coarse_non_smooth_opt_",
			W, lambda,Scale,mu_e[i],mu_e[i],1e-2,1e10,true);
  }
}

void recoverRealDino(){
  
  TRACE_FUN();

  const string data_root = "/home/simba/Workspace/SolidSimulator/data/dino-mtlopt/model/";
  const string data_root_opt ="/home/simba/Workspace/AnimationEditor/Data/dino_scaled/model/";

  MatrixXd W, eig_W, Scale;
  VectorXd lambda, eig_lambda, mu_e;

  { //save real mateial
	TetMesh tetmesh;
	bool succ = tetmesh.load(data_root+"mesh.abq"); assert(succ);
	succ = tetmesh.loadElasticMtl(data_root+"/elastic.elastic"); assert(succ);
	succ = tetmesh.writeElasticMtlVTK("./tempt/dino/"); assert(succ);
  }

  { // load real value
  	bool succ = load(data_root+"/eigenvectors.b",eig_W); assert(succ);
  	succ = load(data_root+"/eigenvalues.b",eig_lambda); assert(succ);
  	cout << "lambda: " << eig_lambda.transpose() << "\n\n";
  }

  const int r = 2;
  Scale = MatrixXd::Ones(r,r);

  W = eig_W.leftCols(r);
  lambda = eig_lambda.head(r);
  mu_e.resize(1);
  mu_e << 1e-6;
  for (int i = 0; i < mu_e.size(); ++i){

	MaterialFitting_EV_MA_K mtlfit_k;
	mtlfit_k.useHessian(true);
	mtlfit_k.setFullSpaceConPen(mu_e[i]);
	recover(&mtlfit_k,data_root, data_root_opt+"/mesh.elastic", 
			"./tempt/dino_ma"+TOSTR(mu_e[i]),
			W, lambda, Scale, mu_e[i], mu_e[i] , 1e-2, 1e16, false);
  }
}

int main(int argc, char *argv[]){

  // recoverRealCoarseSmooth();
  // recoverRealCoarseNonSmooth();
  // recoverCoarseNonSmooth();
  // recoverCoarseSmooth();
  // recoverFineSmooth();
  // recoverFineNonSmooth();
  // recoverRealFineNonSmooth();
  // recoverRealFineSmooth();

  // recoverDino();
  // recoverRealDino();


  // load mesh
  const string data_root = "/home/simba/Workspace/SolidSimulator/data/dino-mtlopt/model/";
  pTetMesh tetmesh = pTetMesh(new TetMesh());
  bool succ = tetmesh->load(data_root+"mesh.abq"); assert(succ);
  succ = tetmesh->loadElasticMtl(data_root+"/mesh.elastic"); assert(succ);
  INFO_LOG("save material");
  succ = tetmesh->writeElasticMtlVTK("./tempt/material_correct"); assert(succ);

  // compute M
  INFO_LOG("compute mass matrix");
  MassMatrix mass;
  SparseMatrix<double> M;
  mass.compute(M,*tetmesh);

  // compute K
  INFO_LOG("compute stiffness matrix");
  INFO_LOG("compute K");
  ElasticForceTetFullStVK elas(tetmesh);
  assert(elas.prepare());
  VectorXd x0(tetmesh->nodes().size()*3);
  tetmesh->nodes(x0);
  SparseMatrix<double> K = elas.K(x0);
  K *= -1.0f;

  // solve eigen value problem
  const SparseMatrix<double> Klower = EIGEN3EXT::getLower(K);
  const SparseMatrix<double> Mlower = EIGEN3EXT::getLower(M);
  const int eigenNum = 20;
  MatrixXd W;
  VectorXd lambda;
  succ = EigenSparseGenEigenSolver::solve(Klower,Mlower,W,lambda,eigenNum); assert(succ);

  // save data
  ofstream out;
  out.open("dino.stiffness");
  out << K.rows()<< "\t" << K.cols();
  for (int k=0; k<K.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(K,k); it; ++it){
	  out << endl;
	  out <<it.row()<< " "<<it.col()<<" " << it.value();
	}
  }
  out.close();
    
  out.open("dino.mass");
  out << M.rows()<< "\t" << M.cols();
  for (int k=0; k<M.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(M,k); it; ++it){
	  out << endl;
	  out <<it.row()<< " "<<it.col()<< " " << it.value();
	}
  }
  out.close();

  out.open("dino.eigenvalues");
  out << lambda.size();
  for (int i = 0; i < lambda.size(); ++i){
    out << endl << lambda[i];
  }
  out.close();

  out.open("dino.eigenvectors");
  out << W.rows()<< "\t" << W.cols();
  for (int r = 0; r < W.rows(); ++r){
	out << endl;
	for (int c = 0; c < W.cols(); ++c){
	  out << W(r,c);
	  if (c != W.cols()-1)
		out << " ";
	}
  }
  out.close();

  return 1;

}
