#include <SparseGenEigenSolver.h>
#include <SparseMatrixTools.h>
#include <MassMatrix.h>
#include <ElasticForceTetFullStVK.h>
#include "MaterialFitting.h"
using namespace UTILITY;
using namespace ELASTIC_OPT;

void MaterialFitting::computeK(){
  
  TRACE_FUN();
  ComputeStiffnessMat elas(tetmesh);
  initShearG(tetmesh->tets().size(),G);
  initLame(tetmesh->tets().size(),Lame);
  elas.setMaterial(G, Lame);
  bool succ = elas.prepare();  assert(succ);
  VectorXd x0;
  tetmesh->nodes(x0);

  Timer timer;
  timer.start();
  K = elas.K(x0);
  K *= -1.0f;
  timer.stop("computeK: ");
}

void MaterialFitting::computeM(){
  
  TRACE_FUN();
  INFO_LOG("compute mass matrix");
  ComputeMassMat mass;
  initDensity(tetmesh->tets().size(),rho);
  Timer timer;
  timer.start();
  mass.computeDiag(M,*tetmesh,rho);
  // mass.compute(M,*tetmesh,rho);
  timer.stop("computeM: ");
}

void MaterialFitting::computeK(SparseMatrix<double> &K)const{
  
  ElasticForceTetFullStVK elas_c(tetmesh);
  const bool succ = elas_c.prepare(); assert(succ);
  VectorXd x0(tetmesh->nodes().size()*3);
  tetmesh->nodes(x0);
  K = elas_c.K(x0);
  K = -K;
}

void MaterialFitting::computeK(SparseMatrix<double> &K,const vector<double> &G,const vector<double> &L)const{

  const vector<double> G0 = tetmesh->material()._G;
  const vector<double> L0 = tetmesh->material()._lambda;
  
  tetmesh->material()._G = G;
  tetmesh->material()._lambda = L;

  ElasticForceTetFullStVK elas_c(tetmesh);
  const bool succ = elas_c.prepare(); assert(succ);
  VectorXd x0(tetmesh->nodes().size()*3);
  tetmesh->nodes(x0);
  K = elas_c.K(x0);
  K = -K;

  tetmesh->material()._G = G0;
  tetmesh->material()._lambda = L0;
}

void MaterialFitting::computeM(DiagonalMatrix<double,-1>&M)const{

  assert(tetmesh);
  MassMatrix mass;
  mass.compute(M,*tetmesh);
}

void MaterialFitting::computeM(DiagonalMatrix<double,-1>&M,const vector<double>&rh)const{
  
  assert(tetmesh);
  MassMatrix mass;
  const vector<double> r0 = tetmesh->material()._rho;
  tetmesh->material()._rho = rh;
  mass.compute(M,*tetmesh);
  tetmesh->material()._rho = r0;
}

void MaterialFitting::computeM(SparseMatrix<double> &M)const{
  assert(tetmesh);
  MassMatrix mass;
  mass.compute(M,*tetmesh);
}

void MaterialFitting::computeM(SparseMatrix<double> &M,const vector<double> &rh)const{

  assert(tetmesh);
  MassMatrix mass;
  const vector<double> r0 = tetmesh->material()._rho;
  tetmesh->material()._rho = rh;
  mass.compute(M,*tetmesh);
  tetmesh->material()._rho = r0;
}

SparseMatrix<double> MaterialFitting::getMatrixForRemovingFixedDOFs()const{

  assert(tetmesh);
  const int nx3 = tetmesh->nodes().size()*3;
  SparseMatrix<double> Pm;
  if (fixednodes.size() > 0){
	EIGEN3EXT::genReshapeMatrix(nx3,3,fixednodes,Pm);
  }else{
	Pm = EIGEN3EXT::eye(nx3, 1.0);
  }
  return Pm;
}

void MaterialFitting::removeFixedDOFs(){
  
  if (fixednodes.size() > 0){

	SparseMatrix<double> Pm;
	EIGEN3EXT::genReshapeMatrix(K.size1(),3,fixednodes,Pm);
	SXMatrix P;
	CASADI::convert(Pm,P);
	assert_eq(P.size1(),K.size1()-fixednodes.size()*3);
	assert_eq(P.size2(),K.size1());
	K = P.mul(K.mul(trans(P)));
	M = P.mul(M.mul(trans(P)));
	W = P.mul(W);
  }
}

void MaterialFitting::getInitValue(VectorXd &init_x)const{
  
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> &rh = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();

  init_x.resize(num_tet*3);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
    init_x[i+num_tet*2] = rh[i];
  }
}

vector<double> MaterialFitting::getShearGResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet);
  return vector<double>(rlst.begin(),rlst.begin()+num_tet);
}

vector<double> MaterialFitting::getLameResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet*2);
  return vector<double>(rlst.begin()+num_tet,rlst.begin()+2*num_tet);
}

vector<double> MaterialFitting::getDensityResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet*3);
  return vector<double>(rlst.begin()+2*num_tet,rlst.begin()+3*num_tet);
}

void MaterialFitting::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);
  addFullSpaceConstraints(objfun);

  const SXMatrix M1 = assembleObjMatrix();
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j)
	  objfun += M1.elem(i,j)*M1.elem(i,j)*ScaleMat(i,j)*ScaleMat(i,j);
  }

  objfun *= scalor;
}

SXMatrix MaterialFitting::assembleObjMatrix(){

  const SXMatrix WtMWL = trans(W).mul(M.mul(W)).mul(lambda);
  const SXMatrix M1 = (trans(W).mul(K.mul(W)) - WtMWL);
  return M1;
}

void MaterialFitting::addFullSpaceConstraints(SX &objfun)const{

  TRACE_FUN();

  if (fullspace_con_pen <= 0.0f){
	return ;
  }

  Timer timer;
  timer.start("start compute kernel");
  const MatrixXd W = CASADI::convert<double>(this->W).transpose();
  const MatrixXd W_kernel = W.fullPivLu().kernel();
  timer.stop("finish compute kernel: ");
  assert_gt(W_kernel.norm(), 0);
  
  timer.start("start compute full con");
  INFO_LOG("wait until i = " << W_kernel.cols()<<"\n\ni = ");
  for (int i = 0; i < W_kernel.cols(); ++i){
	const VectorXd wi = W_kernel.col(i);
	const SXMatrix wis = CASADI::convert(wi);
	objfun += fullspace_con_pen/trans(wis).mul(K.mul(wis)).elem(0,0);
	if (i % 50 == 0){
	  cout << i << ", ";
	  cout.flush();
	}
  }
  timer.stop("\n\nfinish compute full con: ");
  
  // VSX x;
  // initAllVariables(x);
  // for (int i = 0; i < x.size(); ++i){
  //   objfun += fullspace_con_pen/x[i];
  // }
}

void MaterialFitting::computeSmoothObjFunctions(vector<SX> &funs)const{

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  funs.clear();
  funs.reserve(neigh_tet.size()*4);

  const VSX G = getG();
  const VSX Lame = getLame();
  const VSX rho = getRho();
  const VSX E = getE();
  const VSX v = getV();
  for (int i = 0; i < neigh_tet.size(); ++i){

	const double vi = tetmesh->volume(i);
	assert_gt(vi,0);
	for (int k = 0; k < 4; ++k){
	  const int j = neigh_tet[i][k];
	  if (j >= 0){
		const double vj = tetmesh->volume(j);
		assert_gt(vj,0);
		assert_ne(i,j);
		if (mu_neigh_G > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_G*(vi+vj))*(G[i]-G[j]) );
		if (mu_neigh_L > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_L*(vi+vj))*(Lame[i]-Lame[j]) );
		if (mu_neigh_rho > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_rho*(vi+vj))*(rho[i]-rho[j]) );
		if (mu_neigh_E > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_E*(vi+vj))*(E[i]-E[j]) );
		if (mu_neigh_poission > 0)
		  funs.push_back( 0.5f*sqrt(mu_neigh_poission*(vi+vj))*(v[i]-v[j]) );
	  }
	}
  }
}

void MaterialFitting::addSmoothObjfun(SX &objfun)const{
  
  vector<SX> funs;
  computeSmoothObjFunctions(funs);
  for (int i = 0; i < funs.size(); ++i){
    objfun += funs[i]*funs[i];
  }
}

void MaterialFitting::addAverageDensityObjfun(SX &objfun)const{
  
  const vector<double> &desired_rho = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();
  assert_eq(desired_rho.size(), num_tet);
  assert_eq(rho.size(), num_tet);
  double desired_mass = 0.0f;
  SX mass = 0.0f;
  for (int i = 0; i < num_tet; ++i){
    const double v = tetmesh->volume(i);
	assert_gt(v,0.0f);
	desired_mass += desired_rho[i]*v;
	mass += rho[i]*v;
  }
  assert_gt(desired_mass, 0.0f);
  
  objfun += mu_average_rho*(desired_mass-mass)*(desired_mass-mass);
}

void MaterialFitting::addFixedNodesObjfun(SX &objfun)const{
  
  if(fixednodes.size() <= 0){
	return ;
  }

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  const VVec4i &tets = tetmesh->tets();
  assert_eq(neigh_tet.size(), tets.size());
  assert_eq(G.size(), tets.size());
  assert_eq(Lame.size(), tets.size());
  assert_eq(rho.size(), tets.size());
  for (int i = 0; i < tets.size(); ++i){
	int findnodes = 0;
    for (int j = 0; j < 4; ++j){
	  if ( fixednodes.find(tets[i][j])!=fixednodes.end() ){
		findnodes ++;
	  }
	}
	if (4 == findnodes){
	  const double vi = tetmesh->volume(i);
	  for (int k = 0; k < 4; ++k){
		const int j = neigh_tet[i][k];
		if (j >= 0){
		  const double vj = tetmesh->volume(j);
		  assert_gt(vj,0);
		  assert(i!=j);
		  const double scalor = 1e-8; /// @todo how to choose scalor?
		  objfun += 0.25f*(vi+vj)*(G[i]-G[j])*(G[i]-G[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(Lame[i]-Lame[j])*(Lame[i]-Lame[j])*scalor;
		  objfun += 0.25f*(vi+vj)*(rho[i]-rho[j])*(rho[i]-rho[j])*scalor;
		}
	  }
	}
  }
}

void MaterialFitting::scale(){
  
  // scale W
  const double w_norm = CASADI::convert<double>(W).norm();
  assert_gt(w_norm,0);
  W *= 1.0f/w_norm;
  
  // scale Lambda
  assert_gt(lambda.size(),0);
  scaled_lambda = lambda.elem(0,0).getValue();
  for (int i = 0; i < lambda.size1(); ++i){
  	const double la = lambda.elem(i,i).getValue();
  	assert_gt_ext(la,0.0f,i);
    if (la < scaled_lambda)
  	  scaled_lambda = la;
  }
  lambda *= 1.0f/scaled_lambda;

  // scale mass matrix
  scaled_mass = CASADI::convert<double>(M).norm();
  assert_gt(scaled_mass,0.0f);
  M *= 1.0f/scaled_mass;

}

void MaterialFitting::unScale(){

  assert_gt(scaled_mass,0);
  assert_gt(scaled_lambda,0);
  DEBUG_LOG("scaled_mass: "<< scaled_mass);
  DEBUG_LOG("scaled_lambda: "<< scaled_lambda);
  for (int i = 0; i < rlst.size(); ++i){
    rlst[i] = rlst[i]*scaled_mass*scaled_lambda;
  }
}

void MaterialFitting::hessGrad(MatrixXd &H, VectorXd &g)const{
  
  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  fun.init();

  CasADi::SXFunction grad_fun = CasADi::SXFunction(variable_x,fun.jac());
  grad_fun.init();

  CasADi::SXFunction hess_fun = CasADi::SXFunction(variable_x,fun.hess());
  hess_fun.init();

  VectorXd x0(variable_x.size());
  x0.setZero();
  CASADI::evaluate(grad_fun, x0, g);
  CASADI::evaluate(hess_fun, x0, H);

  assert_eq(g.norm(), g.norm());
  assert_eq(H.norm(), H.norm());
}

void MaterialFitting::hess(MatrixXd &H)const{

  // unimpelemented function
  ERROR_LOG("unimpelemented function");
  exit(0);

  FUNC_TIMER();
  
  // init variables
  VSX variable_x;
  initAllVariables(variable_x);

  // compute H for smooth functions
  SparseMatrix<double> smooth_H;
  {
	vector<SX> funs;
	computeSmoothObjFunctions(funs);
	SX smooth_obj = 0.0;
	for (int i = 0; i < funs.size(); ++i){
	  smooth_obj += funs[i]*funs[i];
	}

	CasADi::SXFunction smooth_fun = CasADi::SXFunction(variable_x,smooth_obj);
	smooth_fun.init();
	const SXMatrix sH = smooth_fun.hess();
	CASADI::convert(sH, smooth_H);
  }

  // compute H for ||W^t*K(x)*W - W^t*M0*W*\Lambda||
  {
	const int n = variable_x.size();
	H.resize(n,n);
	H.setZero();

	const VVec4i &tets = tetmesh->tets();
	const VectorUseti &nodeNeighNode = tetmesh->nodeNeighNode();
	const int num_tet = tets.size();
	assert_ge(n, num_tet);
	assert_eq(n % num_tet,0);

	MatrixXd W;
	CASADI::convert(this->W, W);
	const int r = W.cols();
	MatrixXd Mi(r,r), Mj(r,r);

	for (int i = 0; i < n; ++i){

	  Mi.setZero();
	  const int ti = i%num_tet;
	  for (int t = 0; t < 4; ++t){

		if (isFixed(tets[ti][t]))
		  continue;
		const int node_i = tets[ti][t]-fixedNodeBefore(tets[ti][t]);
		assert_in(node_i*3,0,K.size1()-3);
		BOOST_FOREACH(int node_neigh_i, nodeNeighNode[tets[ti][t]]){

		  assert_in(node_neigh_i,0,tetmesh->nodes().size()-1);
		  if (isFixed(node_neigh_i))
			continue;
		  node_neigh_i -= fixedNodeBefore(node_neigh_i);
		  assert_in(node_neigh_i*3,0,K.size2()-3);

		  SXMatrix Ki(3,3);
		  for (int tmp = 0; tmp < 3; ++tmp){
			Ki(0, tmp) = K(node_i*3+0, node_neigh_i*3+tmp);
			Ki(1, tmp) = K(node_i*3+1, node_neigh_i*3+tmp);
			Ki(2, tmp) = K(node_i*3+2, node_neigh_i*3+tmp);
		  }

		  Matrix3d Jm;
		  { // compute Jm
		  	CasADi::SXFunction fun = CasADi::SXFunction(variable_x[i],Ki);
		  	fun.init();
		  	const SXMatrix J = fun.jac();
		  	Jm(0,0) = J.elem(0,0).getValue(); 
		  	Jm(1,0) = J.elem(1,0).getValue(); 
		  	Jm(2,0) = J.elem(2,0).getValue();

		  	Jm(0,1) = J.elem(3,0).getValue(); 
		  	Jm(1,1) = J.elem(4,0).getValue(); 
		  	Jm(2,1) = J.elem(5,0).getValue();

		  	Jm(0,2) = J.elem(6,0).getValue(); 
		  	Jm(1,2) = J.elem(7,0).getValue(); 
		  	Jm(2,2) = J.elem(8,0).getValue();
		  }
		  Mi += W.block(node_neigh_i*3,0,3,r).transpose()*Jm*W.block(node_neigh_i*3,0,3,r);
		}
	  }

	  // for (int j = i; j < n; ++j){
	  // 	const int tj = j%num_tet;
	  // 	H(i,j) = H(j,i) = Map<Matrix<double, 36,1> >(&dKdXm(0,0)).dot(Map<Matrix<double, 36,1> >(&dKdXn(0,0)));
	  // }
	}	

  }

  assert_eq(H.rows(), smooth_H.rows());
  assert_eq(H.cols(), smooth_H.cols());
  H += smooth_H;
}

void MaterialFitting::grad(VectorXd &g)const{

  VSX variable_x;
  initAllVariables(variable_x);
  CasADi::SXFunction fun = CasADi::SXFunction(variable_x,objfun);
  fun.init();

  CasADi::SXFunction grad_fun = CasADi::SXFunction(variable_x,fun.jac());
  grad_fun.init();

  VectorXd x0(variable_x.size());
  x0.setZero();
  CASADI::evaluate(grad_fun, x0, g);
}

void MaterialFitting::print_NNLS_exit_code(const int exit_code)const{
  
  switch (exit_code){
  case 1: INFO_LOG("nnls solve success!");
	break;
  case 2: ERROR_LOG("nnls failed, the dimensions of the problem are bad, either m .le. 0 or n .le. 0.");
	break;
  case 3: ERROR_LOG("nnls failed, iteration count exceeded.");
	break;
  default: ERROR_LOG("nnls failed, unknown reason.");
	break;
  }
}

void MaterialFitting::testFixedNodes(){
  
  if(fixednodes.size() <= 0){
	return ;
  }

  const VVec4i &neigh_tet = tetmesh->faceNeighTet();
  const VVec4i &tets = tetmesh->tets();
  for (int i = 0; i < tets.size(); ++i){
	int findnodes = 0;
    for (int j = 0; j < 4; ++j){
	  if ( fixednodes.find(tets[i][j])!=fixednodes.end() )
		findnodes ++;
	}
	if (4 == findnodes){
	  tetmesh->material()._rho[i] *= 100.0f;
	  tetmesh->material()._G[i] *= 100.0f;
	  tetmesh->material()._lambda[i] *= 100.0f;
	}
  }
}

void MaterialFitting::saveResults(const string filename)const{
  
  const ElasticMaterial<double> mtl_0 = tetmesh->material();
  const vector<double> g = getShearGResult();
  const vector<double> l = getLameResult();
  const vector<double> r = getDensityResult();

  tetmesh->material()._G = g;
  tetmesh->material()._lambda = l;
  tetmesh->material()._rho = r;

  INFO_LOG("save results to: " << filename);

  bool succ = tetmesh->writeElasticMtlVTK(filename); assert(succ);
  succ = tetmesh->writeElasticMtl(filename+".elastic"); assert(succ);
  succ = writeVec(filename+"_x.b", rlst); assert(succ);

  tetmesh->material() = mtl_0;
}

void MaterialFitting::printResult()const{
  
  cout << "\nnumber of variables: " << rlst.size() << endl;
  for (int i = 0; i < rlst.size(); ++i)
    cout << rlst[i]<< " ";
  cout << "\n\n";
}

void MaterialFitting_MA_K::initDensity(const int num_tet, VSX &rho){

  const vector<double> &r = tetmesh->material()._rho;
  assert_eq(num_tet, r.size());
  rho.resize(num_tet);
  for (int i = 0; i < num_tet; ++i)	rho[i] = r[i];
}

void MaterialFitting_MA_K::getInitValue(VectorXd &init_x)const{

  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet*2);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = g[i];
    init_x[i+num_tet] = la[i];
  }
}

void MaterialFitting_Diag_KM::assembleObjfun(){
  
  TRACE_FUN();
  objfun = 0.0f;
  addSmoothObjfun(objfun);
  addAverageDensityObjfun(objfun);
  addFixedNodesObjfun(objfun);

  const MatrixXd I = MatrixXd::Identity(lambda.size(),lambda.size());
  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  const SXMatrix _Diag_KM = trans(W).mul(M.mul(W))-CASADI::convert(I);
  for (int i = 0; i < M1.size1(); ++i){
    for (int j = 0; j < M1.size2(); ++j){
	  objfun += mu_stiff*M1.elem(i,j)*M1.elem(i,j);
	  objfun += mu_mass*_Diag_KM.elem(i,j)*_Diag_KM.elem(i,j);
	}
  }

  objfun *= 0.5f;
}

SXMatrix MaterialFitting_Diag_M::assembleObjMatrix(){
  
  const MatrixXd I = MatrixXd::Identity(lambda.size(),lambda.size());
  SXMatrix M1 = trans(W).mul(M.mul(W))-CASADI::convert(I);
  // for (int i = 0; i < M1.size1(); ++i){
  //   for (int j = 0; j < M1.size2(); ++j){
  // 	  if (i!=j)	M1(i,j) *= 0.1f;
  // 	}
  // }
  return M1;
}

void MaterialFitting_Diag_M::initShearG(const int num_tet, VSX &G){
  
  const vector<double> &r = tetmesh->material()._G;
  assert_eq(num_tet, r.size());
  G.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
	G[i] = r[i];
  }
}

void MaterialFitting_Diag_M::initLame(const int num_tet, VSX &Lame){

  const vector<double> &r = tetmesh->material()._lambda;
  assert_eq(num_tet, r.size());
  Lame.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
	Lame[i] = r[i];
  }
}

void MaterialFitting_Diag_M::getInitValue(VectorXd &init_x)const{
    
  const vector<double> &r = tetmesh->material()._rho;
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = r[i];
  }
}

vector<double> MaterialFitting_Diag_M::getDensityResult()const{

  const int num_tet = tetmesh->tets().size();
  assert_ge(rlst.size(),num_tet);
  return vector<double>(rlst.begin(),rlst.begin()+num_tet);
}

SXMatrix MaterialFitting_Diag_K::assembleObjMatrix(){

  const SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  return M1;
}

void MaterialFitting_EV::getInitValue(VectorXd &init_x)const{
  
  const vector<double> &g = tetmesh->material()._G;
  const vector<double> &la = tetmesh->material()._lambda;
  const vector<double> E = MaterialFitting::getYoungE(g,la);
  const vector<double> v = MaterialFitting::getPoissonV(g,la);
  
  const int num_tet = tetmesh->tets().size();
  init_x.resize(num_tet*2);
  for (int i = 0; i < num_tet; ++i){
    init_x[i] = E[i];
    init_x[i+num_tet] = v[i];
  }
}

SXMatrix MaterialFitting_EV::assembleObjMatrix(){

  SXMatrix M1 = trans(W).mul(K.mul(W))-lambda;
  // for (int i = 0; i < M1.size1(); ++i){
  //   for (int j = 0; j < M1.size2(); ++j){
  // 	  if (i!=j)	M1(i,j) *= 1.0f;
  // 	}
  // }
  return M1;
}

SXMatrix MaterialFitting_Diag_M_ScaleW::assembleObjMatrix(){
  
  const SXMatrix S = CASADI::makeEyeMatrix(scale_W);
  const SXMatrix I = CASADI::makeEyeMatrix(VectorXd::Ones(scale_W.size()));

  const SXMatrix M1 = trans(W).mul(M.mul(W))-S;
  const SXMatrix M2 = sqrt(1e-4)*(I-S);

  SXMatrix M3(S.size1()*2, S.size2());
  for (int i = 0; i < S.size1(); ++i)
    for (int j = 0; j < S.size2(); ++j){
	  M3(i,j) = M1(i,j);
	  M3(S.size2()+i,j) = M2(i,j);
	}
  return M3;
}

void MaterialFitting_Diag_M_ScaleW::getInitValue(VectorXd &init_x)const{
  
  const int num_tet = tetmesh->tets().size();
  const int r = W.size2();
  const vector<double> &rho = tetmesh->material()._rho;
  init_x.resize(num_tet+r);
  for (int i = 0; i < num_tet; ++i)
    init_x[i] = rho[i];
  for (int i = 0; i < r; ++i)
    init_x[num_tet+i] = 1.0f;
}

SXMatrix MaterialFitting_EV_MA_K::assembleObjMatrix(){

  const SXMatrix WtMWL = trans(W).mul(M.mul(W)).mul(lambda);
  const SXMatrix M1 = (trans(W).mul(K.mul(W)) - WtMWL);
  return M1;
}

// save input W, lambda, M , fixed nodes, and mesh.vol
bool MaterialFitting::saveAllInputs(const string mesh_name)const{

  TRACE_FUN();

  const int r = W.size2();
  const int nx3 = tetmesh->nodes().size()*3;

  const MatrixXd W = CASADI::convert<double>(this->W);
  bool succ = write(mesh_name+"W.b_norm_is_"+TOSTR(W.norm()), W); assert(succ);

  cout<< "\n\n" << W << "\n\n";

  VectorXd lambda(r);
  for (int i = 0; i < r; ++i){
    lambda[i] = (this->lambda.elem(i,i).getValue());
  }
  succ=write(mesh_name+"lambda_scaled.b_norm_is_"+TOSTR(lambda.norm()),lambda);assert(succ);

  cout<< "\n\n" << lambda << "\n\n";

  VectorXd M_diag(nx3);
  assert_eq(M.size1(), M.size2());
  assert_eq(M.size1(), nx3);
  for (int i = 0; i < nx3; ++i){
    M_diag[i] = M.elem(i,i).getValue();
  }
  succ = write(mesh_name+"M.b_norm_is_"+TOSTR(M_diag.norm()), M_diag); assert(succ);

  cout<< "\n\n" << M_diag << "\n\n";

  vector<int> fixednodes_vec;
  BOOST_FOREACH(const int ele, fixednodes){
	fixednodes_vec.push_back(ele);
  }
  succ = writeVec(mesh_name+"con_nodes.txt",fixednodes_vec,UTILITY::TEXT); assert(succ);
  succ = tetmesh->write(mesh_name+"mesh.vol"); assert(succ);

  return true;
}
