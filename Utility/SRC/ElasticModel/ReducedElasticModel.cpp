#include "ReducedElasticModel.h"
#include <MassMatrix.h>
#include <JsonFilePaser.h>
using namespace UTILITY;
using namespace SIMULATOR;

bool ReducedElasticModel::init(const string filename){
  
  JsonFilePaser jsonf;
  bool succ = false;
  if (jsonf.open(filename)){
  	succ = jsonf.readMatFile("subspace_basis",B);
  }else{
  	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

bool DirectReductionElasticModel::prepare(){

  bool succ = false;
  ElasticForceTetFullStVK::prepare();
  if (_vol_mesh){

	succ = true;
	_vol_mesh->nodes(rest_shape);
	assert_gt(rest_shape.size(),0);

	UTILITY::MassMatrix mass;
	Eigen::DiagonalMatrix<double,-1> diag_M;
	mass.compute(diag_M, *_vol_mesh);
	assert_eq(diag_M.rows(),B.rows());
	M = B.transpose()*diag_M*B;
  }
  return succ;
}

bool CubaturedElasticModel::init(const string filename){

  bool succ = ReducedElasticModel::init(filename);
  if (!succ){
  	return false;
  }

  JsonFilePaser jsonf;
  if (jsonf.open(filename)){
  	jsonf.readVecFile("cubature_weights",weights);
  	jsonf.readVecFile("cubature_points",sampledTets,UTILITY::TEXT);
	if (weights.size() <=0 || sampledTets.size() <= 0){
	  weights.clear();
	  sampledTets.clear();
	}
  }else{
  	ERROR_LOG("failed to open the initfile: " << filename);
  }
  return succ;
}

bool CubaturedElasticModel::prepare(){

  bool succ = false;
  ElasticForceTetFullStVK::prepare();
  if (_vol_mesh){

	succ = true;
	_vol_mesh->nodes(rest_shape);
	assert_gt(rest_shape.size(),0);

  	UTILITY::MassMatrix mass;
  	Eigen::DiagonalMatrix<double,-1> diag_M;
  	mass.compute(diag_M, *_vol_mesh);
  	assert_eq(diag_M.rows(),B.rows());
  	M = B.transpose()*diag_M*B;

	if (weights.size() <=0 || sampledTets.size() <= 0){

	  const int n = _vol_mesh->tets().size();
	  weights.resize(n);
	  sampledTets.resize(n);
	  for (int i = 0; i < n; ++i){
		weights[i] = 1.0f;
		sampledTets[i] = i;
	  }
	}
  }
  return succ;
}

// f = sum wi*Bi^t*fi(q), for i in S.
bool CubaturedElasticModel::evaluateF(const VectorXd &reduced_u,VectorXd &f){
  
  assert_eq(weights.size(), sampledTets.size());
  assert_eq(rest_shape.size(), B.rows());
  assert_eq(reduced_u.size(), B.cols());
  const VectorXd x = rest_shape+B*reduced_u;

  f.resize(reducedDim());
  f.setZero();

  static mat3x4 f_tet;
  for (size_t i = 0; i < sampledTets.size(); ++i){
	const int tet_id = sampledTets[i];
    force_tet(f_tet, tet_id, x);
	for (int j = 0;  j < 4; ++j){
	  const int vi = _vol_mesh->tets()[tet_id][j];
	  f += B.block(3*vi,0,3,B.cols()).transpose()*((-1.0f*weights[i])*f_tet.col(j));
	}
  }
  return true;
}

// K = sum wi*Bi^t*Ki(q)*Bi, for i in S.
bool CubaturedElasticModel::evaluateK(const VectorXd &reduced_u,MatrixXd &K){
  
  assert_eq(weights.size(), sampledTets.size());
  assert_eq(rest_shape.size(), B.rows());
  assert_eq(reduced_u.size(), B.cols()); 
  const VectorXd x = rest_shape+B*reduced_u;

  K.resize(reducedDim(),reducedDim());
  K.setZero();

  TetDF K_tet;
  for (size_t i = 0; i < sampledTets.size(); ++i){
    forceDerivX_tet(K_tet, sampledTets[i], x);
	for (int j = 0; j < 4; ++j){
	  const int vj = _vol_mesh->tets()[i][j];
	  const Matrix<double,3,-1> &Bj = B.block(3*vj,0,3,B.cols());
	  for (int k = 0; k < 4; ++k){
		const int vk = _vol_mesh->tets()[i][k];
		const Matrix<double,3,-1> &Bk = B.block(3*vk,0,3,B.cols());
		const Matrix3d Ajk = (-1.0f*weights[i])*K_tet.df[j][k];
		K += Bj.transpose()*Ajk*Bk;
	  }
	}
  }
  return true;
}
