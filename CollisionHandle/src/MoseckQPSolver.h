#ifndef _MOSECKQPSOLVER_H_
#define _MOSECKQPSOLVER_H_

// #include </home/simba/mosek/7/tools/platform/linux64x86/h/mosek.h>
#include <mosek.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
using namespace std;
using namespace Eigen;

/**
 * @class MoseckQPSolver an class providing the interface for solving QP problem using moseck.
 * 
 */
class MoseckQPSolver{
	
public:
  MoseckQPSolver(){
	debug = false;
	env = NULL;
	task = NULL;
  }
  bool setConstraints(const SparseMatrix<double,ColMajor> &A, const VectorXd &c,const int NUMVAR); //A*x >= c
  bool solve(const SparseMatrix<double> &Q, const VectorXd &b, VectorXd &x);
  void clear();
  ~MoseckQPSolver(){
	clear();
  }

protected:
  bool setObjectiveFunc(const SparseMatrix<double> &Q, const VectorXd &b);
  bool optimize(VectorXd &x);
	
private:
  bool debug;
  MSKrescodee r;
  MSKenv_t env;
  MSKtask_t task;

  vector<double> qval;
  vector<int> qsubi, qsubj;
};

#endif /*_MOSECKQPSOLVER_H_*/
