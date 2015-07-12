#include <stdio.h> 
#include <Log.h>
#include "MoseckQPSolver.h"

static void MSKAPI printstr(void *handle, MSKCONST char str[]){
  // printf("%s",str);
}

//set A and c for: A*x >= c
bool MoseckQPSolver::setConstraints(const SparseMatrix<double>&A,const VectorXd &c, const int NUMVAR){

  assert(A.isCompressed());
  clear();
  const int NUMCON = A.rows();

  r = MSK_makeenv(&env,NULL);
  assert( r == MSK_RES_OK );

  /* Create the optimization task. */
  r = MSK_maketask(env,NUMCON,NUMVAR,&task);
  assert( r == MSK_RES_OK );

  r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
  assert ( r == MSK_RES_OK );
  
  /* Append 'NUMCON' empty constraints. The constraints will initially have no bounds. */
  r = MSK_appendcons(task, NUMCON);
  assert ( r == MSK_RES_OK );
  
  /* Append 'NUMVAR' variables. The variables will initially be fixed at zero (x=0). */
  r = MSK_appendvars(task, NUMVAR);
  
  /* Optionally add a constant term to the objective. */
  assert ( r ==MSK_RES_OK );
  r = MSK_putcfix(task,0.0);

  /* Set the bounds on variable j. blx[j] <= x_j <= bux[j] */
  for(int j=0; j<NUMVAR && r == MSK_RES_OK; ++j){
	assert(r == MSK_RES_OK);
	r = MSK_putvarbound(task,
						j,           /* Index of variable.*/
						MSK_BK_LO,   /* Bound key.*/
						-MSK_INFINITY,      /* Numerical value of lower bound.*/
						+MSK_INFINITY);     /* Numerical value of upper bound.*/
  }
  assert(r == MSK_RES_OK);

  /* Input column j of A */   
  if(A.nonZeros() > 0){
	const double *aval = A.valuePtr();
	const int *asub = A.innerIndexPtr();
	const int *aptrb = A.outerIndexPtr();
	const int *aptre = A.outerIndexPtr()+1;
	for(int j=0; j < NUMVAR && r == MSK_RES_OK; ++j){
	  assert(r == MSK_RES_OK);
	  r = MSK_putacol(task,
					  j,                 /* Variable (column) index.*/
					  aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
					  asub+aptrb[j],     /* Pointer to row indexes of column j.*/
					  aval+aptrb[j]);    /* Pointer to Values of column j.*/
	}
	assert(r == MSK_RES_OK);
  }
  
  /* Set the bounds on constraints. for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
  for(int i=0; i<NUMCON && r == MSK_RES_OK; ++i){
	r = MSK_putconbound(task,
						i,           /* Index of constraint.*/
						MSK_BK_LO,      /* Bound key.*/
						c[i],        /* Numerical value of lower bound.*/
						+MSK_INFINITY);     /* Numerical value of upper bound.*/
  }
  assert(r == MSK_RES_OK);

  return (r == MSK_RES_OK);
}

bool MoseckQPSolver::solve(const SparseMatrix<double> &Q, const VectorXd &b, VectorXd &x){

  x.resize(b.size());
  setObjectiveFunc(Q, b);
  optimize(x);
  return (r == MSK_RES_OK);
}

void MoseckQPSolver::clear(){

  if (env)
	MSK_deleteenv(&env);
}

bool MoseckQPSolver::setObjectiveFunc(const SparseMatrix<double> &Q, const VectorXd &b){

  assert ( r == MSK_RES_OK );
  const int NUMVAR = b.size();

  // set b
  for(int j=0; j<NUMVAR && r == MSK_RES_OK; ++j){
	/* Set the linear term b_j in the objective.*/  
	if(r == MSK_RES_OK)
	  r = MSK_putcj(task,j,b[j]);
  }
  assert ( r == MSK_RES_OK );

  // set Q

  /*
   * The lower triangular part of the Q
   * matrix in the objective is specified.
   */
  const int NUMQNZ = Q.nonZeros();
  qval.clear();
  qsubi.clear();
  qsubj.clear();
  qval.reserve(NUMQNZ);
  qsubi.reserve(NUMQNZ);
  qsubj.reserve(NUMQNZ);

  for (int k = 0; k < Q.outerSize(); ++k){
	for(const typename SparseMatrix<double>::InnerIterator it(Q, k);it;++it){
	  const int r = it.row();
	  const int c = it.col();
	  if( r >= c ){
		qsubi.push_back(r);
		qsubj.push_back(c);
		qval.push_back(it.value());
	  }
	}
  }

  /* Input the Q for the objective. */
  r = MSK_putqobj(task,qval.size(),&qsubi[0],&qsubj[0],&qval[0]);
  assert ( r == MSK_RES_OK );

  return (r == MSK_RES_OK);
}

bool MoseckQPSolver::optimize(VectorXd &x){

  assert (x.size() > 0);
  assert ( r == MSK_RES_OK );
        
  /* Run optimizer */
  MSKrescodee trmcode;
  r = MSK_optimizetrm(task, &trmcode);

  /* Print a summary containing information about the solution for debugging purposes*/
  if (debug)
	MSK_solutionsummary (task, MSK_STREAM_MSG);

  // get solution
  if ( r == MSK_RES_OK ){

	MSKsolstae solsta;
	MSK_getsolsta (task,MSK_SOL_ITR,&solsta);

	switch(solsta) {

	case MSK_SOL_STA_OPTIMAL:
	case MSK_SOL_STA_NEAR_OPTIMAL:/* Request the interior solution. */
	  MSK_getxx(task, MSK_SOL_ITR, &x[0]);
	  break;

	case MSK_SOL_STA_DUAL_INFEAS_CER:
	case MSK_SOL_STA_PRIM_INFEAS_CER:
	case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
	case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
	  printf("Primal or dual infeasibility certificate found.\n");
	  break;
              
	case MSK_SOL_STA_UNKNOWN:
	  printf("The status of the solution could not be determined.\n");
	  break;

	default:
	  printf("Other solution status.");
	  break;
	}

  }else{

	/* In case of an error print error code and description. */
	char symname[MSK_MAX_STR_LEN];
	char desc[MSK_MAX_STR_LEN];
        
	printf("An error occurred while optimizing.\n");     
	MSK_getcodedesc (r, symname, desc);
	printf("Error %s - '%s'\n",symname, desc);
  }

  return (r == MSK_RES_OK);
}
