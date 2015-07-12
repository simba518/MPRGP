#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <MoseckQPSolver.h>

BOOST_AUTO_TEST_SUITE(MoseckQPSolverTest)

BOOST_AUTO_TEST_CASE(testQPSolver){

  const int num_var = 3;
  const int num_con = 4;

  SparseMatrix<double> Q, A;
  VectorXd b, c, x;

  {// init constraints: A*x >= c
	A.resize(num_con, num_var);
	A.reserve(6);
	for (int i = 0; i < 3; ++i){
	  A.insert(0,i) = 1;
	  A.insert(i+1,i) = 1;
	}
	A.makeCompressed();

	c.resize(num_con);
	c.setZero();
	c[0] = 1;
  }

  {// init objective function: Q, b
	Q.resize(num_var, num_var);
	Q.reserve(5);
	Q.insert(0,0) = 2.0;
	Q.insert(1,1) = 0.2;
	Q.insert(2,2) = 2.0;
	Q.insert(2,0) = -1.0;
	Q.insert(0,2) = -1.0;

	b.resize(num_var);
	b.setZero();
	b[1] = -1.0;
  }

  {// print the qp info
	// const MatrixXd Am = A;
	// const MatrixXd Qm = Q;
	// cout<< "Q: \n" << Qm << endl << endl;
	// cout<< "b: " << b.transpose() << endl;
	// cout<< "A: \n" << Am << endl << endl;
	// cout<< "c: " << c.transpose() << endl;
  }

  {// test QP with constraints.
	x.resize(num_var);
	x.setZero();
	MoseckQPSolver solver;
	solver.setConstraints(A, c, num_var);
	TEST_ASSERT( solver.solve(Q,b,x) );

	// check results
	VectorXd correct_x(3);
	correct_x << 5.97500593165e-05, 5.00000000451, 5.97500593165e-05;
	ASSERT_EQ_SMALL_VEC_TOL(correct_x, x, 3, 1e-10);
  }

  {// test QP without constraints.

	x.resize(num_var);
	x.setZero();

	MoseckQPSolver solver;
	A.setZero();
	c.setZero();
	solver.setConstraints(A, c, num_var);
	TEST_ASSERT( solver.solve(Q,b,x) );

	// check results
	ASSERT_LE((Q*x + b).norm(), 1e-8);
  }

  {// test QP without constraints.

	x.resize(num_var);
	x.setZero();

	MoseckQPSolver solver;
	A.resize(0,0);
	c.resize(0);
	solver.setConstraints(A, c, num_var);
	TEST_ASSERT( solver.solve(Q,b,x) );

	// check results
	ASSERT_LE((Q*x + b).norm(), 1e-8);
  }
  
}

static void MSKAPI printstr_test(void *handle, MSKCONST char str[]){
  // printf("%s",str);
}

BOOST_AUTO_TEST_CASE(testMoseck){

  // example download from: http://docs.mosek.com/7.1/capi/Quadratic_optimization.html
  const int NUMCON = 1;   /* Number of constraints.             */ 
  const int NUMVAR = 3;   /* Number of variables.               */ 
  // const int NUMANZ = 3;   /* Number of non-zeros in A.           */ 
  const int NUMQNZ = 4;   /* Number of non-zeros in Q.           */ 

  double        c[]   = {0.0,-1.0,0.0};

  MSKboundkeye  bkc[] = {MSK_BK_LO};
  double        blc[] = {1.0};
  double        buc[] = {+MSK_INFINITY}; 
    
  MSKboundkeye  bkx[] = {MSK_BK_LO,
                         MSK_BK_LO,
                         MSK_BK_LO};
  double        blx[] = {0.0,
                         0.0,
                         0.0};
  double        bux[] = {+MSK_INFINITY,
                         +MSK_INFINITY,
                         +MSK_INFINITY};
  
  MSKint32t     aptrb[] = {0,   1,   2}, aptre[] = {1,   2,   3}, asub[]  = {0,   0,   0};
  double        aval[]  = {1.0, 1.0, 1.0};
  
  MSKint32t     qsubi[NUMQNZ];
  MSKint32t     qsubj[NUMQNZ];
  double        qval[NUMQNZ];
  
  MSKint32t     i,j;
  double        xx[NUMVAR];

  MSKenv_t      env = NULL;
  MSKtask_t     task = NULL;
  MSKrescodee   r;
  
  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);

  if ( r==MSK_RES_OK ){ 

	  /* Create the optimization task. */
	  r = MSK_maketask(env,NUMCON,NUMVAR,&task);

	  if ( r==MSK_RES_OK )
		{
		  r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr_test);
  
		  /* Append 'NUMCON' empty constraints.
			 The constraints will initially have no bounds. */
		  if ( r == MSK_RES_OK )
			r = MSK_appendcons(task,NUMCON);
  
		  /* Append 'NUMVAR' variables.
			 The variables will initially be fixed at zero (x=0). */
		  if ( r == MSK_RES_OK )
			r = MSK_appendvars(task,NUMVAR);
  
		  /* Optionally add a constant term to the objective. */
		  if ( r ==MSK_RES_OK )
			r = MSK_putcfix(task,0.0);
		  for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
			{
			  /* Set the linear term c_j in the objective.*/  
			  if(r == MSK_RES_OK)
				r = MSK_putcj(task,j,c[j]);
  
			  /* Set the bounds on variable j.
				 blx[j] <= x_j <= bux[j] */
			  if(r == MSK_RES_OK)
				r = MSK_putvarbound(task,
									j,           /* Index of variable.*/
									bkx[j],      /* Bound key.*/
									blx[j],      /* Numerical value of lower bound.*/
									bux[j]);     /* Numerical value of upper bound.*/
  
			  /* Input column j of A */   
			  if(r == MSK_RES_OK)
				r = MSK_putacol(task,
								j,                 /* Variable (column) index.*/
								aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
								asub+aptrb[j],     /* Pointer to row indexes of column j.*/
								aval+aptrb[j]);    /* Pointer to Values of column j.*/
        
			}
  
		  /* Set the bounds on constraints.
			 for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
		  for(i=0; i<NUMCON && r==MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
								i,           /* Index of constraint.*/
								bkc[i],      /* Bound key.*/
								blc[i],      /* Numerical value of lower bound.*/
								buc[i]);     /* Numerical value of upper bound.*/

		  if ( r==MSK_RES_OK )
			{
			  /*
			   * The lower triangular part of the Q
			   * matrix in the objective is specified.
			   */

			  qsubi[0] = 0;   qsubj[0] = 0;  qval[0] = 2.0;
			  qsubi[1] = 1;   qsubj[1] = 1;  qval[1] = 0.2;
			  qsubi[2] = 2;   qsubj[2] = 0;  qval[2] = -1.0;
			  qsubi[3] = 2;   qsubj[3] = 2;  qval[3] = 2.0;

			  /* Input the Q for the objective. */

			  r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);
			}

		  if ( r==MSK_RES_OK )
			{
			  MSKrescodee trmcode;
        
			  /* Run optimizer */
			  r = MSK_optimizetrm(task,&trmcode);

			  /* Print a summary containing information
				 about the solution for debugging purposes*/
			  MSK_solutionsummary (task,MSK_STREAM_MSG);
        
			  if ( r==MSK_RES_OK )
				{
				  MSKsolstae solsta;
          
				  MSK_getsolsta (task,MSK_SOL_ITR,&solsta);
          
				  switch(solsta)
					{
					case MSK_SOL_STA_OPTIMAL:   
					case MSK_SOL_STA_NEAR_OPTIMAL:
					  MSK_getxx(task,
								MSK_SOL_ITR,    /* Request the interior solution. */
								xx);
					  ASSERT_EQ_TOL(xx[0], 5.975006e-05, 1e-5);
					  ASSERT_EQ_TOL(xx[1], 5.000000e+00, 1e-5);
					  ASSERT_EQ_TOL(xx[2], 5.975006e-05, 1e-5);
					  // printf("Optimal primal solution\n");
					  // for(int j=0; j<NUMVAR; ++j)
					  // 	printf("x[%d]: %e\n",j,xx[j]);
              
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
				}
			  else
				{
				  printf("Error while optimizing.\n");
				}
			}
    
		  if (r != MSK_RES_OK)
			{
			  /* In case of an error print error code and description. */      
			  char symname[MSK_MAX_STR_LEN];
			  char desc[MSK_MAX_STR_LEN];
        
			  printf("An error occurred while optimizing.\n");     
			  MSK_getcodedesc (r,
							   symname,
							   desc);
			  printf("Error %s - '%s'\n",symname,desc);
			}
		}
	  MSK_deletetask(&task);
	}

  MSK_deleteenv(&env);

}

BOOST_AUTO_TEST_SUITE_END()
