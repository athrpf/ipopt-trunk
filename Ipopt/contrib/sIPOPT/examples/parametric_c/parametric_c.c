

#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x,
	    Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

/* Main Program */
int main()
{
  const Index n = 5;          /* Number of variables */
  const Index m = 4;          /* Number of constraints */
  const Index nnz_jac_g = 10; /* Number of non-zeros in jac_g */
  const Index nnz_h_lag = 5;  /* Number of non-zeros in h_lag */
  const Index index_style = 1; /* Fortran style */

  Number* x_L = NULL;
  Number* x_U = NULL;
  Number* g_L = NULL;
  Number* g_U = NULL;

  IpoptProblem nlp = NULL;
  enum ApplicationReturnStatus status;
  Number* x = NULL;

  Number* mult_g = NULL;               /* constraint multipliers
             at the solution */
  Number* mult_x_L = NULL;             /* lower bound multipliers
             at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
             at the solution */
  Number obj;                          /* objective value */
  int i; /* counter used in several places */

  /* our user data for the function evalutions. */
  struct MyUserData user_data;

  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  x   = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds and starting points*/
  for (i=0; i<3; i++) {
    x_L[i] = 0.0;
    x_U[i] = 1.0e19;
    x[i] = 5.0;
  }
  for (i=3; i<5; i++) {
    x_L[i] = -1.0e-19;
    x_U[i] = 1.0e19;
    x[i] = 1.0;
  }

  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);

  for (i=0; i<2; i++) {
    g_l[i] = 0.0;
    g_u[i] = 0.0;
  }

  /* Set initial value constraints */
  g_l[2] = g_u[2] = nominal_eta1;
  g_l[3] = g_u[3] = nominal_eta2;

  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &eval_f, &eval_g, &eval_grad_f,
                           &eval_jac_g, &eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  /* allocate space to store the bound multipliers at the solution */
  mult_g = (Number*)malloc(sizeof(Number)*m);
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);


  if (status == Solve_Succeeded) {
    status = sIpoptSolve(nlp, ...);

    printf("\n\nSolution of the primal variables, x\n");
    for (i=0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the ccnstraint multipliers, lambda\n");
    for (i=0; i<m; i++) {
      printf("lambda[%d] = %e\n", i, mult_g[i]);
    }
    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i=0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i=0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);
  }
  else {
    printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
  }

  FreeIpoptProblem(nlp);
  free(x);
  free(mult_g);
  free(mult_x_L);
  free(mult_x_U);

  return (int)status;
}
