/*
  Large Scale Computing

  solve Ax=b  
  Conjugate Gradient Method
  General Sparse Matrix A

  solve_cg solve Ax=b with CG method
  this returns the number of iterations, or negative if failed.

  Note: A must be symmetric and positive-defined matrix
 */
#ifndef __CG_OMP_H__
#define __CG_OMP_H__

#include<iostream>
#include<cmath>
#include <omp.h>

int solve_cg_OMP(int dim,          /* dimension */
	     int nnz,       /* # of non-zeros in the matrix */
	     double *nzval,    /* array of nonzero values */
	     int *colidx,      /* array of column indices of the nonzeros */
	     int *rowptr,      /* the first column of nonzero in each row */
	     /* rowptr[j] stores the location of nzval[] and colidx[], */
	     /* which starts row j. This has dim+1 entry that is nnz */
	     double *x,        /* solition */
	     double *b,        /* right hand side */
	     double tol);      /* tolerance */


#endif // __CG_OMP_H__
