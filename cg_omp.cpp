//
//  Large Scale Computing
//
//  solve Ax=b  
//  Conjugate Gradient Method
//  General Sparse Matrix A
//  OpenMP Version
//
//  Note: A must be symmetric and positive-defined matrix
//
#include "cg_omp.h"

int solve_cg_OMP(int dim, /* dimension */
			int nnz, /* # of non-zeros in the matrix */
			double *nzval, /* array of nonzero values */
			int *colidx, /* array of column indices of the nonzeros */
			int *rowptr, /* the first column of nonzero in each row */
			double *x, /* solition */
			double *b, /* right hand side */
			double tol) { /* tolerance */

	/* Allocate Working Spaces */
	double *r = new double[dim];
	double *p = new double[dim];
	double *Ap = new double[dim];

	int csi = std::max(10,dim / omp_get_max_threads() / 10);
	//std::cout << csi << std::endl;
	/* compute r0 */
#pragma omp parallel for schedule(guided,csi)
	for (int i = 0; i < dim; i++) {
		r[i] = b[i];
		for (int j = rowptr[i]; j < rowptr[i + 1]; j++) {
			r[i] -= nzval[j] * x[colidx[j]];
		}
	}

	double rr = 0.0;
#pragma omp parallel for reduction(+:rr) schedule(guided,csi)
	for (int i = 0; i < dim; i++) {
		p[i] = r[i]; /* p = r */
		rr += r[i] * r[i]; /* rr = r.r */
	}
	double nRhs = rr;

	if (nRhs == 0.0) {
		//for (int i = 0; i < dim; i++) x[i] = 0.0;
		return 0;
	}

	/* cg iteration */
	int count = 0;
	double pAp;
	double rr1;
	while (rr > tol * tol * nRhs) {

		// Ap = A*p
#pragma omp parallel for schedule(guided,csi)
		for (int i = 0; i < dim; i++) {
			Ap[i] = 0.0;
			for (int j = rowptr[i]; j < rowptr[i + 1]; j++) {
				Ap[i] += nzval[j] * p[colidx[j]];
			}
		}

		// alpha = r.r / p.Ap
		pAp = 0.0;
#pragma omp parallel for reduction(+:pAp) schedule(guided,csi)
		for (int i = 0; i < dim; i++) {
			pAp += p[i] * Ap[i];
		}
		double alpha = rr / pAp;

		//Beta
		rr1 = 0.0;
#pragma omp parallel for reduction(+:rr1) schedule(guided,csi)
		for (int i = 0; i < dim; i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			rr1 += r[i] * r[i];
		}
		double beta = rr1 / rr;

#pragma omp parallel for schedule(guided,csi)
		for (int i = 0; i < dim; i++) {
			p[i] = r[i] + beta * p[i];
		}
		rr = rr1;
		count++;
	}

	/* Deallocate Working Spaces */
	delete[] r;
	delete[] p;
	delete[] Ap;

	return count;
}
