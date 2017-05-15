
//  Conjugate Gradient Method
//  Input density on 3d grid - Solve 3D discretized Poisson Equation (Periodic Boundaries) and output potentials
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cg_omp.h"   // use cg parallel

int main(int argc, char **argv) {

//inputs
        int numx,numy,numz,BD,Nf;
	double dx, tol;
	std::ifstream ifile1;
        ifile1.open("inputs_potential.txt");
	ifile1 >> numx >> numy >> numz >> dx >> tol >> Nf >> BD;
	ifile1.close();
	std::ifstream ifile;
        ifile.open("density.txt");
//calculate conversion factor for RHS (SI units)
 double F = 1.602 * dx * dx * 1000/(8.854 * BD * Nf);
 std::cout<<"\nConversion factor = "<<F<<"\n";
//copy density onto rhs array, multiply by dx*dx and convert to SI units
	double *den = new double[(numx+1) * (numy+1) * (numz+1)] ();
	double *rhs = new double[(numx) * (numy) * (numz)] ();

	   for (int k = 0; k <= numz; k++) {
                for (int j = 0; j <= numy; j++) {
                        for(int i=0; i<= numx;i++) {
                   	ifile >> den[i + (numx+1) * (j+(numy+1)*k)] ;
                }
        }
}
           for (int k = 0; k < numz; k++) {
                for (int j = 0; j < numy; j++) {
                        for(int i=0; i< numx;i++) {
                        rhs[i + (numx) * (j+(numy)*k)]=den[i + (numx+1) * (j+(numy+1)*k)]*(-F);
			//-F * den = RHS
                }
        }
}
delete []den;


	/// allocate memory space
	double *phi = new double[numx * numy * numz] ();
	double *sol = new double[numx * numy * numz] (); 
//	double *rhs = new double[numx * numy * numz] ();
	double *nzval = new double[7 * numx * numy * numz] ();
	int *colind = new int[7 * numx * numy * numz] ();
	int *rowptr = new int[numx * numy *numz+ 1] ();


	/* Initialize */
	for (int n = 0; n < numx * numy * numz; n++) {
		phi[n] = 0.0;
	}

	double start = clock();
	/* Setup sparse Matrix A and RHS */
	int nnz = 0;
      for(int k=0; k<numz;k++) {
	for (int j = 0; j < numy ; j++) {
		for (int i = 0; i < numx ; i++) {

			int n = i + numx * (j+numy*k); // index of grid points
			int nd = i + (numx+1) * (j+(numy+1)*k); // index of vector elements

			/* initial guess */
			sol[n] = phi[n];

			/* set first nonzero column of row index(i,j) */
			rowptr[n] = nnz;
//27 cases to set up sparse matrix A for periodic case 
if (k == 0){
                if(j == 0){
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                 nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                 nzval[nnz] = 1.0; colind[nnz] = n + (numx)*(numy-1); nnz++;
                        }
                else if (j == numy-1){
                nzval[nnz] = 1.0; colind[nnz] = n - (numx)*(numy-1); nnz++;
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                                }
                else    {
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                        }

        nzval[nnz] = 1.0; colind[nnz] = n + numx*numy; nnz++;
        nzval[nnz] = 1.0; colind[nnz] = n + numx*numy*(numz-1); nnz++;
        }
else if (k == numz-1){
        nzval[nnz] = 1.0; colind[nnz] = n - numx*numy*(numz-1); nnz++;
        nzval[nnz] = 1.0; colind[nnz] = n - numx*numy; nnz++;
                if(j == 0){
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                nzval[nnz] = 1.0; colind[nnz] = n + (numx)*(numy-1); nnz++;
                        }
                else if (j == numy-1){
                nzval[nnz] = 1.0; colind[nnz] = n - (numx)*(numy-1); nnz++;
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                                }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                                }
                else    {
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                }

        }
else    {
        nzval[nnz] = 1.0; colind[nnz] = n - numx*numy; nnz++;
                if(j == 0){
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                                }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                nzval[nnz] = 1.0; colind[nnz] = n + (numx)*(numy-1); nnz++;
                }
                else if (j == numy-1){
                nzval[nnz] = 1.0; colind[nnz] = n - (numx)*(numy-1); nnz++;
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                                }
                else    {
                nzval[nnz] = 1.0; colind[nnz] = n - numx; nnz++;
                                if (i == 0){
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + numx - 1; nnz++;
                                        }
                                else if (i == numx-1){
                                nzval[nnz] = 1.0; colind[nnz] = n - (numx-1); nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                        }
                                else    {
                                nzval[nnz] = 1.0; colind[nnz] = n - 1; nnz++;
                                nzval[nnz] = -6.0; colind[nnz] = n; nnz++;
                                nzval[nnz] = 1.0; colind[nnz] = n + 1; nnz++;
                                        }
                nzval[nnz] = 1.0; colind[nnz] = n + numx; nnz++;
                        }

        nzval[nnz] = 1.0; colind[nnz] = n + numy*numx; nnz++;
        }

//Setup sparse matrix for periodic case end
		}
	}
}

	/* last element of rowptr si nnz */
	rowptr[(numx) * (numy)* (numz)] = nnz;

	/* solve with CG */
//	const double tol = 1.0e-4;
	int k = solve_cg_OMP((numx) * (numy)* (numz), nnz, nzval, colind, rowptr, sol, rhs, tol);
	if (k < 0) {
		std::cout << "calculation failed\n";
		return 0;
	}
	double tcost = (clock() - start) / CLOCKS_PER_SEC;
	std::cout << "# of Iteration=" << k << std::endl;
	std::cout << "Time cost (CPU) = " << tcost << "(sec)\n";

	// restore
for(int k=0; k<numz ; k++) {
	for (int j = 0; j < numy ; j++) {
		for (int i = 0; i < numx ; i++) {
			int n = i + numx *( j+numy*k); // index of grid points
			//int nn = (i - 1) + (num - 2) *( (j - 1)+(num-2)*(k-1)); // index of vector elements
			phi[n] = sol[n];
		}
	}
}
	// Output Result
	std::ofstream ofile;
	ofile.open("potential.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
ofile << "x,y,z,v" << std::endl;
	for (int i = 0; i < numx; i++) {
		for (int j = 0; j < numy; j++) {
                 	for(int k=0; k<numz;k++) {
			//double x = -1.0 + 2.0 * i / (num - 1);
			//double y = -1.0 + 2.0 * j / (num - 1);
			//double z = -1.0 + 2.0 * k / (num - 1);
			ofile << i*0.05 << "," << j*0.05 << "," << k*0.05 << "," << phi[i + numx * (j+numy*k)] << std::endl;
		}
	}
}

	ofile.close();
	ifile.close();
	delete []phi;
	delete []nzval;
	delete []colind;
	delete []rowptr;
	delete []sol;
	delete []rhs;	
	std::cout << "Done\n";
	return 0;
}
