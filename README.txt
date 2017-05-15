Submission script : jobsc

-------------------------------------------------
Code 
-------------------------------------------------
*density.cpp (Calculate Matrix B) 
*potential.cpp (set up Matrix A) and call Conjugate gradient routine
*cg_omp.cpp (solve AX=B)
*cg_omp.h (header defines the calling routine name)

-------------------------------------------------
Density script (density.cpp) - Calculate B Matrix
-------------------------------------------------
* Read Location of atoms from the simulation box (coord.crd)
* coord.crd format:  X1	Y1 Z1 X2 Y2 Z2................Xn Yn Zn Lx Ly Lz
* n is number of atoms, Lx, Ly and Lz are box lengths
 * Calculate the charge densities on a 3D grid, assigning partial charges to the atoms from PC.txt, BF4.txt, TEA.txt, CNTN_*e.txt, CNTP_*e.txt
* Specific to the CNT system, low reusability
-------------------------------------------------
Potential script (potential.cpp) - set up Matrix A)
-------------------------------------------------
*Solve Poisson Equation to calculate potentials in the simulation box (Periodic Boundary Conditions)
*Convert differential eqn to linear difference equations on a 3d grid to get AX=B.
*Use Conjugate gradient method to solve AX=B(Periodic Boundary Conditions)
*high reusability

-------------------------------------------------
Inputs 
-------------------------------------------------
*inputs_density.txt (mdp style file,self expl)
*inputs_potential.txt (mdp style file,self expl)
*coord.crd (converted gromacs trajectory to crdbox format)

-------------------------------------------------
Outputs
-------------------------------------------------
*Potential.dat by solving Poisson (PBC) using Conjugate gradient method
in the format 
x,y,z,Potential 
with units (nm,nm,nm,Volts)
Solves (24*24*684) mesh size in about 10-15 sec
