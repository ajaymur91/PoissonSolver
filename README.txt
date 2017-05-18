PoissonSolver v1.0 (with OpenMP support)
-------------------------------------------------
When to use this Package?
* You have a trajectory/frame from a Molecular Dynamics package and want to calculate electrostatic potentials in the system. 

How to use this Package?
* Convert the MD equilibrated trajectory or (single frame) to crdbox format and name it coord.crd
* Supply input information through inputs_density.txt and inputa_potential.txt (mdp style file,self explanatory)
* Run the submission script "jobsc". alter this file to suit your batch server system (or computer)..

Disclaimer:
* This is the preliminary build and requires plenty of documentation/clean up. 
* At this point the most reusable part of the code is solving an AX=B system using conjugate gradients optimization. 
* If you are interested in using this for your research or would like to develop this further, catch me at (ajaymur91 at gmail dot com)

-------------------------------------------------
Submission script : jobsc (View jobsc to understand how the code works together)
-------------------------------------------------
Inputs 
*inputs_density.txt (mdp style file,self expl)
*inputs_potential.txt (mdp style file,self expl)
*coord.crd (converted gromacs trajectory to crdbox format)

Outputs
*Potential.dat by solving Poisson (PBC) using Conjugate gradient method
in the format 
x,y,z,Potential 
with units (nm,nm,nm,Volts)
Solves (24*24*684) mesh size in about 10-15 sec
-------------------------------------------------


-------------------------------------------------
Main Code: 
-------------------------------------------------
*cg_omp.cpp (solve AX=B)
*cg_omp.h (header defines the calling routine name)
*density.cpp (Calculate Matrix B) 
*potential.cpp (set up Matrix A) and call Conjugate gradient routine
-------------------------------------------------
Density script (density.cpp) - Calculate B Matrix
-------------------------------------------------
* Input 1: coord.crd 
  coord.crd format:  X1	Y1 Z1 X2 Y2 Z2................Xn Yn Zn Lx Ly Lz
  
* Input 2: input_density.txt (Open input_density.txt for details)
* Read Location of atoms from the simulation box (coord.crd)
* n is number of atoms, Lx, Ly and Lz are box lengths
* Calculate the charge densities on a 3D grid, assigning partial charges to the atoms from PC.txt, BF4.txt, TEA.txt, CNTN_*e.txt, CNTP_*e.txt
* Prints output to density.txt 
* Specific to the CNT system, low reusability

-------------------------------------------------
Potential script (potential.cpp) - set up Matrix A)
-------------------------------------------------
* Input 1: density.txt (output from previous step) 
* Input 2: input_potential.txt (Open input_potential.txt for details)
* Solve Poisson Equation to calculate potentials in the simulation box (Periodic Boundary Conditions)
* Convert differential eqn to linear difference equations on a 3d grid to get AX=B.
* Use Conjugate gradient method to solve AX=B(Periodic Boundary Conditions)
* High reusability
-------------------------------------------------
