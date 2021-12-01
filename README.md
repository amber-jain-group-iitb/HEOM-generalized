# heom_amber
Codes for running hierarichal equation of motion. 
Strumpfer Schulten, JCTC 8, 2808 (2012). 
Help also taken from  Chen, Zheng, Shi, Tan, JCP 131, 094502 (2009)

To run, type make (uses ifort). Run with ./aout. Parallelized the code using openmp. Should take a few seconds. Creates a file rho.out that has the diagonal entries of the density matrix as a function of time.

A simpler (ans slower) compilation is:
gfortran mod_spectra_heom.f90 spectra_heom.f90 -fopenmp
./a.out
