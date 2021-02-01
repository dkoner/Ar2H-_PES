# Ar<sub>2</sub>H<sup>+</sup> PES
Global 3D potential energy surface (PES) for Ar<sub>2</sub>H<sup>+</sup> constructed using CCSD(T)/aug-cc-pVQZ energies
 
**Requirements**

(1) fortran90 compiler

**Compile the PES**

A test program file (test.f90) is given and it can be compiled as

`gfortran -O3 Ar2H+.f90 test.f90`

**Running the executable**

Before running the executable make sure that the coeff.dat file is present in the current directory (or change the file path in the fortran program).

**Cite as**

D. Koner J. Chem. Phys. 154, 054303 (2021)
