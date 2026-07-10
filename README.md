# Xerox PARC (Palo Alto Research Center) DFT plane-wave code

This historical project preserves the DFT plane-wave code written originally
by Karel Kunc and Richard M. Martin at the Xerox PARC (Palo Alto Research Center) during 1979-1981.

The code was further developed by Ole Holm Nielsen at Xerox PARC during 1982-1983.

The code implements Total Energy, Forces, and Stresses
with a pseudopotential plane-wave basis set.

The code is written in Fortran-77, the most modern Fortran compiler at the time.

A run of the code consists of three sequential code steps:

1. The run290 code analyses the crystal structure and generates
   Special k-points and symmetry operations.

2. The run213 code sets up the pseudopotential.

3. The run207 code makes a self-consistent DFT (Density Functional Theory)
   calculation of Total Energy, Forces, and Stresses.

Input files:

* TBD

Output files:

* TBD
