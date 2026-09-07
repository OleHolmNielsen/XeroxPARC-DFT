#!/bin/bash

# set -x
TMPDIR=.

# Install the GNU gfortran compiler and BLAS libraries.
# CHANGE this for your software environment: A Linux system may use:
# sudo dnf install gcc-gfortran make openblas

BINDIR=.
CRYSTAL=bn
FILE=$CRYSTAL.dat

# The lattice structure file fort.2 (Fortran unit 2) with atomic coordinates
cat <<'EOF' >fort.2
 1-cell Boron Nitride
 2      Atoms in the unit cell
 0.0 0.5 0.5      0.5 0.0 0.5     0.5 0.5 0.0  a1,a2,a3
 5         -0.125    -0.125    -0.125          B atom position
 7          0.125     0.125     0.125          N atom position
 3.61                                          Lattice constant (Angstrom)
EOF
#
# Set proper dimensions:
#
NTYPMX=2
NSPIN=1
NDIM1=7000
NDIM2=500
NDIM3=150
NDIM4=800
NDIM6=32768
NDIM8=4
NDIM9=2
NDIM10=100
NDIM13=$NDIM3
NG1MAX=20
NG2MAX=20
NG3MAX=20
NCMPLX=2
# NCMPLX=1
export NTYPMX NSPIN NDIM1 NDIM2 NDIM3 NDIM4 NDIM6 NDIM8 NDIM9 NDIM13 \
NG1MAX NG2MAX NG3MAX NCMPLX
 
#
RUN290=1
RUN213=1
RUN214=0
RUN207=1

line="============================================================================="

if [[ $RUN290 -eq 1 ]]
then
	PROGRAM=run290
	cat <<EOF
$line

Running $PROGRAM

$line
EOF
	make run290
	./run290 <<'EOF'
0 0 0 0 0    (sphere radius)**2, mesh size na1,na2,na3, and epsilon
0 0          Skip more detail
-1           Continue
2 2 2 0 0 0  IQ1,IQ2,IQ3,WVK0
1            Save to file fort.3
0 0 0 0 0 0  Exit
EOF
fi

if [[ $? -ne 0 ]]
then
	echo "Program $PROGRAM exited with errors"
	exit 1
fi
echo $line
read -p "Press Enter to continue"

# ==================================================================================

if [[ $RUN213 -eq 1 ]]
then
	PROGRAM=run213
	cat <<EOF
$line

Running $PROGRAM

$line
EOF
    # Potential 17 = Bell Labs 1982 "periodic table" (Bachelet-Hamann-Schluter,
    # Phys. Rev. B 26, 4199 (1982)).  Both B and N are in the POT17 parameter
    # database, so no external pseudopotential files are needed.
    # The alternative is potential 20 (numerically generated), which reads
    # Fourier-transformed potentials from B.VG and N.VG - files not present here.
    make run213
    ./run213 <<EOF
-1              No: virtual crystal approximation
0 0 0 0 0       Plane-wave cutoff in Rydbergs, NG1,NG2,NG3, EPSILON
0 0             No details
-1              Continue
17              Atom 1 potential-type
17              Atom 2 potential-type
0.8             Exchange-factor for linear screening of the ionic potential
-1              Use a different number of electrons?
0 0             No details
1.8e-6          Potential clean-up lower limit epsilon
1               Repeat the display of total potential for selected G-vectors
0 0             No details
-1              No clean-up of Potential
1               Save results to binary file fort.10
EOF

if [[ $? -ne 0 ]]
then
	echo "Program $PROGRAM exited with errors"
	exit 1
fi
    echo Move the output file to become input file for K207
    mv fort.10 fort.4
echo $line
read -p "Press Enter to continue"
fi

# ==================================================================================

if [[ $RUN214 = 1 ]]
then  
    # fetch fort.1 -mV2 -fTR -t"DSN=rkmk005.$FILE"
    PROGRAM=run214
	cat <<EOF
$line

Running $PROGRAM

$line
EOF
    make $PROGRAM
    ./$PROGRAM

if [[ $? -ne 0 ]]
then
	echo "Program $PROGRAM exited with errors"
	exit 1
fi
    echo Move the output file to become input file for K207
    mv fort.10 fort.4
echo $line
read -p "Press Enter to continue"
fi

# ==================================================================================

if [[ $RUN207 = 1 ]]
then
    if test $NCMPLX = 1; then
        PROGRAM=run207
    else
        PROGRAM=crun207
    fi

cat <<EOF
$line

Running $PROGRAM

$line
EOF
    make $PROGRAM
    ./$PROGRAM <<'EOF'
1            spec. pts.
1            XC is OK
1            semiconductor
1            # electrons OK
1            # eigenvalues OK
12 18 2      Plane-wave cutoffs E1 (A waves) E2 (Lowdin B waves) in Rydberg units (2)
2            iterative diagonalization
15           max # iterations
1E-7         eigenvalue accuracy
0.5          FAC
2            Cycle when recycling eigenvectors
-1           No VNL file
-1           - nor any reading of VNL from a file
1            Test the dimension of the Hamiltonian (A and B waves)
1            # plane waves are OK
1            Band structure display - NDSPL1,NDSPL2 OK
1            V(G) display OK
-1           Change ISWCH
2            ISWCH: total energy only from bands, no forces or stress (see k207aux.for)
1            Rho(r) display OK
-1           Change FFT dimension
32 32 32     FFT mesh dimension
1            stress
1            forces
10           No. of SCF cycles
-1           no initial Rho(r) guess
-1           Wavefunction projection radius - off
-1           no display of mesh points
-1           no further potential display
2 3 1.0E-8 0.9 0.5 2  IDMAT (dielectric matrix), UPDATE mixing (Broyden), parameters
2            Modify any parameter defining the run
1            Yes, change some parameters
1            E1,E2 OK
1            NDSPL1,NDSPL2 OK
1            V(G) display OK
-1           change ISWCH
1            ISWCH: full calculation with possibly energy forces and stress (see k207aux.for)
1            Rho(r) display OK
-1           change FFT
32 32 32     FFT mesh dimension
1            Calculate stress
1            Calculate forces
1            No. of cycles
-1           no initial Rho(r) guess
-1           projection radius - off
-1           no display of mesh points
-1           no further potential display
2 3 1.0E-8 0.9 0.5 2  IDMAT (dielectric matrix), UPDATE mixing (Broyden), parameters
2            modify any parameter defining the run
-1           do not modify parameters
-1           Band structure
3            # k-points
0 0 0        GAMMA point
1 0 0        X-point
.5 .5 .5     L-point
5            Stop, save the last working potential for Fortran unit 10 (fort.10)
-1           No eigenvalues stored on unit 13 (fort.13)
EOF

if [[ $? -ne 0 ]]
then
	echo "Program $PROGRAM exited with errors"
	exit 1
fi
    # Saving the data file:
    # $HOME/archive fort.10 $CRYSTAL.dat
    # rm EVFILE fort.15 fort.14 fort.13
    # rm $PROGRAM
echo $line
fi
