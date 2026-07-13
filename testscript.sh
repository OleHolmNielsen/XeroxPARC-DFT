#!/bin/bash

# set -x
TMPDIR=.

# Trap error signals:
trap "rm -rf $TMPDIR; exit 2" 1 2 3 14 15 19

BINDIR=.
CRYSTAL=bn
FILE=$CRYSTAL.dat

# The structure file with atomic coordinates
cat <<'EOF' >fort.2
 1-cell Boron Nitride
 2
 0.0 0.5 0.5      0.5 0.0 0.5     0.5 0.5 0.0
 5         -0.125    -0.125    -0.125
 7          0.125     0.125     0.125
 3.00
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
NDIM13=$NDIM3
NG1MAX=20
NG2MAX=20
NG3MAX=20
NCMPLX=2
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
0 0 0 0 0
0 0
-1
2 2 2 0 0 0
1
0 0 0 0 0 0
EOF
fi

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
-1
0 0 0 0 0
0 0
-1
17
17
0.8
-1
0 0
1.8e-6
1
0 0
-1
1
EOF
    # Move the output file to become input file for K207
    mv fort.10 fort.4
fi

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
    # Move the output file to become input file for K207
    mv fort.10 fort.4
fi

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
6 18 2       cutoffs (Ry)
2            iterative diagonalization
15           max # iterations
1E-7         eigenvalue accuracy
0.5          FAC
2            Cycle when recycling eigenvectors
-1           No VNL file
-1           - nor any reading of VNL
1            test dim Hamiltonian
1            - OK
1            NDSPLi OK
1            V(G) display OK
-1           change ISWCH
    2        new ISWCH
1            Rho(r) display OK
-1           change FFT
32 32 32     FFT
1            stress
1            forces
    4        # cycles
-1           no initial Rho(r) guess
-1           projection radius - off
-1           no display of mesh points
-1           no further potential display
2 3 1.0E-8 0.9 0.5 2      mixing parameters
2            switch
1            modify parameters
1            E1,E2 OK
1            NDSPLi OK
1            V(G) display OK
-1           change ISWCH
    1        new ISWCH
1            Rho(r) display OK
-1           change FFT
32 32 32     FFT
1            stress
1            forces
    1        # cycles
-1           no initial Rho(r) guess
-1           projection radius - off
-1           no display of mesh points
-1           no further potential display
2 3 1.0E-8 0.9 0.5 2      mixing parameters
2            switch
-1           do not modify parameters
-1           Band structure
3            # k-points
0 0 0        GAMMA point
1 0 0        X-point
.5 .5 .5     L-point
5            Stop, keep potential
-1           No eigenvalues kept
EOF
    # Saving the data file:
    # $HOME/archive fort.10 $CRYSTAL.dat
    # rm EVFILE fort.15 fort.14 fort.13
    # rm $PROGRAM
fi
