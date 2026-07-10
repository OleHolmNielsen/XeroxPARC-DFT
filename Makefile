#
# This Makefile is adapted for Linux systems.
# Author: Ole.H.Nielsen@fysik.dtu.dk
# Homepage: https://github.com/OleHolmNielsen/XeroxPARC-DFT/
#
# Maintain the compiled subroutines and executables.
# This version does NOT use libraries
#

# Add file types to the known list of suffixes (file extensions)
.SUFFIXES: .for .start

# Run silently.  Use "set -x" to display commands executed
# .SILENT:
# Ignore errors: trap them yourself !
# .IGNORE:

# File to receive the listings
LISTING = listing

# The fortran compiler:
# FC = f77
# Intel Fortran
# FC = ifx
# GNU Fortran
FC = gfortran
# FFLAGS=-O
FFLAGS=-cpp -Dx86_64 -g -std=legacy --warn-no-argument-mismatch

# The loader
LDR = $(FC) $(FFLAGS)

# Any extra stuff that you want to link in (libraries etc.)
# Generic UNIX:
# OBJS = sciblas.o unixtra.o unixus.o error.o
# UNIX with BLAS library
# OBJS = unixtra.o unixus.o error.o
OBJS =
LIBS =

#
# The default target:
#
none:

LIBS: $(LIBS)

#
# Program K290
#

OBJ_290 = gcode.o     gshell.o   sorting.o  sppt2.o   \
          charge.o   day.o     group1.o  lenstr.o  reclat.o \
          usage.o machine.o  consts.o   structur.o k290.o $(OBJS)
# lookup.o   sorting.o  

run290: run290.o $(OBJ_290)
	$(LDR) -o $@ $@.o $(OBJ_290) $(LIBS)

#
# Program K213
#

OBJ_213 = alphaz.o    eps1.o      gcode.o     gshell.o   \
          sfact.o     sorting.o   atomrdwr.o \
          charge.o   day.o      group1.o   ldftyp.o  lenstr.o   \
          paramete.o potentia.o potredwr.o reclat.o  structur.o \
          machine.o  consts.o   usage.o    erf.o       k213.o 

run213: run213.o $(OBJ_213)
	$(LDR) -o $@ $@.o $(OBJ_213) $(LIBS)

#
# Program K214
#

OBJ_214 = eps1.o      gcode.o     lenstr.o   lookup.o    \
          consts.o   charge.o   day.o      potredwr.o \
          machine.o  usage.o    k214.o 

run214: run214.o $(OBJ_214)                                    \
    ; set -x                                                  \
    ; $(LDR) -o $@ $@.o $(OBJ_214) $(LIBS)

#
# Program K207 (real version)
#

# OBJ_207 = alphaz.o    apwsum.o     cfft.o    dblas.o     compit.o   \
# IBM AIX version with BLAS built in:
OBJ_207 = alphaz.o    apwsum.o     cfft.o    dsdot.o    compit.o   \
          diagon.o    diaham.o     efermi.o  eigen.o     eisherm.o  \
          eps1.o      eval.o       exch4.o   gcode.o     gkcut.o    \
          gshell.o    ionion.o     kinetic.o lookup.o    lowper.o   \
          lumat.o     matel.o      projct.o  putnl.o     realit.o   \
          ro5.o       rosym4.o     rowham.o  rpwsum.o    setup.o    \
          sfact.o     sorting.o    unicon.o  utilmem.o   vnlkka.o   \
          vnlkknum.o  chkinv.o                           \
          atomrdwr.o charge.o   day.o     excorr.o   forces.o  \
          gaussq.o   ldftyp.o   lenstr.o  mixvg.o    k207aux.o \
          potget.o   rwev.o     strloc.o  potentia.o potredwr.o\
          consts.o   tsep.o     usage.o    vnlred.o  k207.o    \
          machine.o  spptrd.o   erf.o

run207: run207.o $(OBJ_207)
	$(LDR) -o $@ $@.o $(OBJ_207) $(LIBS)

# Apollo specials (owing to compiler bugs):
apollo_special_1: 
	$(MAKE) FFLAGS=-g k207.o ionion.o mixvg.o


#
# Program K207 (complex version)
#

OBJ_207c = alphaz.o    apwsum.o     cfft.o     compit.o    dblas.o    \
           cdiagon.o   cdiaham.o    efermi.o   eigen.o     eisherm.o  \
           eps1.o      ceval.o      exch4.o    gcode.o     gkcut.o    \
           gshell.o    ionion.o     ckinetic.o lookup.o    clowper.o  \
           clumat.o    matel.o      projct.o   cputnl.o    realit.o   \
           cro5.o      rosym4.o     crowham.o  rpwsum.o    setup.o    \
           sfact.o     sorting.o    unicon.o   utilmem.o   cvnlkka.o  \
           cvnlkknum.o chkinv.o                            \
           atomrdwr.o charge.o   day.o      excorr.o   forces.o  \
           gaussq.o   ldftyp.o   lenstr.o   cmixvg.o   k207aux.o \
           potget.o   crwev.o    strloc.o   potentia.o potredwr.o\
           consts.o   tsep.o     usage.o    vnlred.o   k207.o    \
           machine.o  spptrd.o   erf.o
 
crun207: crun207.o $(OBJ_207c)  \
    ; set -x                                       \
    ; $(LDR) -o $@ $@.o $(OBJ_207c) $(LIBS)

# Apollo specials (owing to compiler bugs):
apollo_special_2: 
	$(MAKE) FFLAGS=-g k207.o ionion.o cmixvg.o
 

# The object files
objects: $(OBJ_290) $(OBJ_213) $(OBJ_207) $(OBJ_207c)

#
# Create main program Fortran files with correct dimensions
#
 
# Force the preocessing of these main programs:
run290.for run213.for run214.for run207.for crun207.for: FRC
FRC:

utilmem.for: FRC
 
# Rule to build main programs:
# Replace dimension parameters NDIM?? etc. by desired values from
# the file dimensions.sh
.start.for:                                      \
    ; echo Building $@ from $< ..             \
    ; . ./dimensions.sh                                \
    ; if test -z "$$NTYPMX"                      \
    ; then echo Error: dimensions are undefined  \
    ;      exit 2                                \
    ; fi                                         \
    ; export TMP=`mktemp`               \
    ; echo "."                                \
    ; echo "s/\$$NTYPMX/$$NTYPMX/g" >>$$TMP      \
    ; echo "s/\$$NDIM13/$$NDIM13/g" >>$$TMP      \
    ; echo "s/\$$NSPIN/$$NSPIN/g"   >>$$TMP      \
    ; echo "s/\$$NDIM1/$$NDIM1/g"   >>$$TMP      \
    ; echo "s/\$$NDIM2/$$NDIM2/g"   >>$$TMP      \
    ; echo "s/\$$NDIM3/$$NDIM3/g"   >>$$TMP      \
    ; echo "s/\$$NDIM4/$$NDIM4/g"   >>$$TMP      \
    ; echo "s/\$$NDIM6/$$NDIM6/g"   >>$$TMP      \
    ; echo "s/\$$NDIM8/$$NDIM8/g"   >>$$TMP      \
    ; echo "s/\$$NDIM9/$$NDIM9/g"   >>$$TMP      \
    ; echo "s/\$$NDIM10/$$NDIM10/g" >>$$TMP      \
    ; echo "s/\$$NG1MAX/$$NG1MAX/g" >>$$TMP      \
    ; echo "s/\$$NG2MAX/$$NG2MAX/g" >>$$TMP      \
    ; echo "s/\$$NG3MAX/$$NG3MAX/g" >>$$TMP      \
    ; echo "s/\$$NCMPLX/$$NCMPLX/g" >>$$TMP      \
    ; echo "."                                \
    ; rm -f $@                                   \
    ; sed -f $$TMP $< >$@                        \
    ; echo "."                                \
    ; rm -f $$TMP                                \
    ; echo " Done."


#
# Create the subroutines that handle complex Hamiltonians
# using the public-domain "patch" program
#

cdiagon.for: diagon.for cdiagon.diff
	echo Building $@ using $?
	patch <cdiagon.diff diagon.for -o $@
cdiaham.for: diaham.for cdiaham.diff
	echo Building $@ using $?
	patch <cdiaham.diff diaham.for -o $@
ceval.for: eval.for ceval.diff
	echo Building $@ using $?
	patch <ceval.diff eval.for -o $@
ckinetic.for: kinetic.for ckinetic.diff
	echo Building $@ using $?
	patch <ckinetic.diff kinetic.for -o $@
clowper.for: lowper.for clowper.diff
	echo Building $@ using $?
	patch <clowper.diff lowper.for -o $@
clumat.for: lumat.for clumat.diff
	echo Building $@ using $?
	patch <clumat.diff lumat.for -o $@
cmixvg.for: mixvg.for cmixvg.diff
	echo Building $@ using $?
	patch <cmixvg.diff mixvg.for -o $@
cputnl.for: putnl.for cputnl.diff
	echo Building $@ using $?
	patch <cputnl.diff putnl.for -o $@
cro5.for: ro5.for cro5.diff
	echo Building $@ using $?
	patch <cro5.diff ro5.for -o $@
crowham.for: rowham.for crowham.diff
	echo Building $@ using $?
	patch <crowham.diff rowham.for -o $@
crun207.for: run207.for crun207.diff
	echo Building $@ using $?
	patch <crun207.diff run207.for -o $@
crwev.for: rwev.for crwev.diff
	echo Building $@ using $?
	patch <crwev.diff rwev.for -o $@
cvnlkka.for: vnlkka.for cvnlkka.diff
	echo Building $@ using $?
	patch <cvnlkka.diff vnlkka.for -o $@
cvnlkknum.for: vnlkknum.for cvnlkknum.diff
	echo Building $@ using $?
	patch <cvnlkknum.diff vnlkknum.for -o $@

#
# The datafiles
#

C.VG Si.VG B.VG N.VG:                              \
    ;- set -xv                                     \
    ; FILE=$@                                      \
    ; fetch $@ -mV2 -fTR -t"DSN=rkmk005.$$FILE"    \
    ; if test ! -f $@ ; then                       \
        echo "File $@ not found"                   \
    ;   exit 2                                     \
    ; fi

#
# Subroutines necessary for the "stressfield" program
#
stressfield: alphaz.o atomrdwr.o cfft.o charge.o chkinv.o clgor.o   \
	consts.o day.o efermi.o exch4.o excorr.o forces.o   \
	gkcut.o ionion.o k207aux.o ldftyp.o lenstr.o lookup.o       \
	potentia.o potget.o potredwr.o ro5.o rosym4.o rwev.o       \
	sciblas.o sfact.o sorting.o spptrd.o unixtra.o unixus.o       \
	machine.o  usage.o vnlred.o erf.o

testeis: $$@.o eigen.o eisherm.o unixus.o \
		unixtra.o day.o usage.o machine.o lenstr.o
	set -x; $(FC) $(FFLAGS) -o $@ $@.o eigen.o eisherm.o unixus.o \
		unixtra.o day.o usage.o machine.o lenstr.o

testion: testion.o ionion.o usage.o unixus.o \
		unixtra.o day.o machine.o lenstr.o consts.o erf.o
	set -x; $(FC) $(FFLAGS) -o testion \
		testion.o ionion.o usage.o \
		unixus.o unixtra.o day.o machine.o \
		lenstr.o consts.o erf.o

#
# Maintenance
#
clean:
	rm -f *.f *.bak *.orig *~ run290 run290.for run213 run207 crun207 core listing
	rm -f *.o

#
# Compilation rules (Generic UNIX)
#

.for.o:
	$(FC) $(FFLAGS) -c $*.for

FLINT_FLAGS = -gftspa
RUN290_SRC = run290.for gcode.for     gshell.for   sorting.for sppt2.for   \
	charge.for   day.for     group1.for  lenstr.for  reclat.for \
	machine.for  consts.for   structur.for usage.for   k290.for

RUN213_SRC = run213.for alphaz.for    eps1.for      gcode.for     gshell.for   \
	sfact.for    sorting.for  atomrdwr.for \
	charge.for   day.for      group1.for   ldftyp.for  lenstr.for   \
	paramete.for potentia.for potredwr.for reclat.for  structur.for \
	machine.for  consts.for   usage.for    erf.for       k213.for 

RUN207_SRC = run207.for \
	alphaz.for   apwsum.for     cfft.for    dblas.for     compit.for   \
	diagon.for   diaham.for     efermi.for  eigen.for     eisherm.for  \
	eps1.for     eval.for       exch4.for   gcode.for     gkcut.for    \
	gshell.for   ionion.for     kinetic.for lookup.for    lowper.for   \
	lumat.for    matel.for      projct.for  putnl.for     realit.for   \
	ro5.for      rosym4.for     rowham.for  rpwsum.for    setup.for    \
	sfact.for    sorting.for    unicon.for  utilmem.for   vnlkka.for   \
	vnlkknum.for chkinv.for                           \
	atomrdwr.for charge.for   day.for     excorr.for   forces.for  \
	gaussq.for   ldftyp.for   lenstr.for  mixvg.for    k207aux.for \
	potget.for   rwev.for     strloc.for  potentia.for potredwr.for\
	consts.for   tsep.for     usage.for    vnlred.for  k207.for    \
	machine.for  spptrd.for   erf.for \
	sciblas.for  unixtra.for 

RUN207C_SRC = run207.for \
	alphaz.for    apwsum.for     cfft.for     compit.for    dblas.for    \
	cdiagon.for   cdiaham.for    efermi.for   eigen.for     eisherm.for  \
	eps1.for      ceval.for      exch4.for    gcode.for     gkcut.for    \
	gshell.for    ionion.for     ckinetic.for lookup.for    clowper.for  \
	clumat.for    matel.for      projct.for   cputnl.for    realit.for   \
	cro5.for      rosym4.for     crowham.for  rpwsum.for    setup.for    \
	sfact.for     sorting.for    unicon.for   utilmem.for   cvnlkka.for  \
	cvnlkknum.for chkinv.for                            \
	atomrdwr.for charge.for   day.for      excorr.for   forces.for  \
	gaussq.for   ldftyp.for   lenstr.for   cmixvg.for   k207aux.for \
	potget.for   crwev.for    strloc.for   potentia.for potredwr.for\
	consts.for   tsep.for     usage.for    vnlred.for   k207.for    \
	machine.for  spptrd.for   erf.for \
	sciblas.for  unixtra.for 

flint: flint_run290

flint_run290:
	echo Doing flint on RUN290_SRC
	flint $(FLINT_FLAGS) $(RUN290_SRC)

flint_run213:
	echo Doing flint on RUN213_SRC
	flint $(FLINT_FLAGS) $(RUN213_SRC)

flint_run207:
	echo Doing flint on RUN207_SRC
	flint $(FLINT_FLAGS) $(RUN207_SRC)

flint_run207c:
	echo Doing flint on RUN207C_SRC
	flint $(FLINT_FLAGS) $(RUN207C_SRC)
