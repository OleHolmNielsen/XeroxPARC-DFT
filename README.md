# Xerox PARC (Palo Alto Research Center) DFT plane-wave code

Purpose of this project
===========================

This historical project preserves the [Density Functional Theory](https://en.wikipedia.org/wiki/Density_functional_theory) (DFT)
[plane-wave](https://en.wikipedia.org/wiki/Plane_wave) 
[pseudopotentials](https://en.wikipedia.org/wiki/Density_functional_theory#Pseudo-potentials)
code written originally
by [Karel Kunc](http://www-ext.impmc.upmc.fr/~kunc/) (<krl.kunc@gmail.com>)
and [Richard M. Martin](https://people.physics.illinois.edu/rmartin/) (<rmartin@illinois.edu>)
working at the Xerox [PARC](https://en.wikipedia.org/wiki/PARC_(company)) (*Palo Alto Research Center*) during 1979-1981.

This new scientific development started with two papers by Wendel and Martin cited below.
Kunc and Martin developed in the following papers the
self-consistent DFT calculation methods for obtaining *Total Energy* and *Forces*.
Nielsen and Martin developed the theory of the *Stress tensor*.
Minimization of the *Total Energy* for obtaining crystal structures etc. 
using various optimization algorithms was to some extent a hand-held process at the time
where computers were orders of magnitude slower than today.

The present historical project does *not* aim to give a complete coverage of the scientific field,
nor does it attempt to survey subsequent developments and the works of other groups.
For such an overview the reader is referrred to Richard Martin's book
cited in the *Publications* section.

It must be emphasized that this code was used for cutting-edge research during the 1980ies,
but is today more like a [time capsule](https://en.wikipedia.org/wiki/Time_capsule)
which preserves the code for posterity.
The present code is at the very root of plane-wave pseudopotential calculations,
and many newer codes may have been inspired by this work.

Author of this project: [Ole Holm Nielsen](https://dcwww.fysik.dtu.dk/~ohnielse/) (<Ole.H.Nielsen@fysik.dtu.dk>).

Project homepage: https://github.com/OleHolmNielsen/XeroxPARC-DFT/

Developers of the code
---------------------------

The original code was designed and written by Kunc and Martin at PARC.
It was further developed by Ole Holm Nielsen at PARC during 1982-1983,
and later at [NORDITA](https://nordita.org/) and [University of Copenhagen](https://nbi.ku.dk/english/) until about 1990.

[Richard Needs](https://www.phy.cam.ac.uk/profile/prof-richard-needs/),
working in the [Volker Heine](https://www.phy.cam.ac.uk/profile/prof-volker-heine/) group at the
[Cavendish Laboratory](https://www.phy.cam.ac.uk/) of [University of Cambridge](https://www.cam.ac.uk/),
continued the work at Xerox PARC and Cavendish after Ole Holm Nielsen.

Later on there have been other users of the code at multiple institutions, particularly in the 
[Richard Martin](https://physics.illinois.edu/people/directory/profile/rmartin) group
at [University of Illinois Urbana-Champaign](https://illinois.edu/).
Probably many of these users have made independent improvements to the code over the years,
leading to unavoidable divergences and multiple versions of the code.
This all happened before the [Internet](https://en.wikipedia.org/wiki/Internet) became ubiquitous,
and software [version control](https://en.wikipedia.org/wiki/Version_control)
changed the way we maintain computer codes.

Publications
------------------

Richard Martin's book
[Electronic Structure: Basic Theory and Practical Methods, 2nd Edition](https://www.cambridge.org/dk/universitypress/subjects/physics/condensed-matter-physics-nanoscience-and-mesoscopic-physics/electronic-structure-basic-theory-and-practical-methods-2nd-edition)
constitutes an authoritative overview of the field up to 2020,
covering also historical methods.

The following publications from the PARC group document the theory and computational methods
behind the original codes in the present project,
and they highlight some of the first and groundbreaking scientific and computational methods created in this field.
As is well-known, DFT calculations have become a sprawling scientific endeavour in numerous fields of science,
and they occupy large amounts of supercomputer time world-wide.

In chronological order the fundamental papers from the PARC group are:

* *Charge Density and Structural Properties of Covalent Semiconductors*,
  H. Wendel and Richard M. Martin,
  [Phys. Rev. Lett. 40, 950, 1978](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.40.950).
  Presents a practical scheme to calculate the charge density and total energy
  of crystals as a function of atomic displacements.

* *Theory of structural properties of covalent semiconductors*,
  H. Wendel and Richard M. Martin,
  [Phys. Rev. B 19, 5251, 1979](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.19.5251).
  These papers laid the foundations for the PARC DFT codes for solving the pseudopotential
  [Hamiltonian](https://en.wikipedia.org/wiki/Hamiltonian_(quantum_mechanics))
  [eigenvalue problem](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors) 
  and obtain the total energy within a plane-wave basis set,
  but the calculations were not fully self-consistent. 

* *Relaxation of Polar Ge-GaAs (100) Interfaces: Self-Consistent Calculations of Total Energy*,
  K. KUNC, R. M. MARTIN,
  in [Proc. 15th Int. Conf. Phys. Semiconductors, Kyoto, 1980 (S. Tanaka et al., eds.)](https://www.jps.or.jp//books/jpsjs/49A/jpsj.49sa-toc.html):
  [J. Phys. Soc. Japan 49 - Suppl. A, 1117 - 1120 (1980)](https://www.jps.or.jp//books/jpsjs/49A/jpsj.49sa.1117.pdf).
  This paper is the first publication of fully self-consistent calculations made with the PARC code.
  The calculations optimized the structure using the total energy and also derived band offsets.  

* *Ab initio Force Constants of Germanium*. K. KUNC, R. M. MARTIN,
  in *Proc. 1st Int. Conf. Phonon Physics, Bloomington, Indiana, 1981 (W. E. Bron, ed.)*:
  [J. Phys. (Paris) 42 - Suppl. C6, 649 - 651 (1981)](https://jphyscol.journaldephysique.org/articles/jphyscol/abs/1981/06/jphyscol198142C6189/jphyscol198142C6189.html)
  ([DOI](https://doi.org/10.1051/jphyscol:19816189)).
  This paper presents interatomic force constants for calculating phonon spectra.
  To our knowledge it is the first published use of
  [Hellmann-Feynman forces](https://en.wikipedia.org/wiki/Hellmann%E2%80%93Feynman_theorem)
  evaluated from DFT calculations.

* *Density-functional calculation of static and dynamic properties of GaAs*,
  K. Kunc and Richard M. Martin,
  [Phys. Rev. B 24, 2311(R), 1981](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.24.2311).
  To our knowledge this paper presents the first "direct" calculations of structure and lattice dynamics of a semiconductor
  using local forms of the ionic pseudopotentials and the density functional.

* *Ab Initio Force Constants of GaAs: A New Approach to Calculation of Phonons and Dielectric Properties*,
  K. Kunc and Richard M. Martin,
  [Phys. Rev. Lett. 48, 406, 1982](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.406).
  It is shown that self-consistent calculations of the electronic charge density in large periodic cells
  containing a single displaced atom provide all the information needed for ab initio determination
  of force constants, phonon dispersion curves, effective charges, and the static dielectric constant. 

The theory of the quantum mechanical *Stress Tensor* was developed in these papers:

* *First-Principles Calculation of Stress*, O. H. Nielsen and Richard M. Martin,
  [Phys. Rev. Lett. 50, 697, 1983](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.50.697)
  This is the first derivation of the *Stress Theorem* for calculating the 
  macroscopic stress tensor in quantum mechanics, and this theory is applied in DFT calculations.

* *Stresses in semiconductors: Ab initio calculations on Si, Ge, and GaAs*,
  O. H. Nielsen and Richard M. Martin
  [Phys. Rev. B 32, 3792, 1985](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.32.3792).
  This paper has explicit formulas for the *Stress Tensor* terms in DFT
  when expressed in a plane-wave basis set.

For calculating [static dielectric properties](https://en.wikipedia.org/wiki/Relative_permittivity)
there was a series of papers by Kunc, Resta and Tosatti:

* *External fields in the Self-Consistent Theory of Electronic states:  A new method for direct evaluation of Macroscopic and Microscopic dielectric response*,
  K. Kunc and R. Resta,
  [Phys. Rev. Letters 51, 686 (1983)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.51.686).

* *Direct evaluation of the inverse dielectric matrix in semiconductors*,
  K. Kunc and E. Tosatti,
  [Phys. Rev. B 29, 7045 (R) (1984)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.29.7045).

* *Self-consistent theory of electronic states and dielectric response in semiconductors*,
  R. Resta, K. Kunc,
  [Phys. Rev. B34, 7146 (1986)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.34.7146).

Other parallel developments
-----------------------------------

To our knowledge the main other development of similar DFT codes
was done independently in the [Marvin Cohen](https://en.wikipedia.org/wiki/Marvin_L._Cohen) 
group at [University of California, Berkeley](https://en.wikipedia.org/wiki/University_of_California,_Berkeley).
The following selected papers reported DFT calculations:

* *Momentum-space formalism for the total energy of solids*,
  J Ihm, A Zunger and M L Cohen,
  [J. Phys. C: Solid State Phys. 12 4409, 1979](https://iopscience.iop.org/article/10.1088/0022-3719/12/21/009).

* *Microscopic Theory of the Phase Transformation and Lattice Dynamics of Si*,
  M. T. Yin and Marvin L. Cohen, [Phys. Rev. Lett. 45, 1004, 1980 ](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.45.1004).

* *Theory of static structural properties, crystal stability, and phase transformations: Application to Si and Ge*,
  M. T. Yin and Marvin L. Cohen,
  [Phys. Rev. B 26, 5668, 1982](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.5668).

Prerequisites to pseudopotential DFT calculations
======================================================

We list some important prerequisite methods which were available for DFT calculations in 1979-1982.

Pseudopotentials 
-------------------------------

State-of-the-art [pseudopotentials](https://en.wikipedia.org/wiki/Density_functional_theory#Pseudo-potentials)
available in 1979-1982 are implemented in the code module
[K213](k213.for) (*Preparation of an input file with pseudopotential for running K207*) including the
following forms of pseudopotentials:
[Appelbaum-Hamann](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.8.1777),
[Berkeley ionic](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.1738),
[Hamann-Schlüter-Chiang norm-conserving pseudopotentials](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.43.1494),
and Bell Labs [Pseudopotentials that work: From H to Pu](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.4199).

To set things in perspective, what we **did not have back then** include more modern pseudopotentials including for example
[Kleinman-Bylander pseudopotentials](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.1425),
[Ultra-soft pseudopotentials](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.7892) (USPP),
and [Projector augmented wave method](https://en.wikipedia.org/wiki/Projector_augmented_wave_method) (PAW).

Exchange-correlation (XC) energy functionals
---------------------------------------------------

In 1979-1982 very few [Exchange-correlation](https://en.wikipedia.org/wiki/Local-density_approximation)
energy functionals were in existence!
Besides the old [Slater X-alpha](https://journals.aps.org/pr/abstract/10.1103/PhysRev.81.385) method
and [Wigner interpolation](https://journals.aps.org/pr/abstract/10.1103/PhysRev.46.1002),
the [K207](k207.for) module implements the
[Ceperley-Alder](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.45.566)
correlation calculated by a [Quantum Monte Carlo](https://en.wikipedia.org/wiki/Quantum_Monte_Carlo) method,
which was by far the most accurate functional available in 1980.

To set things in perspective, what we **did not have back then** include modern
[Exchange-correlation](https://en.wikipedia.org/wiki/Local-density_approximation) methods 
which are described in Chapters 8 and 9 of Richard Martin's book listed in *Publications* above.

Crystal symmetry analysis
----------------------------

The [Crystal structure](https://en.wikipedia.org/wiki/Crystal_structure) is very important,
and symmetry analysis can reduce significantly the amount of computations for, e.g.,
the covalent semiconductors which we mainly studied with the present DFT code.

The lattice symmetry analysis is implemented in the [K290](k290.for) module
(*Analyses the crystal structure and generates Special k-points and symmetry operations*)
which generates input files for the [K207](k207.for) module.

The main [K207](k207.for) DFT calculation module is distinguished
by working with **any kind of lattice symmetry**.
It can work with either *Real* or *Complex*
[Fortran number](https://www.tutorialspoint.com/fortran/fortran_numbers.htm)
matrices and wavefunctions as required by the lattice symmetry.

To our knowledge the ability of the [K207](k207.for) module for 
dealing with low-symmetry lattices was unique in the DFT field back in 1979-1982.

Reciprocal lattice integration by "special k-points"
-----------------------------------------------------------

Integration over the [Reciprocal lattice](https://en.wikipedia.org/wiki/Reciprocal_lattice)
[Brillouin zone](https://en.wikipedia.org/wiki/Brillouin_zone) (BZ)
can be approximated by a weighted sum over a discrete set of "special k-points". 

In the present code, special k-points in the BZ for an arbitrary lattice follows the method of
[Monkhorst and Pack, *Phys. Rev. B* **13** (1976) 5188](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.13.5188)
and [MacDonald, *Phys. Rev. B* **18** (1978) 5897](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.18.5897).
This is implemented in the [K290](k290.for) module.

Symmetry operations are then applied to the "special k-points" 
in the *irreducible* BZ so that integrals cover the entire first BZ.

About the code
=====================

The code in this historical project implements
[Density Functional Theory](https://en.wikipedia.org/wiki/Density_functional_theory) (DFT)
calculations of *Total Energy*, *Forces*, and *Stresses* using 
[pseudopotentials](https://en.wikipedia.org/wiki/Density_functional_theory#Pseudo-potentials)
in a [plane-wave](https://en.wikipedia.org/wiki/Plane_wave) basis set.

The code is written in [Fortran-77](https://en.wikipedia.org/wiki/Fortran),
which was the most modern standard Fortran compiler at the time of writing.
Some quirks in the code are due to limitations of [Fortran-77](https://en.wikipedia.org/wiki/Fortran):

* Arrays had to be allocated with static dimensions in the [main program entry-point](https://en.wikipedia.org/wiki/Entry_point)
  (in the ``run2xx.for`` files) since dynamically allocatable arrays were not yet in the Fortran standard.
  Therefore the scripts must define static dimension variables such as ``NDIM1, NDIM2, NDIM3``, etc.,
  to generate ``run2xx.for`` files from the ``run2xx.start`` files [run290.start](run290.start),
  [run213.start](run213.start), and [run207.start](run207.start).

* The code works both with matrices of either *Real* or *Complex*
  [Fortran number](https://www.tutorialspoint.com/fortran/fortran_numbers.htm) 
  (the former working much faster than the latter)
  depending on the symmetry properties of the crystal.
  Therefore a number of ``xxx.diff`` files are used to generate the *Complex* code version 
  from the *Real* version.
  The Fortran compiler will likely issue warnings about ``passed COMPLEX(4) to REAL(4)``
  (or similar) due to the "dirty" programming tricks employed in the code.
  Remember that nice code [Preprocessors](https://en.wikipedia.org/wiki/Preprocessor)
  did not exist at the time.

Downloading and building the code
---------------------------------------

You can download a copy of the code from the [GitHub](https://en.wikipedia.org/wiki/GitHub)
repository [main page](https://github.com/OleHolmNielsen/XeroxPARC-DFT/tree/main) 
by clicking on the green ``<> Code`` button near the top.
If you have installed the [git command](https://www.kernel.org/pub/software/scm/git/docs/git-clone.html)
the download is simply:
```
git clone https://github.com/OleHolmNielsen/XeroxPARC-DFT.git
```

We have tested the code on a [Linux](https://en.wikipedia.org/wiki/Linux)
Intel-based [x86-64](https://en.wikipedia.org/wiki/X86-64) system.
Some standard [Linux](https://en.wikipedia.org/wiki/Linux) tools are needed for building the code. 
The user will have to configure the computer with these software tools (or similar):

* A [Fortran-77](https://en.wikipedia.org/wiki/Fortran) compatible compiler.
  For example, in 2026 we can compile the code using the
  Open Source [GNU Fortran](https://en.wikipedia.org/wiki/GNU_Fortran),
  but other compilers may also work.
  There is a [GFortran installation guide](https://fortran-lang.org/learn/os_setup/install_gfortran/).

* The Linux [make command](https://en.wikipedia.org/wiki/Make_(software)) 
  together with the [Makefile](Makefile) designed for GNU Fortran compilation
  is used to build object modules and executable codes.
  The [patch](https://en.wikipedia.org/wiki/Patch_(Unix)) command is used here.

* The [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms)
  library for linear algebra operations is required.
  An Open Source library is [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS).

Most Linux distributions can install these tools using their [Package Manager](https://en.wikipedia.org/wiki/Package_manager) system.
For an RPM-based Linux distribution such as *Red Hat Enterprise Linux* or *AlmaLinux*
you can install the following packages:
```
sudo dnf install git gcc-gfortran make patch openblas
```
On Ubuntu or similar Linux distributions the required packages can be installed like this:
```
sudo apt install git make build-essential patch gfortran libopenblas-dev
```

Running the set of codes
=============================

A run of the code consists of several sequential code steps *K290*, *K213*, and *K207*
described in the sections below.
An example of a script running the codes is in the file [testscript.sh](testscript.sh).
This script uses the Linux/UNIX [make](https://en.wikipedia.org/wiki/Make_(software))
command to build the codes, then runs them with the input data given in the script.

Conventions for Fortran unit names and array dimensions
-------------------------------------------------------------

In Fortran files are managed through [logical units identifiers](https://docs.oracle.com/cd/E19957-01/805-4940/6j4m1u7oj/index.html)
which can be from 0 to 99,
and the files on disk are named as ``fort.NN`` where ``NN`` is the unit number.

The file [consts.for](consts.for) defines Fortran unit numbers as well as 
a number of [physical constants](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.41.375)
used throughout the code.
```
2 ...... INPUT OF CRYSTAL STRUCTURE
3 ...... INPUT OF SPECIAL POINTS INFORMATION
4 ...... INPUT OF POTENTIAL FOR SELFCONSISTENCY
5 ...... ON-LINE INPUT
6 ...... ON-LINE OUTPUT
7 ...... ON-LINE OUTPUT, LESS IMPORTANT INFORMATION
8 ...... ON-LINE OUTPUT, EVEN LESS IMPORTANT INFORMATION
10 ..... OUTPUT OF POTENTIAL FOR SELFCONSISTENCY
11 ..... FILE FOR TEMPORARY STORAGE OF NONLOCAL MATRIX ELEMENTS
12 ..... BRIEF SUMMARY OF OUTPUT (FOR LONG-DISTANCE COMPUTING)
20+I ... I=1,2,...NTYPMX UNITS ASSOCIATED WITH ATOMIC FILES
```
The **array dimensions** used must be specified in the file [dimensions.sh](dimensions.sh)
which is used by the [Makefile](Makefile) 
to define static dimension variables such as ``NDIM1, NDIM2, NDIM3``, etc.,
and generate ``run2xx.for`` files from the ``run2xx.start`` files [run290.start](run290.start),
[run213.start](run213.start), and [run207.start](run207.start).

Typical values and meaning of the **array dimensions** are as follows:
```
NCMPLX=1 for real matrices, 2 for complex ones
NTYPMX=2 Max number of atomic types
NSPIN=1 FOR NON-POLARIZED, 2 FOR SPIN-POLARIZED
NDIM1=7000 VARIABLE DIMENSION OF THE ARRAYS IGLIST AND G2LIST (NUMBER OF THE GENERATED POINTS IN REC. SPACE)
NDIM2=500  VARIABLE DIMENSION OF THE TABLES NTAB, G2TAB
NDIM3=150
NDIM4=800
NDIM6=32768
NDIM8=4
NDIM9=2
NDIM13=same as NDIM3
NG1MAX=20
NG2MAX=20
NG3MAX=20
```

Lattice structure file
--------------------------

Firstly, create a structure file (Fortran unit 2 file ``fort.2``) describing the
[Bravais lattice](https://en.wikipedia.org/wiki/Bravais_lattice)
and the atoms in the unit cell:

* Lattice vectors a1, a2, a3
* Atomic numbers and coordinates
* Unit cell size (Angstrom)

for example:

```
 1-cell Boron Nitride
 2
 0.0 0.5 0.5      0.5 0.0 0.5     0.5 0.5 0.0
 5         -0.125    -0.125    -0.125
 7          0.125     0.125     0.125
 3.00
```

K290: Analyses the crystal structure and generates Special k-points and symmetry operations
-----------------------------------------------------------------------------------------------

History:

- Written on September 12, 1979.
- IBM-retouched on October 27, 1980.
- Generation of special points modified on May 26, 1982 by Ole Holm Nielsen.

Purpose:

- Determine the symmetry operations for the lattice structure in the input file.
- Setting up "special k-points" (see above).
- Creating a "crystallographic" file for band structure calculations.

Notes:

- Testing their efficiency and preparation of the "structural" file for running the self-consistent band structure programs.

- In cases where the point group of the crystal does not contain inversion,
  inversion is artificially added in order to make use of the hermiticity of the Hamiltonian.

K213: Preparation of an input file with pseudopotential for running K207
-----------------------------------------------------------------------------------------------

History:

- Written on July 24, 1981, based on `K97`.
- Modified in April 1982 for non-local potentials by Ole Holm Nielsen.

Purpose::

- Preparation of an input file with pseudopotential for running `K95` (or similar pseudopotential programs), including:

  - generation of reciprocal lattice vectors,
  - establishment of a synoptic table of reciprocal space,
  - calculation of atomic pseudopotentials,
  - combination of the total ionic pseudopotential from the individual atomic ones,
  - definition of the initial pseudopotential for the first run of the self-consistent loop.

- These pseudopotentials were implemented:

  - Appelbaum-Hamann (`POT13`)
  - Berkeley ionic (`POT15`)
  - Hamann-Schlüter-Chiang non-local (`POT16`)
  - Bell Labs 1982 "Periodic Table" non-local (`POT17`)

The initial pseudopotential:

- The initial pseudopotential for the self-consistent loop is chosen to be the ionic one, screened by the dielectric function of a free electron gas.

K214: Potential File Processing
-----------------------------------------------------------------------------------------------

Purpose:

- Read a file with ionic and screened potentials and modify it to create an improved starting point for a new calculation with a different lattice structure.

Method: 

- Add the screening potential from a previous self-consistent calculation
- Combine with V(IONIC) 
- Include a screened term of DELTA-V(IONIC)

K207 self-consistent DFT (Density Functional Theory) calculation of Total Energy, Forces, and Stresses
----------------------------------------------------------------------------------------------------------

Purpose:

- Selfconsistent calculation of band structure and total energy, Hellmann-Feynman forces and stress for semiconductors and metals.

Input files

* Input files from *K290* and *K213*.

Output files

* ToDo
