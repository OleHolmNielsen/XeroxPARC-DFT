Xerox PARC (Palo Alto Research Center) DFT plane-wave code
============================================================

This historical project preserves the [DFT](https://en.wikipedia.org/wiki/Density_functional_theory)
[plane-wave](https://en.wikipedia.org/wiki/Plane_wave) 
[pseudo-potential](https://en.wikipedia.org/wiki/Density_functional_theory#Pseudo-potentials)
code written originally
by [Karel Kunc](http://www-ext.impmc.upmc.fr/~kunc/) (<krl.kunc@gmail.com>)
and [Richard M. Martin](https://people.physics.illinois.edu/rmartin/) (<rmartin@illinois.edu>)
at the Xerox [PARC](https://en.wikipedia.org/wiki/PARC_(company)) (*Palo Alto Research Center*) during 1979-1981.

The present version of the original code was further developed by Ole Holm Nielsen at Xerox PARC during 1982-1983.

Author of this page: [Ole Holm Nielsen](https://dcwww.fysik.dtu.dk/~ohnielse/) (<Ole.H.Nielsen@fysik.dtu.dk>).

Project homepage: https://github.com/OleHolmNielsen/XeroxPARC-DFT/

Publications
------------------

The following publications document the theory behind the codes:

* *Theory of structural properties of covalent semiconductors*,
  H. Wendel and Richard M. Martin,
  [Phys. Rev. B 19, 5251, 1979](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.19.5251).

* *Density-functional calculation of static and dynamic properties of GaAs*,
  K. Kunc and Richard M. Martin,
  [Phys. Rev. B 24, 2311(R), 1981](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.24.2311).

* *Ab Initio Force Constants of GaAs: A New Approach to Calculation of Phonons and Dielectric Properties*,
  K. Kunc and Richard M. Martin,
  [Phys. Rev. Lett. 48, 406, 1982](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.406).

* *Theory of static structural properties, crystal stability, and phase transformations: Application to Si and Ge*,
  M. T. Yin and Marvin L. Cohen,
  [Phys. Rev. B 26, 5668, 1982](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.5668).

* *Stresses in semiconductors: Ab initio calculations on Si, Ge, and GaAs*,
  O. H. Nielsen and Richard M. Martin
  [Phys. Rev. B 32, 3792, 1985](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.32.3792).

About the code
-----------------

The code in this historical project implements Total Energy, Forces, and Stresses
with a pseudopotential plane-wave basis set.
It must be emphasized that this code was used for cutting-edge research during the 1980ies,
but is now more like a [time capsule](https://en.wikipedia.org/wiki/Time_capsule)
which preserves the code for posterity.

The code is written in [Fortran-77](https://en.wikipedia.org/wiki/Fortran),
which was the most modern standard Fortran compiler at the time of writing.
Some quirks in the code are due to limitations of [Fortran-77](https://en.wikipedia.org/wiki/Fortran):

* Arrays had to be allocated with static dimensions in the [main program entry-point](https://en.wikipedia.org/wiki/Entry_point)
  (in the ``run2xx.for`` files) since dynamically allocatable arrays were not yet in the Fortran standard.
  Therefore the scripts must define static dimension variables such as ``NDIM1, NDIM2, NDIM3``, etc.,
  to generate ``run2xx.for`` files from the ``run2xx.start`` files.

* The code works both with either Real and Complex matrices (the former being much faster than the latter)
  depending on the symmetry properties of the crystal.
  Therefore a number of ``xxx.diff`` files are used to generate the Complex code version 
  from the Real version.
  The Fortran compiler will likely issue warnings about ``passed COMPLEX(4) to REAL(4)``
  (or similar) due to the dirty programming tricks employed in the code.
  Remember that nice code [Preprocessors](https://en.wikipedia.org/wiki/Preprocessor)
  did not exist at the time.

Running the set of codes
=============================

A run of the code consists of several sequential code steps described in the sections below.
An example of a script running the codes is in the file [testscript.sh](testscript.sh).

Firstly, create a structure file (Fortran unit 2) describing the
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
- Generation of special points modified on May 26, 1982 by OHN.

Purpose:

- Setting up special k-points and creating a "crystallographic" file for band structure calculations.

Generation of special points in k-space for an arbitrary lattice, following the method of:

- [Monkhorst and Pack, *Phys. Rev. B* **13** (1976) 5188](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.13.5188)
- Modified by [MacDonald, *Phys. Rev. B* **18** (1978) 5897](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.18.5897)
- Modified also by Ole Holm Nielsen ("symmetrization")

Notes:

- Testing their efficiency and preparation of the "structural" file for running the self-consistent band structure programs.

- In cases where the point group of the crystal does not contain inversion, inversion is artificially added in order to make use of the hermiticity of the Hamiltonian.

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

Initial pseudopotential:

- The initial pseudopotential for the self-consistent loop is chosen to be the ionic one, screened by the dielectric function of a free electron gas.

Potentials currently available:

- Appelbaum-Hamann (`POT13`)
- Berkeley ionic (`POT15`)
- Hamann-Schlüter-Chiang non-local (`POT16`)
- Bell Labs 1982 "Periodic Table" non-local (`POT17`)

Bell Labs 1982 "Periodic Table" non-local pseudopotentials
.............................................................

The subroutine `PARM17` (file `paramete.for`) hardcodes pseudopotential parameters for ten elements.

In `bhs_table4.md` and the CSV-file `bhs_table4.csv`
we show a Machine-readable transcription of **Table IV** of the publication:

> G. B. Bachelet, D. R. Hamann and M. Schlüter,
> *"Pseudopotentials that work: From H to Pu"*,
> Phys. Rev. B **26**, 4199 (1982). doi:10.1103/PhysRevB.26.4199

with the correction from the erratum

> Phys. Rev. B **29**, 2309 (1984). doi:10.1103/PhysRevB.29.2309

The article is © 1982 American Physical Society; the tabulated parameters are
numerical data reproduced here for use with the POT17 pseudopotential in this code.

This file covers all 94 (H through Pu).

K214: Potential File Processing
-----------------------------------------------------------------------------------------------

Purpose:

- Read a file with ionic and screened potentials and modify it to create an improved starting point for a new calculation with a different lattice structure.

Method: 

- Add the screening potential from a previous self-consistent calculation
- Combine with V(IONIC) 
- Include a screened term of DELTA-V(IONIC)

K207 self-consistent DFT (Density Functional Theory) calculation of Total Energy, Forces, and Stresses
-----------------------------------------------------------------------------------------------

Purpose:

- Selfconsistent calculation of band structure and total energy, Hellmann-Feynman forces and stress for semiconductors and metals.

ToDo: Further explanations

Input files
...................

* TBD

Output files
...................

* TBD
