Xerox PARC (Palo Alto Research Center) DFT plane-wave code
============================================================

This historical project preserves the DFT plane-wave code written originally
by Karel Kunc and Richard M. Martin at the Xerox PARC (Palo Alto Research Center) during 1979-1981.

The code was further developed by Ole Holm Nielsen at Xerox PARC during 1982-1983.

Author of this page: Ole.H.Nielsen@fysik.dtu.dk

Homepage: https://github.com/OleHolmNielsen/XeroxPARC-DFT/

The code implements Total Energy, Forces, and Stresses
with a pseudopotential plane-wave basis set.

The code is written in Fortran-77, the most modern Fortran compiler at the time.

A run of the code consists of three sequential code steps in the sections below.

K290: Analyses the crystal structure and generates Special k-points and symmetry operations
=================================================================================================

History:

- Written on September 12, 1979.
- IBM-retouched on October 27, 1980.
- Generation of special points modified on May 26, 1982 by OHN.

Purpose:

- Playing with special points and creation of a "crystallographic" file for band structure calculations.

Generation of special points in k-space for an arbitrary lattice, following the method of:

- [Monkhorst and Pack, *Phys. Rev. B* **13** (1976) 5188](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.13.5188)
- Modified by [MacDonald, *Phys. Rev. B* **18** (1978) 5897](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.18.5897)
- Modified also by Ole Holm Nielsen ("symmetrization")

Notes:

- Testing their efficiency and preparation of the "structural" file for running the self-consistent band structure programs.

- In cases where the point group of the crystal does not contain inversion, inversion is artificially added in order to make use of the hermiticity of the Hamiltonian.

K213: Preparation of an input file with pseudopotential for running K207
============================================================================

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

K214: Potential File Processing
============================================================================

Purpose:

- Read a file with ionic and screened potentials and modify it to create an improved starting point for a new calculation with a different lattice structure.

Method: 

- Add the screening potential from a previous self-consistent calculation
- Combine with V(IONIC) 
- Include a screened term of DELTA-V(IONIC)

K207 self-consistent DFT (Density Functional Theory) calculation of Total Energy, Forces, and Stresses
============================================================================================================

Purpose:

- Selfconsistent calculation of band structure and total energy, Hellmann-Feynman forces and stress for semiconductors and metals.

Input files
===============

* TBD

Output files
===============

* TBD
