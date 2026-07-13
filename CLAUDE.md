# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A historical DFT plane-wave code (Total Energy, Forces, Stresses with pseudopotentials), written 1979-1983 at Xerox PARC by Karel Kunc, Richard M. Martin, and Ole Holm Nielsen. It is Fortran-77 being resurrected on modern Linux/gfortran. Treat the physics and the FORTRAN idioms (fixed-form, 6-space indent, uppercase, COMMON blocks, no dynamic allocation) as intentional — match the surrounding style rather than modernizing.

## Build and run

```sh
make run290      # K290: crystal structure -> symmetry ops + special k-points
make run213      # K213: reciprocal lattice + pseudopotential -> starting potential
make run207      # K207: self-consistent DFT (real Hamiltonian)
make crun207     # K207: complex-Hamiltonian variant, built via patch/*.diff
make clean       # removes *.o, generated *.for, executables
./testscript.sh  # end-to-end smoke test: builds and runs the pipeline on 1-cell BN
```

`make` alone does nothing (default target `none:` is empty). There are no unit tests; `testscript.sh` is the test — it writes the structure file, exports dimensions, builds each program, and feeds it a heredoc of interactive answers. It runs K290 and K213 to completion on 1-cell BN (producing `fort.3` and the potential `fort.4`), then `exit 0`s, since the K207 stage is not yet runnable.

`FC`/`FFLAGS` are set at the top of the Makefile (`gfortran -cpp -Dx86_64 -g -std=legacy`). The `-cpp` matters: sources contain `#ifdef` blocks for the historical machines (sun, apollo, hppa, x86_64).

## Repository state

Only part of the original source has been committed. `k207.for` (the driver) exists, but most of its subroutines listed in `OBJ_207` (`cfft`, `diagon`, `forces`, `excorr`, `mixvg`, `lookup`, `apwsum`, `sciblas`, `unixtra`, …), the `run207.start` main program, all the `c*.diff` patch files for the complex build, `k214.for`, and the `*.VG` pseudopotential data files are **absent**. So `make run290` and `make run213` work today; `run207`, `crun207`, `run214`, `stressfield`, `testeis`, `testion`, and the `flint` targets do not. When a target fails on a missing `.for`, check whether the file was ever committed before assuming a build bug.

## Pseudopotentials

`PARM17` (`paramete.for`) is the parameter database for POT17, the Bachelet-Hamann-Schlüter "periodic table" non-local pseudopotential (Phys. Rev. B **26**, 4199 (1982); the paper is at `../PhysRevB.26.4199.pdf`). It covers ten elements — `IATOMS = 1,5,6,7,13,14,15,31,32,33` (H, B, C, N, Al, Si, P, Ga, Ge, As) — and all ten have been verified value-by-value against Table IV of the paper, including the end-to-end check that the code's silicon V(r) reproduces the paper's Table V to <0.002 a.u. Do not "fix" these numbers.

Two things that look like omissions but are faithful to the paper: H/B/C/N/Al/Si carry `LSPIN = -1` because Table IV has no spin-orbit row for them (only P, Ga, Ge, As do); and the A-coefficients `PARM17` derives from the tabulated cᵢ are huge (~10⁵) with heavy cancellation, which is inherent to the BHS Cholesky transform, not a bug. An element outside `IATOMS` returns `IFOUND = 0` and K213 simply re-prompts — adding one means extending `IATOMS`/`NUMATM` and the computed `GOTO`.

Because B and N are in the POT17 database, `testscript.sh` selects potential 17 and needs no external data files. The alternative, potential 20 (numerically generated), reads Fourier-transformed potentials from `B.VG`/`N.VG`, which are **not** in the repo.

## How the three programs fit together

Sequential pipeline, communicating through numbered Fortran units (files `fort.N`), not through APIs:

1. **K290** reads the crystal structure (`fort.2`), finds the point group, and generates special k-points (Monkhorst-Pack, MacDonald-modified, plus OHN "symmetrization"). Writes the structural file `fort.3`. Inversion is artificially added when the point group lacks it, so the Hamiltonian stays hermitian.
2. **K213** reads `fort.3` plus atomic pseudopotential files (unit `20+i` per atom type, e.g. `B.VG`, `N.VG`), builds reciprocal-lattice shells and the total ionic pseudopotential, and screens it with a free-electron dielectric function to get the starting potential. Writes `fort.10`.
3. **K207** self-consistent loop: band structure, total energy, Hellmann-Feynman forces, stress. Reads the previous stage's potential as `fort.4` — hence `mv fort.10 fort.4` between stages in `testscript.sh`.

Unit numbers are assigned once in `consts.for` (`COMMON /FILES/`) and documented in its header: 2 = structure, 3 = K290 output, 4 = potential input, 5/6 = terminal in/out, 10 = potential output, 11 = nonlocal matrix elements, 12 = summary, 13/14 = k-point scratch. Every program calls `DAY` then `CONSTS` first to initialize `/MACH1/`, `/CONST/`, and `/FILES/`.

All three programs are driven interactively from stdin; runs are scripted by piping a heredoc of answers (see `testscript.sh` for the exact question order).

## Array dimensions are a build-time code generation step

There is no dynamic allocation. The main programs are templates named `*.start` (`run290.start`, `run213.start`) containing `$NDIM1`-style placeholders inside `PARAMETER` statements. The Makefile's `.start.for` rule sources `dimensions.sh` and `sed`-substitutes those variables to produce `run290.for` etc., which are then compiled. Generated `.for` files are regenerated every build (forced via `FRC`) and deleted by `make clean` — never edit or commit them; edit the `.start` file.

To change problem size (number of atoms, plane waves, FFT grid), edit `dimensions.sh` — but note `testscript.sh` duplicates the same variable block inline and that copy wins for the test run. Keep the two in sync.

Gotcha: the sed rules run in order and `s/$NDIM1/.../` precedes `s/$NDIM10/.../`, so `$NDIM10` (used by `run290.start`, and *not* defined in `dimensions.sh`) is silently expanded as `$NDIM1` followed by a literal `0` — with `NDIM1=7000` the special-point limit becomes 70000. This currently compiles and runs, but is not what the placeholder names imply.

## The complex-Hamiltonian variant

`crun207` is not a separate source tree: the `c*.for` files (`cdiagon`, `cdiaham`, `ceval`, `ckinetic`, `clowper`, `clumat`, `cmixvg`, `cputnl`, `cro5`, `crowham`, `crun207`, `crwev`, `cvnlkka`, `cvnlkknum`) are produced at build time by applying `c*.diff` with `patch` to their real-valued counterparts. To change behavior shared by both, edit the real source; to change only the complex path, update the `.diff`. `NCMPLX` in `dimensions.sh` selects which binary `testscript.sh` intends to run (1 = real `run207`, 2 = complex `crun207`).

## Machine portability

`machine.for` hardcodes `MACTYP` (12 = Linux/x86) and `usage.for`/`unixus.c` gate CPU-time and page-fault reporting per machine type with `#ifdef`. The old branches (IBM, CRAY, VAX, CDC, Apollo, Fujitsu) are dead but retained deliberately as history — leave them in place.
