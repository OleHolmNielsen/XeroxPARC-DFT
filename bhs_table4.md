# BHS pseudopotential parameters (`bhs_table4.csv`)

Machine-readable transcription of **Table IV** of

> G. B. Bachelet, D. R. Hamann and M. Schlüter,
> *"Pseudopotentials that work: From H to Pu"*,
> Phys. Rev. B **26**, 4199 (1982). doi:10.1103/PhysRevB.26.4199

with the correction from the erratum

> Phys. Rev. B **29**, 2309 (1984). doi:10.1103/PhysRevB.29.2309

The article is © 1982 American Physical Society; the tabulated parameters are
numerical data reproduced here for use with the POT17 pseudopotential in this code.

These are the same parameters that `PARM17` (`paramete.for`) hardcodes for ten
elements. This file covers all 94 (H through Pu).

## Format

One row per channel; 586 rows, 94 elements.

| column | meaning |
|---|---|
| `element`, `Z`, `Zv` | symbol, atomic number, valence charge |
| `channel` | `core`, `l0`, `l1`, `l2`, `l3`, or `so1`/`so2`/`so3` (spin-orbit, by L index) |
| `alpha1..alpha3` | Gaussian exponents (`alpha3` empty for `core`) |
| `c1..c6` | coefficients (`core` rows use only `c1`, `c2`) |

Row counts: 94 `core`/`l0`/`l1`/`l2`, 40 `l3`, 79 `so1`, 66 `so2`, 25 `so3`.

Note that the `c` coefficients are *not* the Gaussian amplitudes used directly in
the potential: they are the Cholesky-normalised form. `PARM17` transforms them
into the A-coefficients (which are large, ~10^5, with heavy cancellation — this is
inherent to the transform, not an error).

## Erratum

The erratum corrects one row: **Cs (Z=55), l=3**. The published alphas
(0.36, 0.41, 0.48) and their six coefficients did not correspond to the potential
curve in Fig. 4 of the paper. This file contains the corrected values
(alphas 0.17, 0.31, 0.51). No other element is affected.

## Verification

The scanned PDF's text layer is corrupt (letters substituted for zeros, commas for
decimal points, dropped minus signs), so every value was read from the page images,
then checked:

- **Core sum rule** `c1 + c2 = 1` holds for all 94 elements.
- **Cross-check against `PARM17`**: for the ten elements the Fortran already
  contains (H, B, C, N, Al, Si, P, Ga, Ge, As), 346 values compared, 0 mismatches.
  Those Fortran values in turn reproduce the paper's Table V silicon potential
  (real-space mesh) to better than 0.002 a.u.
- **Structure**: Z = 1..94 with no gaps; six coefficients on every non-core row.

The 84 elements *not* in `PARM17` rest on the visual transcription plus the sum-rule
and structural checks. Those fully validate each `core` row, but would not catch a
transposed digit in, say, `c4` of an `l`-row. Treat them as good but not
independently confirmed.

## Two values to be aware of

Transcribed **as printed**, not silently corrected:

- `Pt` `so1` alphas `(1.26, 0.40, 0.49)` — breaks the ascending-alpha pattern; the
  neighbouring Au and Ir so1 rows start at ~0.26, so this looks like a typo in the
  original (probably 0.26), but it is what the page shows.
- `Po` `l3` alphas `(3.28, 2.58, 2.99)` — also non-ascending, but the paper does
  print genuinely non-monotonic alphas elsewhere (e.g. Ti `l2` = 4.47, 2.03, 14.24),
  so this is probably correct as printed.

Table VI of the paper (alternative Zn, Cd, Hg potentials constructed without the
3d/4d shells in the core) is **not** included; this file has the Table IV versions
of those three elements.
