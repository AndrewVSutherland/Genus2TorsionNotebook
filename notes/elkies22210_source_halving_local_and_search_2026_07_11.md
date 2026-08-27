# Elkies `[2,2,2,10]` source-halving covers: local and bounded search

Date: 2026-07-11.

## Target and exact covers

On Elkies' Clebsch--Klein surface

```text
sum_i r_i = 0,
sum_i r_i^3 = 0,
C_r: y^2 = x*product_i(x-r_i^2),
```

the smooth open has torsion containing `[2,2,2,10]`, of order `80`.
If any nonzero rational `2`-class has a rational half, then the torsion
contains

```text
[2,2,2,20],  order 160.
```

There are two `S_5`-orbits of nonzero `2`-classes.  Put `a_i=r_i^2`.
The exact Stoll--Zarhin radicands are derived and certified in

```text
code/elkies22210_halving_covers.m.
```

For orbit `01`, represented by `{0,a_1}`, they reduce to

```text
a_1-a_2, a_1-a_3, a_1-a_4, a_1-a_5.
```

For orbit `12`, represented by `{a_1,a_2}`, they are

```text
-(a_1-a_3)(a_1-a_4)(a_1-a_5),
 (a_3-a_2)(a_1-a_4)(a_1-a_5),
 (a_4-a_2)(a_1-a_3)(a_1-a_5),
 (a_5-a_2)(a_1-a_3)(a_1-a_4).
```

All four radicands must be rational squares.  Merely requiring their
quotients to be squares omits the common squareclass and gives false
positives.  The universal script checks all 15 classes on Elkies' printed
source against Magma's exact `IsDivisibleBy`, and includes both a positive
control and a deliberate common-squareclass negative control.  Each generic
cover has degree `16` over the Clebsch--Klein surface.

## Exact finite-field diagnostic

The script

```text
code/elkies22210_source_halving_finite.m
```

enumerates the projective Clebsch--Klein surface, restricts to

```text
r_i != 0,  r_i^2 pairwise distinct,
```

and tests the four exact radicands.  Negative and positive samples are
independently checked using `AbelianGroup(J)`.

Run through `101`:

```text
magma -b prime_min:=5 prime_max:=101 validate_limit:=2 \
  code/elkies22210_source_halving_finite.m
```

The marked counts for fixed representatives of the two `S_5`-orbits are:

The compact machine-readable summary is
`data/elkies22210_source_halving_finite_p101_summary.txt`.

| `p` | open CK points | marked `01` | marked `12` |
|---:|---:|---:|---:|
| 5 | 0 | 0 | 0 |
| 7 | 0 | 0 | 0 |
| 11 | 24 | 0 | 0 |
| 13 | 0 | 0 | 0 |
| 17 | 0 | 0 | 0 |
| 19 | 120 | 0 | 0 |
| 23 | 240 | 24 | 0 |
| 29 | 480 | 0 | 24 |
| 31 | 504 | 0 | 24 |
| 37 | 720 | 48 | 60 |
| 41 | 1104 | 72 | 120 |
| 43 | 960 | 24 | 72 |
| 47 | 1320 | 72 | 72 |
| 53 | 1800 | 48 | 48 |
| 59 | 2640 | 96 | 132 |
| 61 | 2784 | 144 | 108 |
| 67 | 3240 | 144 | 144 |
| 71 | 3744 | 216 | 168 |
| 73 | 3960 | 168 | 288 |
| 79 | 5040 | 408 | 240 |
| 83 | 5400 | 288 | 384 |
| 89 | 6360 | 264 | 276 |
| 97 | 7560 | 264 | 324 |
| 101 | 8544 | 336 | 444 |

Thus every rational point on either cover is forced to the Clebsch--Klein
boundary at `11` and `19`.  More precisely:

```text
orbit 01: boundary at 11,19,29,31;
orbit 12: boundary at 11,19,23.
```

Here boundary means that some `r_i` vanishes or two `r_i^2` collide after
projective reduction.  This is not yet a `Q_p` impossibility proof: bad or
nonintegral boundary disks require their own blowups.  The zero open counts
at `5,7,13,17` come from the Clebsch--Klein open itself, not from a
cover-specific distinction.

## Complete integral-projective height search

The exact C++ search is

```text
code/elkies22210_source_halving_search.cpp.
```

It uses the quadratic Clebsch--Klein enumerator from the earlier Richelot
sweep, removes permutation and global-sign symmetry by ordering the five
absolute values, and tests every marked class using the exact four
Stoll--Zarhin radicands.  It is complete for primitive integral projective
tuples with `max |r_i| <= H`.

At `H=100` it reproduces exactly the `94` unique source curves from the old
Magma sweep.  The completed run at `H=2000` gives:

```text
triples                         5,325,336,000
square discriminants                 519,165
integral quadratic roots             285,794
open CK representatives               23,855
unique primitive source curves        10,290
orbit 01 exact points                       0
orbit 12 exact points                       0
```

Among these `10,290` curves, the necessary forced-boundary filters retain

```text
11,19 boundary:              6,027
orbit 01 at 11,19,29,31:    1,214
orbit 12 at 11,19,23:       3,109.
```

None of these boundary-compatible curves satisfies the exact rational
radicands.

## Orbit-12 boundary lifts at `p=11,19,23`

The normalized classifier

```text
code/elkies22210_orbit12_boundary_lifts.py
```

fixes the first `p`-adic unit among `r_1,...,r_5` to `1`, so its five
charts include the cases in which either or both coordinates in the marked
pair `{r_1^2,r_2^2}` are nonunits.  It solves the two Clebsch--Klein
equations exactly modulo `p^k` and applies the four exact orbit-12 square
conditions.  A zero radicand modulo `p^k` is retained but labelled `deep`,
since such a truncated zero need not have even finite valuation.

Run with

```text
python3 code/elkies22210_orbit12_boundary_lifts.py --p 11 --max-k 3
```

The exact counts are:

| level | modulus | cover states | resolved | deep |
|---:|---:|---:|---:|---:|
| 1 | `11` | 129 | -- | -- |
| 2 | `121` | 3,699 | 0 | 3,699 |
| 3 | `1,331` | 182,289 | 39,960 | 142,329 |

All 129 normalized level-1 states have Clebsch--Klein Jacobian rank `2`;
before imposing the cover condition each therefore has `11^2=121` tangent
lifts.  After the cover filter, the numbers of level-2 children per
level-1 state are

```text
1 child:     3 states
11 children: 105 states
121 children: 21 states.
```

All 38 level-1 boundary strata survive through level 3, but only 19 have
a descendant whose four radicand valuations have resolved.  The marked
unit/nonunit distribution among the 39,960 resolved descendants is

```text
UU 21,780,  UN 7,260,  NU 7,290,  NN 3,630.
```

In particular, the computation certifies a genuine open `Q_11` branch.
One normalized representative modulo `11^3=1331` is

```text
r = (1,242,959,9,120).
```

It satisfies both Clebsch--Klein equations modulo `1331`.  Its four
orbit-12 radicands, displayed with square roots modulo `1331`, are

```text
(484,1089,1089,1197) = (22^2,33^2,33^2,322^2) mod 1331,
valuation signature = (2,2,2,0).
```

The fixed-chart Clebsch--Klein Jacobian has rank `2` modulo `11`, so
multivariate Hensel lifting gives a `Z_11` solution through this residue
class.  The finite even radicand valuations and square leading units remain
unchanged in that class, hence all four radicands are squares in `Q_11`.
Moreover every `r_i` and every difference `r_i^2-r_j^2` is already nonzero
modulo `1331` (the ten differences are
`1,41,1251,242,40,1250,241,1210,201,322`), so the lifted genus-2 curve is
on the smooth open.  Thus `p=11` forces bad reduction but is **not** a local
obstruction to the orbit-12 cover.

The same classifier resolves the other two forced orbit-12 primes already
at level 2:

| `p` | states mod `p` | states mod `p^2` | resolved mod `p^2` | deep mod `p^2` |
|---:|---:|---:|---:|---:|
| 19 | 234 | 13,410 | 2,166 | 11,244 |
| 23 | 252 | 19,128 | 3,174 | 15,954 |

For `p=19`, an explicit normalized seed modulo `19^2=361` is

```text
r = (1,248,98,3,11),
G = (23,137,137,175) = (97^2,55^2,55^2,135^2) mod 361.
```

Its four radicands are units, the fixed-chart CK Jacobian has rank `2`
modulo `19`, and its ten square differences modulo `361` are

```text
228,144,353,241,277,125,13,209,97,249.
```

For `p=23`, an explicit normalized seed modulo `23^2=529` is

```text
r = (1,392,118,8,10),
G = (285,331,400,262) = (49^2,233^2,20^2,43^2) mod 529.
```

Again all four radicands are units, the CK Jacobian has rank `2` modulo
`23`, and the ten square differences modulo `529` are

```text
276,360,466,430,84,190,154,106,70,493.
```

In both cases the unit square conditions and the smooth-open inequalities
are stable throughout the Hensel disk.  Consequently **all three forced
orbit-12 primes `11,19,23` are locally soluble**.  Further attempts to turn
one of these boundary primes into a local obstruction are therefore not the
right next step.

## Compatible-prime/CRT global sieve

The complete integral-projective search now applies a rigorous compatible
local sieve to every fixed marked pair.  After clearing projective
denominators, each exact orbit-12 radicand is an integer square; hence it
must be a square residue modulo

```text
11^3 = 1331,  19^2 = 361,  23^2 = 529.
```

The degree-six radicands change by a sixth power under projective scaling,
so the test is invariant under local unit normalization.  Truncated zero
radicands are retained as square residues, making this a necessary sieve
even in unresolved deep boundary disks.  Testing the same marked pair at
the three primes is equivalent to a CRT condition modulo

```text
1331 * 361 * 529 = 254,179,739.
```

This is implemented in

```text
code/elkies22210_source_halving_search.cpp.
```

The complete rerun at projective height `2000` gives:

| sieve stage | marked pairs | source curves with a surviving pair |
|---|---:|---:|
| all primitive CK curves | 102,900 | 10,290 |
| square residues mod `11^3` | 2,170 | 986 |
| also mod `19^2` | 51 | 36 |
| also mod `23^2` | **0** | **0** |

Thus no marked pair in the entire height-`2000` box even survives the three
compatible local square tests.  This strengthens the previous zero exact-hit
result and explains why enlarging the same blind box has poor expected value.
It is a bounded global congruence result, not a Hasse obstruction: the three
explicit Hensel disks above prove local solubility prime by prime.

For a concrete next-stage lattice target, combine the three displayed seeds
coordinate by coordinate, keeping the marked pair `{r_1^2,r_2^2}` fixed.
The script

```text
code/elkies22210_orbit12_crt_seed.py
```

returns the centered class

```text
(1,-1964314,-56831410,71155269,-12359546) mod 254179739.
```

This representative happens to satisfy `sum r_i=0` exactly and satisfies
`sum r_i^3=0` modulo the full CRT modulus.  Its cubic sum is not zero, so it
is deliberately recorded only as a seed for an exact lattice lift, not as a
curve or cover point.

## Direct orbit-01 chart

For orbit `01`, scale `r_1=1` and parametrize the four unit circles by

```text
r_j = (1-t_j^2)/(1+t_j^2),
z_j = 2*t_j/(1+t_j^2).
```

Then `1-r_j^2=z_j^2` is automatic.  The two Clebsch--Klein equations reduce
exactly to

```text
sum_j 1/(1+t_j^2) = 3/2,
product_j (1-t_j^2) = 16.
```

The meet-in-the-middle script

```text
code/elkies22210_orbit01_unit_circle_search.py
```

hashes the sum and product invariants of pairs.  At rational parameter
height `50` it tested

```text
1,546 t-values,
1,193,512 admissible pairs,
0 complementary pair joins,
0 points.
```

### Boundary lifting at `3` and `7`

The affine `t`-equations have no solutions modulo `3` or modulo `7`, but
this does **not** give a local obstruction: a rational point may reduce to a
chart at infinity.  The proper `(P1)^4` audit is implemented in

```text
code/elkies22210_orbit01_boundary_lifts.py.
```

It homogenizes `E4-E1-5` and `E3-E2+10`, enumerates all finite/infinity
charts, and lifts their exact tangent systems.  The state counts are

```text
p=3:  56 mod 3  -> 1,512 mod 9  ->    13,608 mod 27,
p=7: 248 mod 7  ->28,616 mod 49 -> 1,402,184 mod 343.
```

More importantly, quantitative Hensel witnesses prove that genuine
smooth-open branches survive, rather than only components of the projective
boundary.  For `p=3`, infinity mask `1` and chart variables

```text
(9,7,10,22)
```

give equation values `(2319921,185938497)`.  On Jacobian columns `(0,1)`,
the determinant has `3`-adic valuation `2`, while both equations have
valuation at least `5`.  The Newton correction therefore has valuation at
least `3`, strictly larger than the maximum valuation `2` of every
smooth-open factor.  The limiting root remains off the boundary.

For `p=7`, mask `1` and variables

```text
(7,246,212,3)
```

give equation values `(24473311709,1068847960173)`.  Columns `(1,2)` have
unit determinant; the equation valuations and guaranteed correction
valuation are at least `3`, again strictly above the maximum open-factor
valuation `2`.  Hence the orbit-01 cover has smooth-open points over both
`Q_3` and `Q_7`; neither prime is a local obstruction.

## Recommendation

The source-halving idea remains mathematically natural, but it is not an
open-density search: every rational point is forced into several simultaneous
bad-reduction disks.  All forced orbit-12 primes are locally soluble, so the
local-obstruction phase is complete.  The useful pivot is now global: use a
rational parametrization of the Clebsch cubic surface and impose the three
compatible boundary signatures on its two parameters, or explicitly CRT-lift
selected open Hensel disks and solve the exact CK equations in those lattice
cosets.  The new congruence sieve should remain the inexpensive front end for
either search.  A larger blind Clebsch--Klein height box is lower value.
