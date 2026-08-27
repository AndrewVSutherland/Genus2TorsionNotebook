# `M(2,2,2,8)` and rational `3`-torsion

The `M(2,2,2,8)` model is

```text
C: y^2 = f(x) = x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
```

with

```text
s2(a,b,c,d)^2 = 4abcd.
```

## Cubic-contact condition

Let

```text
f = x^5 + e1 x^4 + e2 x^3 + e3 x^2 + e4 x,
q = x^2 + Ux + V,
h = m x^3 + N x^2 + R x + S.
```

A nonzero degree-2 divisor class represented by `q` is killed by `3` if and only if there is a cubic `h` such that

```text
h(x)^2 - f(x) = m^2 q(x)^3.
```

Geometrically, the cubic `y=h(x)` has triple contact with `C` at the two points of the divisor.  Equivalently,

```text
div(y-h) = 3D - 6 infinity.
```

This captures rational `3`-torsion: in the odd-degree genus-2 model, a nontrivial odd-torsion class cannot have reduced degree `1`, since that would require a rational function with a single pole of order `3` at infinity.

The coefficient equations are:

```text
2mN - 1 = 3m^2U
N^2 + 2mR - e1 = 3m^2(U^2 + V)
2mS + 2NR - e2 = m^2(U^3 + 6UV)
R^2 + 2NS - e3 = 3m^2(U^2V + V^2)
2RS - e4 = 3m^2UV^2
S^2 = m^2V^3
```

with the nondegeneracy conditions

```text
m != 0,  U^2 - 4V != 0,  gcd(q,f)=1.
```

Thus the rational `3`-torsion condition is an explicit algebraic cover of the `M(2,2,2,8)` K3 surface.  One can eliminate `N,R,S` from the first three equations and keep the last three equations in `m,U,V` over the K3 parameters.

## Eliminated section equations

The script

```text
code/m2228_three_torsion_equations.m
```

pulls these equations back to the known rational section families on `M(2,2,2,8)`.

Since

```text
S^2 = m^2 V^3
```

and `gcd(q,f)=1`, we may write

```text
V = v^2,   S = m v^3,   L = 1/m.
```

Let `M=L^2` and

```text
A = 2M e1 + 6(U^2+v^2) - (M+3U)^2.
```

Then `N,R,S` are recovered from

```text
N = (L^2 + 3U)/(2L),
R = A/(4L),
S = v^3/L,
```

and the remaining three equations are

```text
(M+3U)A + 8v^3 - 4e2M - 4U^3 - 24Uv^2 = 0
A^2 + 16(M+3U)v^3 - 16e3M - 48(U^2v^2+v^4) = 0
A v^3 - 2e4M - 6Uv^4 = 0
```

with

```text
L v (U^2-4v^2) != 0,   gcd(x^2+Ux+v^2, f)=1.
```

The full pulled-back polynomials for the named section families are written to

```text
data/m2228_three_torsion_section_equations.txt
```

A summary run

```text
magma -b mode:="summary" code/m2228_three_torsion_equations.m
```

gave:

```text
k3_section_P:         t-degrees 12,24,12; terms 114,367,71
k3_section_P_plus_T1: t-degrees 10,16,12; terms 50,131,35
k3_section_P_plus_T2: t-degrees 12,24,12; terms 92,291,65
k3_section_P_plus_T3: t-degrees 14,24,14; terms 100,315,64
filip_projective:     t-degrees 14,20,24; terms 38,81,31
```

## Local triple-contact residues

The same script can solve the eliminated equations over finite fields:

```text
magma -b mode:="local" code/m2228_three_torsion_equations.m
```

For the named section families, the output has a strong boundary flavor.  For example:

```text
k3_section_P:
  p=13: good 6, bad 7, good_triple_contact 0
  p=19: good 6, bad 13, good_triple_contact 0
  p=23: good 6, bad 17, good_triple_contact 0
  p=29: good 18, bad 11, good_triple_contact 0
  p=31: good 12, bad 19, good_triple_contact 0

k3_section_P_plus_T1:
  p=13: good 6, bad 7, good_triple_contact 0
  p=19: good 6, bad 13, good_triple_contact 0
  p=23: good 6, bad 17, good_triple_contact 0
  p=31: good 12, bad 19, good_triple_contact 0

filip_projective:
  p=13: good 6, bad 7, good_triple_contact 0
  p=19: good 6, bad 13, good_triple_contact 0
  p=23: good 6, bad 17, good_triple_contact 0
  p=31: good 12, bad 19, good_triple_contact 0
```

Thus on these section families, a rational `3`-torsion point would have to reduce to the bad/boundary section residues at several primes.  This matches the height-search output below: the obstruction is not merely numerical point-count noise; it is visible directly on the triple-contact cover.

## Initial finite-field filter

The script

```text
code/m2228_three_torsion_search.m
```

prints the contact equations and applies the necessary condition

```text
3 divides #J(F_p)
```

for every tested good prime `p`.  This is only a necessary condition, but it is a very cheap first filter before solving the cubic-contact equations.

Named section-family run:

```text
magma -b height:=25 mode:="named" code/m2228_three_torsion_search.m
```

Result:

```text
total 3995
usable 3974
candidates 0
```

Section-multiple run:

```text
magma -b height:=20 max_multiple:=4 mode:="multiples" code/m2228_three_torsion_search.m
```

Result:

```text
total 7665
usable 7618
candidates 0
```

So the first low-height section families and the first four multiples of the elliptic-fibration section have no rational `3`-torsion candidate even at the necessary finite-field level.


## Full K3 surface sieve

The more promising route is to leave the named rational curves and test the full `M(2,2,2,8)` K3 surface.  The script

```text
code/m2228_three_torsion_surface_sieve.m
```

has two modes.  In `finite` mode, it precomputes all coefficient tuples `(e1,e2,e3,e4)` over `F_p` satisfying the eliminated triple-contact equations, then scans all good affine K3 points

```text
s2(a,b,c,d)^2 = 4abcd.
```

A run

```text
magma -b mode:="finite" code/m2228_three_torsion_surface_sieve.m
```

gave:

```text
p=11: good_k3_points 240,   contact_points 240
p=13: good_k3_points 288,   contact_points 0
p=17: good_k3_points 1920,  contact_points 768
p=19: good_k3_points 1728,  contact_points 432
p=23: good_k3_points 5280,  contact_points 528
p=29: good_k3_points 12768, contact_points 2016
p=31: good_k3_points 16560, contact_points 6480
p=37: good_k3_points 31104, contact_points 8640
```

Thus the full K3 surface is not locally empty for the triple-contact cover in general, but `p=13` gives a genuine good-reduction obstruction: every rational example must reduce to the boundary at `13`.

In `exact_file` mode, the script reads existing primitive K3 tuple files, keeps only tuples with bad reduction at `13`, applies good-prime point-count bounds for rational `3`-torsion, and only then calls Magma's exact `TorsionSubgroup`.

Runs on existing tuple files gave:

```text
data/surface_tuples_B2000.txt:
  forced13 75, no3_bound 75, survived_bound 0, exact 0, hits 0

data/surface_tuples_B5000_bad11_23.txt:
  forced13 53, no3_bound 53, survived_bound 0, exact 0, hits 0

data/surface_tuples_B10000_strata_non_d_zero_11_23.txt:
  forced13 94, no3_bound 94, survived_bound 0, exact 0, hits 0

data/surface_tuples_B10000_strata_dcollisions_11_23.txt:
  forced13 90, no3_bound 90, survived_bound 0, exact 0, hits 0
```

I also tried to enumerate all height-`10000` and height-`5000` tuples with only the forced `13`-boundary condition, and then the `13:D` component at heights `10000` and `3000`.  Those runs were too weakly filtered to finish quickly, so they were interrupted and no data rows were kept.

Conclusion: the full-surface route is more promising than the named curves because it has good triple-contact points modulo most primes.  However, the `p=13` obstruction forces a boundary search.  A useful next step is not an undifferentiated `13:any` enumeration, but a sharper normalized `13`-adic boundary analysis or a component-wise enumerator with progress/chunking.


## Normalized `13`-adic boundary analysis

The script

```text
code/m2228_three_torsion_padic13.py
```

analyzes the bad-reduction boundary of the full K3 triple-contact cover at `p=13`.  Unlike the good-reduction finite-field sieve, this allows bad cover reduction: `v` may vanish, `q` may have repeated reduction, and `q` may meet `f` modulo `13`.  The only normalization kept is that `L=1/m` is a unit, which follows from the first contact equation.

A run

```text
python3 code/m2228_three_torsion_padic13.py --output data/m2228_three_torsion_padic13_report.txt
```

gave `3264` mod-`13` boundary solutions of the K3 plus triple-contact equations.  Their base signatures are only zero-boundary signatures:

```text
Z1, Z2, Z3, Z4: 720 each
triple-zero/one-unit charts: 24 each
mixed Zi+Ejk charts: 24 each
pure collision charts Eij: 0
```

All `3264` solutions lift to `13^2`, but the first-order behavior is sharply constrained:

```text
has_p2_lift                  3264
can_make_v_exact1            2592
can_make_zero_coords_exact1    96
can_satisfy_all_first_order    96
```

The `96` first-order live solutions are exactly the four triple-zero/one-unit charts:

```text
Z1 Z2 Z3 with d unit
Z1 Z2 Z4 with c unit
Z1 Z3 Z4 with b unit
Z2 Z3 Z4 with a unit
```

For the simple-zero charts `Zi` and the mixed charts `Zi+Ejk`, the zero coordinate is forced to remain `0 mod 13^2`.  Thus those branches are deeper zero charts, not first-order branches.

Conclusion: the normalized `13`-adic boundary analysis reduces the promising full-surface `3`-torsion search to the triple-zero/one-unit charts.  A plain height enumeration over `13:any` is much too broad; the next finite computation should target these four charts directly, with local coordinates such as

```text
a = 13 A, b = 13 B, c = 13 C, d unit,
A+B+C = 0 mod 13
```

and the three analogous permutations.


## Targeted triple-zero `13`-adic chart search

The normalized `13`-adic boundary analysis says that the only first-order live branches of the full-surface `3`-torsion cover are the four triple-zero/one-unit charts, e.g.

```text
a = 13 A,  b = 13 B,  c = 13 C,  d unit,
A + B + C = 0 mod 13.
```

I implemented a targeted enumerator for these charts:

```text
code/enumerate_m2228_triplezero13.cpp
```

Writing `C = -A-B+13K`, the K3 discriminant condition for solving for the unit coordinate becomes

```text
-A*B*(13K-A)*(A+B)*(13K-B)*(13K-A-B) is a square.
```

The enumerator scans this square condition directly, reconstructs the corresponding signed K3 tuple, canonicalizes the curve tuple by absolute values and sorting, and writes primitive curve tuples for

```text
C: y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
```

Runs:

```text
code/enumerate_m2228_triplezero13 80 
  data/m2228_3torsion_triplezero13_H80.txt 5000 -80 80

checked=293655
exact_square=23080
unique=19
```

```text
code/enumerate_m2228_triplezero13 300 
  data/m2228_3torsion_triplezero13_H300.txt 20000 -300 300

checked=16301269
exact_square=327028
unique=69
```

```text
code/enumerate_m2228_triplezero13 1000 
  data/m2228_3torsion_triplezero13_H1000.txt 50000 -1000 1000

checked=611861896
exact_square=3620535
unique=190
```

```text
code/enumerate_m2228_triplezero13 2000 
  data/m2228_3torsion_triplezero13_H2000.txt 200000 -2000 2000

checked=4908957887
exact_square=14480390
unique=434
```

Each generated file was passed to the existing exact tuple-file sieve:

```text
magma -b mode:="exact_file" 
  tuple_file:="data/m2228_3torsion_triplezero13_H2000.txt" 
  code/m2228_three_torsion_surface_sieve.m
```

The largest run gave:

```text
tuples 434
forced13 391
no3_bound 391
survived_bound 0
exact 0
hits 0
```

So the targeted chart search found no rational `3`-torsion candidate: every forced-`13` tuple was killed by good-prime Jacobian order bounds before exact torsion computation.  This is stronger negative evidence than the earlier broad height files, because it works directly on the `13`-adic branches singled out by the local triple-contact analysis.
