# The (p=31) boundary of the q-square ([2,2,2,24]) cover

Date: 2026-07-18

## Outcome

This continuation did not find a rational ([2,2,2,24]) Jacobian, but it
resolved a substantial part of the previously undifferentiated (p=31)
boundary and pushed the boundary-conditioned rational search far beyond the
old height-187 box.

* All 196 normalized finite boundary bases modulo 31 were lifted exhaustively
  to (31^2).  Every internal class still lying on its collision divisor was
  then lifted exhaustively to (31^3).
* On the one-node stable fibres, 144 marked presentations lie in a chart in
  which the complete generalized Jacobian test can be performed.  Their
  elliptic-normalization projections all have exact order 12 and none is
  divisible by 2.  Hence all 144 are rigorously eliminated, before even using
  the torus squareclass.
* Complete smooth finite-Jacobian target masks were added at
  (43,47,53,59,61,67,71).
* Two boundary-conditioned Monte Carlo/CRT runs made 700 million trials with
  coordinate heights up to 30,000.  They produced no nontrivial exact
  double-square point.  This is a large targeted experiment, not an
  exhaustive height-30,000 theorem.
* An independent elliptic-fibration search examined 374 external-boundary
  fibres and 1,447 internal-search fibres.  It also found no nontrivial
  double-square point on the required (p=31) boundary.

The unresolved locus is now sharply identified: marked divisors meeting the
node, degenerate CRT representatives of the same phenomenon, the four
(R=S=0) sheet points, and the toric zero/pole chambers.

## 1. Boundary identities

Put

\[
 x=A^2,\qquad y=B^2,\qquad z=(xy)^{-1},\qquad
 R=x+y+z-3,\qquad S=x^{-1}+y^{-1}+z^{-1}-3,
\]

and (w=S/R=(d/s)^2) on the open sheet.  Direct calculation gives

\[
 S-xR=-\frac{(x-1)^3}{x},
\]

and its two permutations.  Consequently the three external collision
divisors (x=w,y=w,z=w) reduce respectively to
(x=1,y=1,z=1), with *cubic* thickness at 31.

On an internal collision, for example (x=y),

\[
 R=\frac{(x-1)^2(2x+1)}{x^2},\qquad
 S=\frac{(x-1)^2(x+2)}{x}.
\]

Thus the residual double-square test is that (2x+1) and (x(x+2))
be squares.  These identities explain the very different lift depths in the
next section.

## 2. Exhaustive lifts to (31^2) and (31^3)

The finite boundary consists of 196 labelled ((A,B))-bases:

| stratum | bases mod 31 |
|---|---:|
| (R=S=0) sheet intersection | 4 |
| each internal collision (x=y,x=z,y=z) | 8 |
| each external collision (x=w,y=w,z=w) | 56 |

Every base has (31^2=961) normalized lifts in ((A,B)).  All 188,356
lifts remain compatible with the two modular square equations.  Their depth
distribution is:

| stratum | compatible mod (31^2) | exact depth 1 | depth at least 2 |
|---|---:|---:|---:|
| each internal stratum | 7,688 | 7,440 | 248 |
| each external stratum | 53,816 | 0 | 53,816 |
| four sheet points together | 3,844 | -- | 3,844 with (v(R),v(S)\ge2) |

The 744 deep internal (31^2)-classes have 714,984 children modulo (31^3).
For each of the three internal strata the counts are

```text
depth exactly 2     230,640
depth at least 3      7,688
```

For the external strata, the cubic identity proves depth at least 3 for
every lift.  It was also checked computationally on 161,448 canonical
(31^3)-children, with zero failures.

This is an exhaustive base-cover/depth calculation.  It is not a claim that
the singular fibre itself is an ordinary finite Jacobian.

## 3. Complete test on open one-node stable charts

Every one of the 192 non-sheet collision bases has exactly one node.  If

\[
 f(X)=(X-r)^2g(X),
\]

then the normalization is the elliptic curve (E:Y^2=g(X)), and the identity
component of the generalized Jacobian is an extension

\[
 1\longrightarrow T\longrightarrow \operatorname{Pic}^0(C_0)
 \longrightarrow E\longrightarrow 0,

\]

where (T) is the split or nonsplit one-dimensional torus.

For a marked divisor disjoint from the node, the program does all of the
following.

1. Divide its (y)-coordinate by (X-r) and map it to (E(\mathbf F_{31})).
2. Enumerate every ordinary elliptic half.
3. Evaluate the tangent-line cocycle at the two points over the node.
4. Test its squareclass in (\mathbf F_{31}^{*}) in the split case, or in
   the norm-one group of order 32 in the nonsplit case.

This is the complete divisibility criterion in the generalized Jacobian,
not only a test of the elliptic projection.  The result is

```text
one-node marked presentations             2,304
complete generalized-Jacobian tests         144
normalization order distribution       12: 144
ordinary elliptic halves                       0
generalized-Jacobian divisible                  0
marked representative hits node               900
CRT determinant/line chart degenerates       1,260
```

All 144 complete cases are on the three internal collision strata, and all
their nodes are nonsplit.  The direct (D_3+D_4) normalization calculation
independently reproduces the same 144 order-12, nondivisible projections.
The remaining 2,160 presentations have (D_3), (D_4), or their reduced sum
meeting the node.  Their answer genuinely depends on the first-order blow-up
chart and is not inferred from a singular `Jacobian` computation.

## 4. Additional complete good-prime masks

The smooth affine q-square cover was exhaustively enumerated at seven more
primes.  Both square-root sheets and all three marked pairings were included;
divisibility was tested in an explicit finite abelian group.

| (p) | smooth sheet points | marked presentations | target ((A,B)) pairs | target presentations |
|---:|---:|---:|---:|---:|
| 43 | 768 | 2,304 | 48 | 288 |
| 47 | 1,152 | 3,456 | 108 | 648 |
| 53 | 1,728 | 5,184 | 36 | 216 |
| 59 | 2,112 | 6,336 | 72 | 432 |
| 61 | 2,688 | 8,064 | 48 | 288 |
| 67 | 3,072 | 9,216 | 60 | 360 |
| 71 | 3,456 | 10,368 | 108 | 648 |

As at the earlier primes, every smooth marked presentation has exact order
12.  The masks retain bad-reduction boundaries rather than treating them as
smooth negative evidence.

## 5. Boundary-conditioned rational/CRT search

The compiled search does not draw a blind pair ((A,B)).  In its normalized
lane it first chooses one of the complete 196 (p=31) boundary residue
classes.  A separate toric lane forces a numerator or denominator divisible
by 31.  It then applies complete target-or-boundary masks at the other primes
before the exact integer-square test.

The two principal runs were:

| trials | coordinate height | generated valid pairs | 7-prime survivors | exact double-square | exact nontrivial |
|---:|---:|---:|---:|---:|---:|
| 200,000,000 | 10,000 | 187,364,686 | 51,843 | 54 | 0 |
| 500,000,000 | 30,000 | 468,350,289 | 130,408 | 11 | 0 |

Every exact square point encountered was on one of the six exact collision
components and therefore defined a singular curve.  A final ten-prime
five-million-trial validation left 11 modular survivors and again no exact
point.  Candidate files contain only their headers.

These runs sample rational coordinates of the stated heights; they do not
enumerate every such fraction pair.  Their rigorous content is the exact
verification of every sampled survivor and the completeness of each modular
mask on its smooth finite locus.

## 6. Elliptic-fibration search

For fixed (A), the first square equation is the genus-one quartic

\[
 W^2=B^4+(A^2-3)B^2+A^{-2},\qquad W=B\rho,
\]

with the rational boundary point (B=1).  Converting this quartic to an
elliptic curve gives a way to generate exact first-square points rather than
waiting for two random integer squares.

The external run restricted (A\equiv\pm1\pmod {31}), converted 374 smooth
fibres, searched 23,044 elliptic points (including short multiples), and
recovered 11,819 distinct first-square points.  Of these, 379 also satisfy
the second square equation; all are exact collision points.

An internal run converted 1,447 fibres and searched 46,980 elliptic points.
After imposing an internal (p=31) collision residue, it likewise found no
nontrivial double-square point.

## 7. Reproducible artifacts

Code:

* `code/target_22224_p31_boundary_stable.m`
* `code/target_22224_p31_boundary_lifts.py`
* `code/target_22224_p31_boundary_crt_random.cpp`
* `code/target_22224_p31_boundary_fibration.m`
* `code/target_22224_qsquare_crt_sieve.m` (finite masks)

Principal logs/data:

* `results/target_22224_p31_boundary_stable.log`
* `results/target_22224_p31_boundary_lifts.log`
* `results/target_22224_p31_boundary_lifts_summary.tsv`
* `results/target_22224_p31_boundary_crt_random.log`
* `results/target_22224_p31_boundary_crt_random_h30000.log`
* `results/target_22224_p31_boundary_fibration.log`
* `results/target_22224_p31_boundary_fibration_internal.log`
* `results/target_22224_p31_boundary_crt_mask_p{43,47,53,59,61,67,71}_masks.tsv`

## 8. Next local step

The next (p=31) computation should introduce the first-order branch
separation explicitly on the 2,160 node-hit presentations.  For internal
thickness-one classes the regular model has no nontrivial component group,
so one needs the renormalized specialization of the marked divisor.  For the
external cubic-thick classes the component group has odd order three, so its
2-map is invertible; the remaining issue is again the renormalized
generalized-Jacobian coordinate.  The four sheet points and toric chambers
should be treated separately rather than mixed into that one-node chart.

In parallel, the negative fibration experiment suggests that leaving the
q-square slice and searching the full direct contact cover is now at least as
important as spending more trials on this particular slice.
