# Second-stage \([2,2,2,24]\) attack: exhaustive finite cover and boundary-aware height 187 sieve

Date: 2026-07-18

## Result

The second-stage computation found no rational \([2,2,2,24]\) realization,
but it produced a substantially stronger and correctly compactified search
result than the first height-20 experiment.

1. The q-square halving cover was exhaustively enumerated over
   \(\mathbf F_p\) for \(p=29,31,37,41\), with both \(\rho\)-sheets, both
   \(\sigma\)-sheets, and all three pairings defining the marked order-4
   class.
2. Every smooth/open marked class tested has exact order 12.  Finite abelian
   group decompositions independently determine whether it is divisible by
   2, including halves outside the affine \((M,N)\) chart.
3. Prime 31 has **no divisible marked class at any smooth affine slice
   point**.  This is not asserted as a global obstruction: it forces a
   rational target to reduce onto a base/bad-reduction boundary at 31.
4. A compiled, symmetry-complete, boundary-aware scan reached free-coordinate
   height 187.  It checked 95,427,692 ordered rational pairs, retained every
   reduction that lay either in a finite target mask or on a relevant
   boundary, and found no exact nontrivial double-square survivor.
5. The record point is an explicit positive control for the exact-square
   arithmetic: it occurs in the height-187 box, but is rejected because its
   reduction is smooth and nontarget at 31 (and also at 37).

Thus the final exact verifier had zero candidates and zero hits.  The result
is a negative bounded search on the q-square slice, not a nonexistence theorem
for the full \(A(2,2,2,12)\) family.

## 1. What was enumerated

Recall the slice

\[
 C=(AB)^{-1},\qquad
 A^2+B^2+C^2-3=\rho^2,
\]
\[
 A^{-2}+B^{-2}+C^{-2}-3=\sigma^2,
\]

with

\[
 s=(2\rho)^{-1},\quad (a,b,c)=s(A,B,C),\quad
 d=2s^2\sigma,
\]
\[
 t=s^2,\quad u=2t,\quad v=3t^2+d^2/2.
\]

For each nonzero \((A,B)\in\mathbf F_p^2\), the finite program does all of
the following.

* It retains both roots of each nonzero square \(R=\rho^2\) and
  \(S=\sigma^2\).
* It checks the direct contact identity and the curve discriminant.
* It uses all three partitions

  \[
  (a,b)\mid(c,d),\qquad (a,c)\mid(b,d),\qquad
  (b,c)\mid(a,d)
  \]

  in the marked order-4 construction.
* It checks that the CRT determinant and its leading coefficient are nonzero.
* It computes the exact order of every marked \(D_{12}\) in
  \(J(\mathbf F_p)\).
* It writes \(J(\mathbf F_p)\) as an explicit finite abelian group and tests
  whether the coordinate vector of \(D_{12}\) belongs to \(2J(\mathbf F_p)\).
* It retains every affine square-quartic solution \((M,N)\), reconstructs its
  Mumford half, and verifies

  ```text
  2*H = D12,  Order(D12) = 12,  Order(H) = 24.
  ```

The signs are exhaustive.  Varying the signs of \(A,B\) with
\(C=(AB)^{-1}\), followed by the \(\rho\)-sheet involution, realizes all
eight sign patterns on \((a,b,c)\); the \(\sigma\)-sheet gives both signs of
\(d\).  The raw labelled counts therefore deliberately contain symmetry
duplicates.  The log separately records unique marked Mumford classes and
unique halves.

## 2. Exact finite counts

The following table distinguishes raw labelled presentations from
symmetry-deduplicated classes.

| \(p\) | double-square \((A,B)\) | sheet-boundary \((A,B)\) | smooth slice sheets | open marked presentations | unique marked classes | target \((A,B)\) | target slice points | target presentations | unique target classes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29 | 280 | 4 | 384 | 1,152 | 64 | 12 | 24 | 72 | 4 |
| 31 | 244 | 4 | 192 | 576 | 32 | 0 | 0 | 0 | 0 |
| 37 | 408 | 4 | 768 | 2,304 | 128 | 12 | 24 | 72 | 4 |
| 41 | 520 | 52 | 768 | 2,304 | 128 | 36 | 72 | 216 | 12 |

Here “sheet boundary” means \(R=0\) or \(S=0\), so the affine formula for
\(s\) or the normalized fourth branch parameter leaves the open chart.
There were no contact failures and no CRT-chart failures on any smooth slice
sheet.

Every open marked presentation has order 12:

| \(p\) | order tests | order 12 | finite-group decompositions | divisible presentations |
|---:|---:|---:|---:|---:|
| 29 | 1,152 | 1,152 | 64 | 72 |
| 31 | 576 | 576 | 32 | 0 |
| 37 | 2,304 | 2,304 | 128 | 72 |
| 41 | 2,304 | 2,304 | 128 | 216 |

### Affine versus projective halving coordinates

Every smooth curve here has full rational 2-torsion, so a divisible marked
class has exactly 16 group-theoretic halves.  Comparing that number with the
affine \((M,N)\) solutions gives:

| \(p\) | affine retained halves | unique affine halves | projective \((M,N)\)-boundary halves | all group halves |
|---:|---:|---:|---:|---:|
| 29 | 720 | 40 | 432 | 1,152 |
| 31 | 0 | 0 | 0 | 0 |
| 37 | 1,152 | 64 | 0 | 1,152 |
| 41 | 3,240 | 180 | 216 | 3,456 |

There were no finite \(s_4=0\) solutions.  At 29 the affine chart sees 10 of
the 16 halves of each target presentation; at 37 it sees all 16; at 41 it
sees 15.  The missing multiplicities are retained in the boundary TSV files.

This comparison is essential at 31.  The absence of affine \((M,N)\) points
alone would not rule out a projective-chart escape.  The independent finite
abelian-group calculation does: none of the 576 smooth/open marked classes is
in \(2J(\mathbf F_{31})\).  Therefore a projective \((M,N)\) representation
cannot rescue a smooth affine base at 31.

## 3. Tangent ranks and Hensel scope

For each retained affine half, the program evaluates

\[
 \frac{\partial(E_1,E_0)}{\partial(M,N)}.
\]

All ranks are 2:

```text
p=29: [rank 0, rank 1, rank 2] = [0,0,720]
p=37: [rank 0, rank 1, rank 2] = [0,0,1152]
p=41: [rank 0, rank 1, rank 2] = [0,0,3240]
```

The two base equations have nonzero derivatives \(-2\rho\) and
\(-2\sigma\) in their sheet variables.  Consequently the full four-equation
Jacobian in \((A,B,\rho,\sigma,M,N)\) has rank 4 at every retained affine
target point.  These points are nonsingular and admit ordinary Hensel lifting
after choosing lifts of the two free base directions.

## 4. What the empty \(p=31\) target mask means

Of the 244 double-square \((A,B)\)-pairs modulo 31, four have
\(R=S=0\):

```text
(A,B) = (1,1), (1,-1), (-1,1), (-1,-1).
```

The other 240 split into 192 branch-collision pairs and 48 smooth pairs.  The
collision signatures are

```text
A^2=B^2 : 8       A^2=C^2 : 8       B^2=C^2 : 8
A^2=D^2 : 56      B^2=D^2 : 56      C^2=D^2 : 56,
```

where

\[
 D^2=(d/s)^2=(\sigma/\rho)^2=S/R.
\]

Each smooth pair has four sheet points and three pairings, producing the 576
open marked presentations.  None is divisible by 2.

This forces, but does not eliminate, the following p=31 boundary cases for a
rational target.

### A or B has a zero or pole

Let

\[
 \alpha=v_{31}(A),\quad\beta=v_{31}(B),\quad
 \gamma=v_{31}(C)=-\alpha-\beta.
\]

Outside the unit chamber, \(\rho\) is led by the most polar member of
\((A,B,C)\), while \(\sigma\) is led by the inverse of the member with largest
valuation.  After scaling the direct model, at least one nonzero branch
parameter specializes to zero and collides with the branch point at zero.
Thus this toric boundary can carry p-adic generic points, but only through
bad reduction of this direct branch model.  The rational sieve retains every
candidate for which a numerator or denominator of \(A\) or \(B\) is divisible
by 31.

### rho or sigma has positive valuation

After the isomorphism \(x=s^2X\), the normalized nonzero branch parameters
are

\[
 A,\ B,\ C,\ \sigma/\rho.
\]

If exactly one of \(\rho,\sigma\) vanishes modulo 31, the fourth parameter is
zero or infinite and the normalized special fiber is singular.  When both
vanish, the only residue points are the four listed above, where
\(A^2=B^2=C^2=1\), again forcing branch collisions.  These strata may still
support p-adic lifts with bad reduction, so the sieve retains them.

### Affine branch collision

The six collision divisors displayed above are retained.  A rational curve
on one of these residue divisors can be smooth over \(\mathbf Q\) while having
bad reduction at 31.

### Halving-chart boundary

There is no finite \(s_4=0\) solution over a smooth base.  More importantly,
the finite-group calculation proves that there is no group-theoretic half at
all over a smooth base, so neither finite nor projective \((M,N)\) coordinates
can escape the p=31 obstruction there.  A halving-boundary lift combined with
one of the singular base strata remains part of the unresolved bad-reduction
boundary; the bounded sieve retains the base stratum rather than discarding
it.

Accordingly, the phrase “p=31 obstruction” in this computation always means:

> a rational target must have p=31 boundary or bad reduction in the q-square
> coordinates.

It does not mean that the projective compactification has no rational or
p-adic points.

## 5. Boundary-aware rational sieve

The C++ scanner generates all positive reduced fractions with numerator and
denominator at most \(H\), orders the triple by

\[
 A\le B\le C=(AB)^{-1},
\]

and applies all six permutations and all product-\(+1\) sign patterns before
consulting a finite target mask.  Thus the modular mask does not accidentally
select only one labelling of the marked order-4 class.

At each prime, a pair is retained if either

* its reduction is smooth/open and one of the symmetry-related labelled
  pairs occurs in the exact finite target mask; or
* it lies on a zero/pole, sheet, or branch-collision boundary.

It is rejected if a purported rational square becomes a nonsquare, or if its
reduction is smooth/open but absent from the complete finite-group target
mask.  This is intentionally more conservative than requiring affine target
residues at 29, 37, and 41: bad-reduction candidates at those primes are also
retained.

The exact square test is particularly fast.  For

\[
 A=n_1/d_1,\qquad B=n_2/d_2,
\]

both square expressions have the common square denominator

\[
 (n_1d_1n_2d_2)^2.
\]

The program therefore needs only two signed 128-bit integer-square tests.

### Height 20

```text
positive rational values                 255
ordered pairs                         13,420
exact collision/trivial pairs            260
four-prime modular survivors                7
exact nontrivial double-square survivors    0
```

### Height 187

```text
positive rational values              21,407
ordered pairs                      95,427,692
exact collision/trivial pairs         21,464
four-prime modular survivors         982,078
exact nontrivial double-square survivors    0
```

The per-filter counts below are sequential: a pair rejected at an earlier
prime is not presented to a later filter.

| prime | nonsquare rejects | smooth/nontarget rejects | boundary passes | target passes |
|---:|---:|---:|---:|---:|
| 29 | 53,858,814 | 5,124,979 | 31,296,538 | 5,125,897 |
| 31 | 23,695,578 | 1,715,944 | 11,010,913 | 0 |
| 37 | 6,907,096 | 1,116,299 | 2,619,061 | 368,457 |
| 41 | 1,918,820 | 86,620 | 733,610 | 248,468 |

The empty p=31 target column is expected: every pair reaching that prime must
reduce to its retained boundary union.

### Record positive control

In canonical order, the record has

\[
 (A,B,C)=\left(\frac{13}{187},\frac{17}{7},\frac{77}{13}\right).
\]

The compiled exact arithmetic recognizes both square conditions.  Its mask
statuses are

```text
p=29 boundary
p=31 smooth_no_target
p=37 smooth_no_target
p=41 target.
```

Thus the known record is genuinely inside the height-187 enumeration and is
discarded for the mathematically correct reason, not through a height,
ordering, or integer-square bug.

The final candidate TSV contains only its header.  The exact Magma verifier
therefore reports

```text
candidate rows                 0
unique smooth presentations   0
exact divisibility tests       0
hits                           0.
```

Had a row survived, it would have been expanded into all 48 pairing/sign
presentations and tested for exact marked order 12, divisibility by 2, full
rational torsion, and a Frobenius root-power simplicity certificate.

## 6. Scope and next step

This establishes:

> There is no positive q-square-slice target with reduced free coordinates
> \(A,B\) of numerator/denominator height at most 187, after quotienting by
> the ordering \(A\le B\le C\).

The computation includes every sign and pairing of the direct marked class,
and it retains rather than suppresses bad-reduction residues at all four
primes.  It does not cover negative/complex orderings independently (their
sign data is already included in the marked presentations), other components
of the full direct \(A(2,2,2,12)\) family, or arbitrary points at larger
height.

The natural continuation is no longer a blind box enlargement.  It is a
blow-up analysis of the six p=31 collision divisors and the toric valuation
chambers, followed by a p-adic search on those bad-reduction components.  Any
q-square realization must enter one of them modulo 31.

## 7. Reproducible artifacts and commands

Code:

* `code/target_22224_qsquare_crt_sieve.m`
* `code/target_22224_qsquare_crt_sieve.cpp`
* `code/target_22224_qsquare_crt_exact.m`

Principal results:

* `results/target_22224_qsquare_crt_sieve.log`
* `results/target_22224_qsquare_crt_sieve_p29_points.tsv`
* `results/target_22224_qsquare_crt_sieve_p31_points.tsv`
* `results/target_22224_qsquare_crt_sieve_p37_points.tsv`
* `results/target_22224_qsquare_crt_sieve_p41_points.tsv`
* the corresponding `_masks.tsv` and `_boundary.tsv` files
* `results/target_22224_qsquare_crt_sieve_h20.log`
* `results/target_22224_qsquare_crt_sieve_h187.log`
* `results/target_22224_qsquare_crt_sieve_h187_candidates.tsv`
* `results/target_22224_qsquare_crt_exact_h187.log`

Commands from the repository root:

```bash
magma -b PrimeList:=29,31,37,41 \
  output_stem:=results/target_22224_qsquare_crt_sieve \
  code/target_22224_qsquare_crt_sieve.m

c++ -O3 -std=c++17 -o /tmp/target_22224_qsquare_crt_sieve \
  code/target_22224_qsquare_crt_sieve.cpp

/tmp/target_22224_qsquare_crt_sieve 187 \
  results/target_22224_qsquare_crt_sieve \
  results/target_22224_qsquare_crt_sieve_h187_candidates.tsv \
  results/target_22224_qsquare_crt_sieve_h187.log

magma -b \
  candidate_file:=results/target_22224_qsquare_crt_sieve_h187_candidates.tsv \
  log_file:=results/target_22224_qsquare_crt_exact_h187.log \
  code/target_22224_qsquare_crt_exact.m
```

