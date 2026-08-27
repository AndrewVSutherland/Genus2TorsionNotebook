# Memory-bounded modular triage of the `b2=0` M(12)+5 cover

## Scope and resource policy

This follows the degree-one-`B` task in
[`m12_general5_order60.md`](m12_general5_order60.md) without forming the
six-variable characteristic-zero Groebner basis.

Every finite enumeration was capped at 2 GiB virtual memory and 180 seconds.
Every Magma ideal computation was capped externally at 8 GiB and internally
at 6 GB, with a maximum wall timeout of 300 seconds.  In practice the finite
enumerations took between a fraction of a second and 20 seconds.  The only
successful primary decomposition, on one fixed finite-field fiber, took less
than a second after formula construction.

The reproducible code and preserved output are:

```text
code/m12_general5_b2zero_modular_points.py
code/m12_general5_b2zero_modular_triage.m
data/m12_general5_b2zero_modular_triage_2026_07_13.log
```

## The two loci

For `b1 != 0`, use the root-of-`B` coordinates from
[`m12_general5_b2zero_rootquotient.md`](m12_general5_b2zero_rootquotient.md).
The sign quotient is defined by four equations in

```text
(b,w,c,d,e).
```

Their degrees and term counts modulo 7 are

```text
(11,99), (13,244), (12,106), (23,1610).
```

The generic open conditions are imposed directly: the compact base is open,
`d != 0`, `e^2 != 4`, `ell3 != 0`, `F` is a squarefree quintic, and
`gcd(Z^2+e*Z+1,F(dZ-c))=1`.  The signed order-5 cover exists over a quotient
point precisely when

```text
tau = -5*(2-e)^3/(8*ell3)
```

is a nonzero square.

The constant-`B` locus is not obtained by saturating the linear chart.  It is
treated separately by putting `s=b0^2` in

```text
A^2 - s*F = (x^2+U*x+V)^5.
```

After forcing the five high coefficients of `A`, this gives five equations
in `(b,w,U,V,s)`.  A quotient point lifts to constant `B` exactly when the
nonzero value `s` is a square.

## Exact finite-field census

The generic quotient results are:

| `p` | raw zeros | cheap open | full open | signed points | Jacobian ranks | rational base fibers |
|---:|---:|---:|---:|---:|---|---|
| 7  | 612   | 1  | 1  | 2  | one rank 4 | one fiber of size 1 |
| 11 | 2460  | 2  | 2  | 2  | two rank 4 | two fibers of size 1 |
| 13 | 4125  | 3  | 3  | 2  | three rank 4 | three fibers of size 1 |
| 17 | 9341  | 5  | 5  | 4  | five rank 4 | five fibers of size 1 |
| 19 | 13063 | 10 | 10 | 10 | ten rank 4 | eight fibers of size 1, one of size 2 |

Here `signed points` counts both signs of `lambda`, so it is twice the number
of quotient points with square `tau`.  The sequence `2,2,2,4,10` agrees with
the independent counts in the earlier root-quotient implementation.

Every fully open point has the expected maximal equation-Jacobian rank 4.
Thus the generic quotient is locally alive and smooth of dimension one at
all 21 displayed points.  It is, however, very sparse at the small primes.
Except for one base at `p=19`, distinct rational quotient points have distinct
`(b,w)` projections.  This rational-point fiber statement is not a geometric
degree computation.

The constant-`B` results are:

| `p` | raw zeros | cheap/full open | signed points | rank |
|---:|---:|---:|---:|---|
| 7  | 343  | 0 | 0 | -- |
| 11 | 1332 | 1 | 2 | 5 |
| 13 | 2197 | 0 | 0 | -- |

The unique full-open quotient point modulo 11 is

```text
(b,w,U,V,s) = (10,1,7,5,3),
```

and is reduced and isolated.  At 7 and 13 all raw solutions disappear as
soon as the cheap chart factors are inverted.  Consequently any rational
constant-`B` solution would have to enter a denominator or boundary chart at
both primes.  This is not a global nonexistence proof.

## Bounded component experiment

The four-variable quotient obtained after eliminating `e` was worse than the
five-variable chart:

```text
K degrees/terms:       (30,2431), (20,546), (30,1970)
Res(q,F) degree/terms: (70,62650)
```

The modulo-7 raw dimension call made no progress in about 150 seconds and was
interrupted below the 8 GiB cap.  On the five-variable chart the resultant
falls to degree 24 with 1969 terms, but even the first global saturation was
not cheap enough to justify continuing under this triage budget.  No global
primary-decomposition or component-degree claim is made.

As a controlled zero-dimensional test, fix the unique open modulo-7 base
`(b,w)=(5,6)`.  Full saturation and primary decomposition give one component,
dimension zero, with a Groebner basis consisting of five independent linear
polynomials.  It is the reduced rational point

```text
(b,w,c,d,e) = (5,6,2,1,4).
```

Magma's `Degree` intrinsic rejected the nonhomogeneous affine ideal, but the
linear basis and exhaustive point count certify degree one for this fiber.
This does not lift a component over `Q`; it only verifies that the finite
fiber machinery is seeing the expected reduced point.

## Rational reconstruction check

The unique signed quotient residues at `p=7,11,13` were combined modulo
1001.  At the natural numerator/denominator bound 22, only `w=-9/2`,
`d=-4/17`, and `e=-8/19` reconstruct; `b`, `c`, and `tau` do not.  Adding
either signed residue at `p=17` and raising the natural bound to 92 still
does not reconstruct all coordinates simultaneously.  Hence there is no
low-height rational point singled out by these small-prime signed orbits.
This CRT exercise is heuristic because finite points at different primes
need not be reductions of the same rational point.

## Decision

* **Squarefree/coprime open of generic linear `B`: locally go, monolithic
  Groebner no-go.**  This explicitly saturated open quotient is smooth of
  dimension one at every counted point, and its signed cover is locally alive
  at every tested prime.  The current coordinates do not expose a recurring
  low-degree component under the memory/time cap.  A next attempt should use
  structural projection, function-field elimination, or modular interpolation
  one equation at a time, not global saturation of the raw ideal.

* **Constant `B`: deprioritize.**  It is forced into boundary/bad-reduction
  charts at 7 and 13 and has only one isolated open quotient point at 11.
  There is no evidence here for a rational component.

* **No normalization was attempted.**  The only component recovered was the
  reduced point in a fixed modulo-7 fiber.  Treating it as a characteristic-
  zero component would be unjustified.

The cyclic-60 route therefore remains locally viable on the tested
squarefree/coprime open of the generic linear-`B` curve, but this triage found
no rational order-60 example and no small component ready for reconstruction.
The excluded repeated-support and shared-Weierstrass boundary loci were not
decided by these finite counts and remain separate live or obstructed tasks.
