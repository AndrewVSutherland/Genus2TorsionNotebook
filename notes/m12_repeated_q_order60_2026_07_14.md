# Repeated support on compact `M(12)`: exact point-order-10 reduction

Date: 2026-07-14.

## Exact confluent identity

Write the compact model as

\[
 F=LQ,\qquad L=b+(2b-1)x,
\]

where

\[
 Q=LH^2+4b(1+x)^2(wL-x^2),\qquad H=x+w(1+bx).
\]

For a non-Weierstrass point `P=(u,-L(u)C(u))`, the repeated Mumford
support `q=(x-u)^2` is valid on the open set

```text
b*w*(b-1)*(2*b-1)*s*L(u)*C(u) != 0
```

when a quadratic `C` and nonzero `s` satisfy the polynomial identity

\[
 Q(x)=L(x)C(x)^2+2s(x-u)^5.                 \tag{1}
\]

This is stronger than the boundary norm heuristic.  Indeed

\[
 (y-LC)(y+LC)=2sL(x-u)^5,
\]

and `y+LC` has divisor

\[
 5P+W_L-6\infty.
\]

Therefore

\[
 5(P-\infty)=W_L-\infty,
\]

so `P-infinity` has exact order 10.  The repeated class
`[(x-u)^2,v]=2(P-infinity)` has exact order 5.  Thus the repeated-support
and shared-Weierstrass boundaries are two descriptions of the same
point-order-10 locus.

The marked compact `M(12)` class `R` has order 12 and
`6R=W_L-infinity` when the rational 2-rank is one.  Hence `R` together
with `P-infinity` generates a cyclic group of order 60; equivalently an
appropriate sum has order 60.

An exact norm reconstruction is also immediate from (1): put

\[
 B=C/s,\qquad A=(x-u)^5+LC^2/s.
\]

Then `A^2-B^2 F=(x-u)^10` identically.

## Direct finite-field census

`code/m12_repeated_q_order60.m` tests the relation by Jacobian arithmetic,
not merely by the norm identity.  Exact order-10 point counts on smooth
compact fibers were:

| p | exact order-10 points | examples with irreducible quartic cofactor |
|---:|---:|---:|
| 7 | 2 | yes |
| 11 | 4 | yes |
| 13 | 2 | no in the printed sample |
| 17 | 6 | yes |
| 19 | 12 | yes |
| 23 | 10 | no in the printed sample |
| 29 | 16 | yes |
| 31 | 18 | yes |
| 37 | 30 | yes |
| 61 | 42 | yes |

Thus the full quadratic-function locus is locally alive and has many
fibers compatible with rational 2-rank one.

## Linear-`B` repeated slice

Specializing the earlier root-of-`B` quotient to `e=-2` gives the
repeated support `(Z-1)^2`.  After imposing a square signed parameter and
checking the resulting points in the Jacobian, the open quotient counts
were

```text
p=7:  0
p=11: 1
p=13: 0
```

The lone `p=11` point has exact point order 10, but its quintic factors
as `1+2+2`, so it has extra rational 2-torsion over the residue field.
The linear-`B` sublane is therefore not a useful global search lane; the
full quadratic function is essential.

## Reciprocal parameter and resultant

After moving `W_L` to `z=0`, write

```text
G = z*C(z)^2 + k*(z-r)^5.
```

The square cover is rationally reduced by

```text
r  = -b/t^2,
c2 = 2*(b-1)*t^5/b,
k  = -c2^2,
```

with the top coefficients determining `c1,c0`.  The two residual
equations have shapes

```text
N1: total degree 28, degrees (b,w,t)=(12,4,18), 143 terms
N3: total degree 31, degrees (b,w,t)=(9,2,24),  74 terms.
```

Their `w`-resultant factors as

```text
t^30 * b^4 * (b-1)^10 * (b-1/2)^8
* ((5*t^8-1)*b^2 - 10*t^8*b + 5*t^8)^4
* R(b,t).
```

Here `R` is irreducible over `Q`, has total degree 80, bidegree `(22,58)`,
has 446 terms, and is even in `t`.  The quadratic factor has discriminant

\[
 20t^8,
\]

so it has no rational point for nonzero rational `t`.  Consequently every
open rational candidate lies on the single irreducible factor `R`.

## Bounded rational search

For each nonzero reduced rational `t` of height at most 50, the exact
degree-22 polynomial `R(b,t)` was solved over `Q`; this checks **all**
rational `b`, without bounding its height.  The 3,094 slices produced no
rational root.  Eighteen additional representative small slices were
factored completely and likewise produced only boundary roots, the
pointless quadratic above, and an irreducible factor of degree at least
12.

No cyclic-60 curve was obtained.  The exact construction is locally alive,
but the remaining problem is rational points on a high-degree irreducible
curve rather than a low-degree parametrization.
