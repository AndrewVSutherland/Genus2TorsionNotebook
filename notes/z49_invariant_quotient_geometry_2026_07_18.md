# Finite geometry of the order-49 sign quotient

The iterated contact-7 incidence equations in `(a,b,u,v,r)` are invariant
under

```
(u,v) -> (-u,-v).
```

This note examines the quotient itself rather than performing another
rational-height search.

## Two affine quotient charts

On `u != 0`, put

```
w=u^2,  s=v/u.
```

A monomial `u^i*v^j` becomes `w^((i+j)/2)*s^j`; every occurring `i+j`
is even.  On `v != 0`, the root-centered chart is

```
t=v^2,  q=-u/v,
B^2=t*(x-q)^2.
```

Here `u^i*v^j` becomes `(-1)^i*t^((i+j)/2)*q^i`.  The root chart sharply
reduces the visible slope degree: the four equations have total degrees
`12,10,9,7`, compared with `24,21,18,15` in the slope chart.  The missing
`v=0` locus is the already understood constant-`B` chart `s=0`.

## Complete finite enumeration

For fixed `(a,b,t,q,r)`, the norm identity is equivalent to

```
A^2 = t*(x-q)^2*f + (x-1)*(x-r)^7,
```

where `A` is a monic quartic.  Its coefficients are recovered from degrees
7 through 4, after which degrees 0 through 3 are checked.  This gives a
short exact enumeration and includes nonsquare values of `t`, unlike an
enumeration in `(u,v)`.  The constant chart uses the same test with
`A^2=w*f+(x-1)*(x-r)^7`.

The complete open quotient counts over `F_p` are

```
p             5   11   13   17   19   23
quotient      2    5    6    9   13   15
signed lifts  4   10   12   18   26   30
```

The first four signed counts reproduce the independent incidence counts
given for the original system.  Interestingly, every quotient point at all
six primes has square `t` or `w`, so no nonsquare twist occurs in this
sample.  This observation is not asserted as a theorem in characteristic
zero.

## Projection interpolation

For each of `(q,t)`, `(s,b)`, and `(s,r)`, form the evaluation matrix of
all plane monomials of total degree at most `d`.  At both `p=19` and
`p=23`, the ranks for `d=1,2,3` are exactly

```
3, 6, 10,
```

the full monomial counts.  Thus the finite quotient data show no linear,
conic, or cubic plane relation in any proposed projection.  The same data
have full rank for bidegrees `(2,2)`, `(1,3)`, and `(3,1)`.

At `p=23` each projection has 14 distinct points.  The degree-4 matrix has
rank 14 among 15 monomials, so its one-dimensional nullspace is the
tautological quartic interpolant through 14 points.  Without agreement at
another prime or an elimination certificate, it should not be interpreted
as a rational projection equation.

Consequently a useful plane projection, if one exists in these coordinates,
has total degree at least 4 and lies outside all three tested small
bidegree rectangles.  This is a finite-reduction diagnostic rather than a
global elimination proof: no flatness/saturation certificate was computed.

## Components and genus

The open point counts are too small and omit unknown boundary points, so
they do not determine a zeta function, component count, or genus.  They are
compatible with a single connected high-degree curve, but do not prove it.
A defensible genus computation would first require a projective saturated
model or a certified finite map with controlled boundary fibers.

Artifacts:

```
code/z49_invariant_quotient_finite.py
data/z49_invariant_quotient_finite.txt
```

The largest run used about 85 MB, below the 180 MB cap.
