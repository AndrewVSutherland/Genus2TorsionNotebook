# Contact-6 `[6,12]`: exact dual halving on the affine mod-5 cones

Date: 2026-07-10.

## Result

The corrected affine blowup audit left 16 certified smooth points on the
necessary dual square covers:

```text
DB=W^2, U+2v=0:  8 cones, tested against exact halving of R3;
DC=W^2, U+2v=0:  8 cones, tested against exact halving of R2.
```

All lie over

```text
(a,b)=(0,0) mod 5.
```

None survives as a smooth point of the valid leading square-quartic halving
chart.  There is an unresolved degenerate boundary at `m=0`; it is not an
exact half and is not counted as a surviving cone.

## Exact specialization

On the distinguished Richelot dual, write

```text
Delta = a*b+3*a+3*b+5,
R1 = (b^2-2*a-3)x^2 +(2*a*b+6*a+6*b+10)x +(a^2-2*b-3),
R2 = 2x^2-(a+3),
R3 = 2-(b+3)x^2.
```

For the class `Ri`, put

```text
S = Ri*(m*x+n)^2 - Delta*Rj*Rk.
```

If `s4=[x^4]S` is nonzero, `S` is a scalar square of a quadratic exactly
when

```text
E1 = 8*s4^2*s1 - s3*(4*s4*s2-s3^2) = 0,
E0 = 64*s4^3*s0 - (4*s4*s2-s3^2)^2 = 0.
```

At `(a,b)=(0,0)` modulo 5,

```text
Delta=0,
R1=R2=R3=2*(x^2+1).
```

Consequently, for both `R3` and `R2`, direct exact expansion gives

```text
s4 = 2*m^2,
E1 = 4*m^5*n,
E0 = -m^6*(m^2+n^2).
```

On the valid chart `s4 != 0`, hence `m != 0`.  Then `E1=0` forces `n=0`,
while

```text
E0=-m^8 != 0.
```

Thus the valid exact-halving fiber is empty over `F_5`.  This obstruction
is independent of the four projective normal directions and the two values
of `L`, so it applies to all eight `DB` cones and all eight `DC` cones.

## The degenerate `m=0` boundary

The common zero set of the two displayed covariants consists only of

```text
(m,n)=(0,n), n in F_5.
```

At every such point

```text
s4=0
```

and the `2 x 2` covariant Jacobian in `(m,n)` has rank zero.  The
square-quartic equivalence explicitly assumes `s4 != 0`, so these five
residues per core cone are degree-drop artifacts of this affine chart, not
halves and not Hensel-smooth exact-halving points.

A simultaneous weighted blowup of `m=0`, the vertical factor
`Delta=0 mod 5`, and the square-cover normal could in principle contain a
pole/degree-drop halving branch.  That tree was not pursued here.  The
precise conclusion is therefore:

```text
certified smooth exact-halving cones on the valid leading chart: 0;
degenerate m=0 boundary: unresolved.
```

## Reproduction

Run

```text
magma -b code/contact6_m612_affine_dual_exact_leading_mod5.m
```

The calculation derives the covariants, checks the displayed identities,
enumerates all `(m,n)` in `F_5^2`, verifies that the valid chart is empty,
and records rank zero on the degenerate boundary.
