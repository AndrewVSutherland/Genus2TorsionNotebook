# Rational points on the genus-6 trigonal quotient

## Outcome

The genus-6 quotient admits a reasonably compact integral bidegree `(10,3)`
model.  Independent projective searches in its `z`- and `rho`-projections,
each through height `25000`, found six rational places represented by five
plane points.  All six are exact boundary points.  None gives a rational
lift to an open point of the original genus-12 root cover.

The two searches each cover `759,920,359` primitive rational coordinate
values.  Thus `1,519,840,718` independently sieved projection values were
tested with 76 projective finite-prime masks.  The only survivors were

```text
z-survivors:    2, 5, 14/3, 32/7
rho-survivors: -1, 0, 1.
```

This is a bounded exclusion, not yet a determination of all rational points
on the genus-6 curve.

## An explicit integral model

Let

```text
d = 5*z^2 - 20*z - 16,
```

and put

```text
PA = 8461353*z^10 - 322506900*z^9 + 5471086032*z^8
     - 54340244352*z^7 + 349528189824*z^6
     - 1519391001600*z^5 + 4514219274240*z^4
     - 9038973468672*z^3 + 11658061873152*z^2
     - 8735943819264*z + 2886221168640,

PB = -445299*z^8 + 13749864*z^7 - 183283224*z^6
     + 1375090368*z^5 - 6337639680*z^4 + 18329174016*z^3
     - 32393066496*z^2 + 31885688832*z - 13343981568,

PC = 7813*z^6 - 187500*z^5 + 1849992*z^4 - 9583360*z^3
     + 27411072*z^2 - 40903680*z + 24776704.
```

The quotient coefficients from the contact-30 family simplify to

```text
A = PA/d^5,       B = PB/d^4,       C = PC/d^3.
```

Consequently an integral plane equation is

```text
F(z,rho) = 2*d^5*rho^3
           + (PA-3*d^5)*rho^2
           + d*(PB+3*d^4)*rho
           + d^2*(PC-d^3) = 0.                 (1)
```

After bihomogenizing separately in `[Z:W]` and `[X:Y]`, equation (1) has
bidegree `(10,3)`.  Its function field is exactly the quotient constructed
in `code/contact30_c3root_trigonal_geometry.m`, and Magma verifies that its
normalization has genus `6`.

The denominator divisor `d=0` has no rational base value because
`disc(d)=720` is nonsquare.  It therefore introduces no rational point into
the cleared affine model.

## Exact rational fibers and normalization

The plane equation is singular at several relevant specializations, so it is
important to decompose fibers in the normalized function field.  The exact
result is:

| `z` | plane `rho` factorization | degrees and ramification upstairs | rational places |
|---:|---|---|---:|
| `2` | `2(rho-1)^3` | `(1,1),(1,2)` | 2 |
| `5` | `2 rho^2(rho+1)` | `(1,1),(1,2)` | 2 |
| `14/3` | `2(rho-1)^3` | `(1,3)` | 1 |
| `32/7` | `2(rho-1)^3` | `(1,1),(2,1)` | 1 |

Thus there are six normalized rational places represented by the five plane
points

```text
(2,1), (5,-1), (5,0), (14/3,1), (32/7,1).
```

The two places over `(2,1)` coalesce in this nonmaximal plane order.  At
`z=32/7`, only one of the three normalized branches is rational; the other
place has residue degree two.

The `rho`-fiber factorizations independently recover exactly the same plane
points:

```text
rho=-1: rational z-root 5;
rho= 0: rational z-root 5;
rho= 1: rational z-roots 2, 32/7, 14/3.
```

At `z=infinity`, the exact fiber is

```text
6250*rho^3 + 8451978*rho^2 - 2217120*rho + 192200,
```

which is irreducible over `Q`.  At `rho=infinity`, equation (1) requires
`d=0`, which has no rational root.  Neither projective infinity chart adds a
rational point.

## The quadratic lift to the genus-12 cover

The quotient invariant was

```text
z = (3*R^2-7)/(3*R-5),
```

so a quotient point lifts rationally only if

```text
3*R^2 - 3*z*R + 5*z - 7 = 0.                  (2)
```

For the four surviving `z`-values:

| `z` | discriminant of (2) | rational `R`-lifts | classification |
|---:|---:|---|---|
| `2` | `0` | `R=1` | `c=0` boundary |
| `5` | `9` | `R=2,3` | `c=0` boundary |
| `14/3` | `0` | `R=7/3` | `c=0` boundary |
| `32/7` | `-108/49` | none | no rational lift |

The last point is also intrinsically a boundary of the descended open locus.
Over `Q(sqrt(-3))`, its two lifts are

```text
R = (16 +/- sqrt(-3))/7.
```

For either lift, `t=(1 +/- sqrt(-3))/2`, `u=t^3=-1`, hence `c=0`; both the
numerator and denominator in the recovery of `q` vanish.  The extra quotient
point is therefore not concealing an open curve over a quadratic field.

## Projective search

The search used the bihomogeneous form of (1).  For each of 76 primes
`7 <= p <= 401`, it formed two exact masks:

1. a projective `z` residue is allowed if some `rho in P1(F_p)` satisfies
   (1);
2. a projective `rho` residue is allowed if some `z in P1(F_p)` satisfies
   (1).

Because this is direct reduction of an integral bihomogeneous equation, the
test remains a necessary condition at bad or singular fibers; no potentially
valid bad disk is discarded.  The exhaustive boxes were

```text
z=a/b:    b>0, gcd(a,b)=1, |a|<=25000, b<=25000;
rho=a/b:  b>0, gcd(a,b)=1, |a|<=25000, b<=25000.
```

Each box contains `759,920,359` primitive rationals.  After all 76 masks,
only the coordinate values listed above remained.  Exact factorization and
normalization then classified all of them as boundary.

This quotient search is genuinely different from the earlier height search
in `R`: it also sees rational quotient points with no rational `R`-lift.  The
new `z=32/7` boundary point is an example.

## Reproduction

```bash
c++ -O3 -std=c++17 code/contact30_g6_quotient_projective_sieve.cpp \
    -o /tmp/contact30_g6_quotient_projective_sieve
/tmp/contact30_g6_quotient_projective_sieve 25000 quiet
magma -b code/contact30_g6_quotient_points_verify.m
```

Artifacts:

- `code/contact30_g6_quotient_projective_sieve.cpp`
- `code/contact30_g6_quotient_points_verify.m`
- `data/contact30_g6_quotient_projective_h25000.txt`
