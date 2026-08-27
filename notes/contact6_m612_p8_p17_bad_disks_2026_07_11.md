# The ten forced `p=17` disks on the `P8` extra-3 cover (2026-07-11)

## Result

The ten projective parameter residues retained by the good-reduction sieve
fall into three exact geometric types:

| parameter keys | `e mod 17` | source event |
|---|---:|---|
| `2,14=-3` | `0` | endpoint `e=0` |
| `1,3,11,15` | `8` | `a=-2`; the two quadratic factors share `x=1` |
| `0,8,12,infinity` | `4` | `a=13`; `DC=0`, so the `C` quadratic has a double root |

The first weighted/local layer gives:

- The two `e=0` disks have no E3 branch.  Their leading `e`-units are
  squares, whereas the E3 system over `F_17` has points precisely for
  nonsquare leading units.  The same squareclass excludes every uniformly
  scaled E3 depth.  Other endpoint weight signatures were not exhaustively
  enumerated at 17.
- The four `a=-2` disks have no affine signed contact point even on the raw
  special-fiber equations.  This is a certified affine obstruction.  A
  nonintegral weighted branch could still exist; no unbounded pole fan was
  attempted.
- The four `DC=0` disks each contain two signed, orthogonal, fixed-parameter
  rank-5 branches.  Hence every punctured disk has points on the desired
  degree-12 support factor and its degree-24 signed lift.

Thus `p=17` does not obstruct the `P8` extra-3 route globally: the `DC=0`
disks support smooth local branches.

## Exact classification

As before, set

```text
A=tau^2+tau-6,
B=tau^2+6,
e=-400*A^2*B^2/(768*A^4-1200*A^2*B^2+1250*B^4),
a=1/e.
```

For `b=0`, the source polynomial is

```text
f=x*(3*x^2+(a-3)*x+2)*(2*x^2-3*x+(a+3)).
```

Its discriminant support includes

```text
a+2,  DB=a^2-6*a-15,  DC=-8*a-15.
```

Direct evaluation on `P^1(F_17)` gives

```text
key 0,8,12,infinity: e=4,  a=13, DC=0;
key 1,3,11,15:       e=8,  a=15=-2, a+2=0;
key 2,14:             e=0.
```

There are no `tau^2+6=0` parameter-pole residues modulo 17 and no pole of
the rational `e` map among these keys.

At every finite collision key the derivative `de/dtau` is nonzero.  In the
order shown above, the exact derivative data are

```text
(0,13), (1,1), (3,12), (8,12),
(11,3), (12,9), (15,1).
```

In reciprocal coordinate `w=1/tau`, the derivative at infinity is `7`.
Therefore all eight collision disks are etale over their corresponding
`e`-disk.  Each contains one exact singular parameter; away from it the
source curve is smooth.

## The two endpoint disks

Substitute `tau=tau0+17*s` and write `e=17^2*E`.  Exact reduction gives

```text
tau0= 2: E= 4*s^2 mod 17,
tau0=-3: E=15*s^2 mod 17.
```

Both `4` and `15` are squares modulo 17, so every primitive layer in either
disk has square leading `E`.

The E3 equations are the same five equations derived in the `p=7` disk
analysis.  Efficient complete enumeration over `F_17` gives points exactly
for

```text
E in {3,5,6,7,10,11,12,14},
```

the eight nonsquare classes.  Each such `E` has two signed points and the
fixed-variable `5 x 5` Jacobian has rank 5.  Consequently neither endpoint
disk meets E3 at its first layer.  More generally, for
`tau=tau0+17^n*s` with `s` a unit, the leading unit is still `4*s^2` or
`15*s^2`, so every uniformly scaled E3 layer is excluded.

This does not classify different, non-E3 valuation signatures.  In
particular, no claim is made about an unbounded endpoint fan at `p=17`.

## The four common-root disks

At `a=-2 mod 17`, the special source factors as

```text
f = x*(x-1)^2*(3*x-2)*(2*x-1).
```

The script enumerates every tuple `(L,U,V)` with

```text
L*V*(U^2-4*V^2) != 0,
```

recovers `N,R` from the first two coefficient equations, and tests all five
coefficients of

```text
H^2-q^3=L^2*f.
```

There are zero raw solutions, before even imposing support-disjointness or
the Weil-pairing condition.  Hence the affine signed contact chart is empty
in all four common-root disks.

This is not yet a full local obstruction.  A hypothetical point could have
nonintegral `(L,U,V,N,R)` and appear only after a weighted blowup of the
node.  That pole fan was deliberately left unresolved rather than inferred
from the affine calculation.

## The four `DC=0` disks

At `a=13 mod 17`, one has

```text
f = 6*x*(x-2)*(x-5)^2*(x-6),
```

where `C=2*(x-5)^2`.  Complete enumeration finds exactly two signed open
solutions:

```text
(L,U,V,N,R)=( 1,10,11,1,5),
(L,U,V,N,R)=(-1,10,11,1,5).
```

Equivalently,

```text
q=x^2+10*x+2,
H=x^3+x^2+5*x+5.
```

The support is disjoint from the singular branch locus.  At both points the
Jacobian of the five exact coefficient equations with respect to
`(L,U,V,N,R)` has rank 5.  Since the parameter map is etale, fixing any
17-adic parameter in a collision disk and applying multivariate Hensel gives
the corresponding signed lift; at a punctured parameter the lifted curve is
smooth and the class has exact order 3.

The explicit Weil-pairing expression from the `p=7` analysis evaluates to
`1` at both signs.  Also, `Q_17` contains no nontrivial cube root of unity,
so any pair of rational 3-classes is necessarily orthogonal.  These are
therefore branches of the degree-12 orthogonal support and its degree-24
signed cover, not of the degree-27 nonorthogonal orbit.

## Status by disk

| disks | certified conclusion | unresolved part |
|---|---|---|
| `tau=2,-3` | no E3 branch at any scaled depth | other weighted endpoint signatures |
| `tau=1,3,11,15` | affine signed contact chart empty | nonintegral/node blowups |
| `tau=0,8,12,infinity` | two smooth signed orthogonal branches | no local obstruction remains |

The positive `DC=0` branches mean that another prime or global geometry is
needed to close `P8`.

## Reproduction

From `torsion_jac` run

```bash
magma -b code/contact6_m612_p8_p17_bad_disks.m
```

The script finishes in about a second and asserts the disk classification,
etaleness, endpoint squareclass split, complete E3 enumeration, affine
common-root emptiness, and the rank/pairing certificates on the `DC=0`
branches.
