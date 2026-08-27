# Contact-6 `[6,12]`: weighted infinity in the cubic-contact variables

Date: 2026-07-10.

## Conclusion

The warning in `contact6_m612_boundary_core_modp.m` is substantive, but the
new weighted boundary is much more concentrated than the raw projective
closure suggests.

- Above every base-infinity direction with `B != 0`, including
  `[A:B:T]=[0:1:0]`, all simultaneous-pole initial forms are empty on the
  contact-open torus over both `F_5` and `F_7` in the bounded valuation-fan
  audit.  The only nontrivial candidate is an `x^5`-forced degeneration, and
  its exact saturated initial ideal is the unit ideal in both
  characteristics.
- The other endpoint, `[1:0:0]`, is different.  On the exact endpoint slice
  `a=1/e, b=0`, there are contact-open rational weighted initial strata at
  both primes.  Modulo `7`, four signatures have respectively `6,4,8,24`
  torus points and full Jacobian rank `5`.  Modulo `5`, one signature has
  `16` torus points; its open initial scheme has dimension `2` and every
  point has Jacobian rank `4`, so those points are smooth on the initial
  scheme but do not yet give a full-rank Hensel certificate for all five
  original equations.
- At every rational survivor at either prime, `DB` is a square and `DC` is a
  nonsquare.  Thus the discriminant square covers leave only the dual class
  `R3` as a possible half on this endpoint.  The `R2` and mixed `R1` lanes
  are locally excluded here.

So weighted poles shift the infinity priority from the ordinary-chart
exception `[0:1:0]` to `[1:0:0]`.  The next exact calculation is the next
weighted coefficient on the mod-`5` signature `E9`, followed by the exact
`R3` halving equations on any lifted branch.


## The uneliminated coefficient model

The eliminated equations in `(L,U,v)` acquire very large excess components
at infinity.  Retain instead the two coefficients `N,R` and write

```text
q = x^2 + U*x + v^2,
H = x^3 + N*x^2 + R*x + v^3.
```

The cubic-contact identity is exactly

```text
H^2 - q^3 = L^2*f.
```

Writing `c_i=[x^i]f`, coefficient comparison gives

```text
E5 = 2*N - 3*U                         - c5*L^2,
E4 = N^2 + 2*R - 3*U^2 - 3*v^2         - c4*L^2,
E3 = 2*v^3 + 2*N*R - U^3 - 6*U*v^2     - c3*L^2,
E2 = R^2 + 2*N*v^3 - 3*U^2*v^2-3*v^4  - c2*L^2,
E1 = 2*R*v^3 - 3*U*v^4                 - c1*L^2.
```

These five equations are equivalent to the three eliminated core equations
away from `2L=0`, and `p=5,7` are safe for this use.  More importantly, their
weighted leading forms retain the polynomial identity: if the `L^2*f` side
is too small, the initial equation is `H_0^2=q_0^3`.  Over an algebraic
closure this forces

```text
q_0=(x-r)^2,  H_0=(x-r)^3,
```

and hence `U_0^2-4v_0^2=0`.  Such a component is outside the contact-open
chart and is removed by saturation.


## Generic base infinity: `B != 0`

Normalize `B=1` and use

```text
[A:B:T]=[t:1:e],  a=t/e,  b=1/e.
```

For a generic direction the coefficient valuations are

```text
v(c5),v(c4),v(c3),v(c2),v(c1) = (-s,-2s,-2s,-2s,-s),
```

where `s=v(e)>0`.  Put

```text
(s,l,u,w,n,r)
 = (v(e),v(L),v(U),v(v),v(N),v(R)).
```

On the simultaneous-pole region `l,u,w<0`, exact minimum comparison in the
five coefficient equations produces three initial support signatures in the
bounded fan scan.  Representatives are

```text
G0  (1,-1,-4,-4,-4,-8),
G1  (1,-1,-3,-3,-3,-6),
G2  (1,-1,-3,-3,-2,-6).
```

`G0` is the pure `H_0^2=q_0^3` component.  The other two are the two faces
of the single `x^5`-forced cone.  If `s=1` and `v(L)=-d`, that cone has

```text
v(U)=v(v)=-(2d+1),
v(R)=-2(2d+1),
v(N) >= -(2d+1),                    d>=1,
```

and leading identity

```text
H_0^2-q_0^3 = (2/E)*ell^2*x^5,
```

where `e=p^s E+...` and `E` is a nonzero residue.  The code retains `E` as
a variable; setting it to `1` would lose weighted-projective residue orbits.

After saturation by

```text
E*ell*U*v*N*R*(U^2-4v^2),
```

all three initial ideals are the unit ideal over both `F_5` and `F_7`.
Thus no contact-open generic simultaneous-pole stratum survives.


## Why this does not contradict the ordinary infinity calculation

The ordinary calculation keeps `(L,U,v)` affine.  In the `B=1` chart its
first compatibility equation is

```text
t*(v^3-1)=0,
```

so `[0:1:0]` is the unique ordinary direction not immediately forced onto
`v^3=1`.  That observation remains correct.

For weighted poles, however, the leading coefficient pattern depends first
on whether `B` vanishes, not on whether `A` vanishes.  The point
`[0:1:0]` still has `B != 0`, so it has exactly the same `x^5`-forced pole
system as every other point in the `B=1` chart; that system is contact-open
empty at `p=5,7`.  At `[1:0:0]`, `B=0`, the valuations of the coefficients
of `f` jump and new initial forms appear.  Thus the ordinary and weighted
priorities concern different strata and are consistent.


## The exceptional weighted endpoint `[1:0:0]`

Use

```text
[A:B:T]=[1:t:e],  a=1/e,  b=t/e,
```

and first take the exact endpoint slice `t=0`, so `b=0`.  Then

```text
(c5,c4,c3,c2,c1)
 = (6, 2/e-15, 22, 1/e^2-15, 2/e+6),
```

with valuation pattern

```text
(0,-s,0,-2s,-s).
```

The bounded exact minimum scan gives the following ten support signatures.
The table reports the geometry of the initial scheme after contact-open
saturation.  Point counts include the leading unit `E` of `e`.

| id | representative `(s,l,u,w,n,r)` | `F_5`: dimension / points / ranks | `F_7`: dimension / points / ranks |
|:--|:--|:--|:--|
| `E0` | `(1,-2,-5,-5,-5,-10)` | empty | empty |
| `E1` | `(4,-1,-3,-3,-3,-6)` | empty | empty |
| `E2` | `(1,-2,-4,-4,-4,-8)` | empty | empty |
| `E3` | `(2,-1,-2,-2,-2,-4)` | empty | `1 / 6 / rank 5` |
| `E4` | `(6,-3,-2,-4,-6,-6)` | `1 / 0 / -` | `1 / 4 / rank 5` |
| `E5` | `(8,-4,-7,-6,-8,-13)` | empty | empty |
| `E6` | `(6,-3,-5,-4,-6,-9)` | `1 / 0 / -` | `1 / 8 / rank 5` |
| `E7` | `(6,-3,-5,-1,-6,-9)` | `1 / 0 / -` | `1 / 24 / rank 5` |
| `E8` | `(1,-2,-4,-4,-3,-8)` | empty | empty |
| `E9` | `(2,-1,-2,-2,-1,-4)` | `2 / 16 / rank 4` | empty |

For `E3,E4,E6,E7` modulo `7`, dimension `1` and rank `5` in the six torus
coordinates `(E,L,U,v,N,R)` mean every listed point is smooth; the one
dimension is the expected weighted scaling direction.  These points give
full-rank starting data for Hensel lifting on the exact `t=0` transformed
systems.

For `E9` modulo `5`, the saturated initial scheme has dimension `2` and rank
`4` at all `16` rational points, so those points are also smooth *on the
initial scheme*.  Because the five-equation Jacobian is not full rank, this
statement is deliberately weaker: one must compute the next weighted
coefficient to see whether the fifth equation obstructs a lift.  No claim of
a full `5`-adic core branch is made yet.

The ten signatures are the finite candidate list seen in the bounded
primitive integral weight scan used tonight; the script checks their exact
initial ideals, not a numerical approximation.  A mixed chart in which
`t` itself has a finite positive weight can refine the endpoint fan and is
not claimed to be exhausted here.


## `DB`, `DC`, and the surviving dual class

On the exact endpoint slice `a=1/e,b=0`, one has

```text
DB = e^-2*(1-6e-15e^2),
DC = e^-1*(-8-15e).
```

For either odd prime, `DB` is a square: `e^-2` is a square and its unit is
congruent to `1`.  Write

```text
e = p^s*E + higher terms,  E in F_p^*.
```

Then `DC` can be a square only if `s` is even and the residue `-8/E` is a
square.  Every rational initial survivor above has even `s`, but:

- on `E9` modulo `5`, `E` is `1` or `4`, and `-8/E` is nonsquare;
- on `E3` modulo `7`, `E` is `1`, `2`, or `4`, and `-8/E` is nonsquare;
- on `E4,E6,E7` modulo `7`, `E=1`, and again `-8/E` is nonsquare.

Therefore every surviving rational leading point has

```text
DB square,  DC nonsquare,  DB*DC nonsquare.
```

Using the exact dual-class labeling,

```text
R1 halved => DB*DC square,
R2 halved => DC square,
R3 halved => DB square,
```

only `R3` remains possible on this endpoint.


## Reproduction and next tests

Run

```text
magma -b code/contact6_m612_weighted_infinity_leading.m
```

The script prints every minimum-index signature, its representative weight,
the saturated open dimension, rational point count, Jacobian ranks, samples,
and the `DC` leading squareclass split by `E`.

The next focused calculations are:

1. substitute the `E9` weights
   `(2,-1,-2,-2,-1,-4)` into the exact endpoint equations over `Z_5` and
   compute the next nonzero weighted coefficient;
2. if any of the `16` initial points passes, impose the exact `R3` halving
   equations (the `DB` square cover is already automatic);
3. treat a finite positive weight for `t` to determine whether mixed paths
   approaching `[1:0:0]` add further signatures;
4. use the full-rank mod-`7` signatures as the cleaner positive controls for
   the weighted Hensel implementation.
