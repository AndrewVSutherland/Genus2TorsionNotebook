# Z/35 next route: blow up the simultaneous contact equations

This continues the Z/35 lane after the open contact-7 plus 5 route and the
shallow `p=3` boundary search looked cold.  I did **not** run another broad
height search in `(a,b)`.  The better continuation is option **B**: analyze the
interrupted simultaneous contact7/contact5 point equations structurally, then
use the small-prime obstruction to choose local charts.

## Setup

The contact-7 family is

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3,
f = (h^2 + (x - 1)^7)/x^2.
```

The old point-contact 5 route imposed

```text
q^2 - f = const*(x-r)^5,        q = c0 + c1*x + c2*x^2.
```

The top equation is

```text
c2^2 - b^2 = 5*r - 7.
```

Writing

```text
d = c2 - b,     e = c2 + b
```

gives

```text
r  = (d*e + 7)/5,
b  = (e-d)/2,
c2 = (e+d)/2.
```

The previous enumerator then divided by `b` and `c2` to solve for `a` and
`c0`, leaving two residual equations

```text
N0(d,e,c1) = N1(d,e,c1) = 0.
```

So the route is not really a 3-parameter search; it is a curve in the
`(d,e,c1)` space, with bad charts where `b=0` or `c2=0`.

## Probe

Code:

```text
code/agent_Z35_next_route_probe.m
```

Run:

```text
magma code/agent_Z35_next_route_probe.m
```

The residual equations have:

```text
N0 degree/terms 16 165
N1 degree/terms 11 49
```

Eliminating `c1` gives:

```text
Res_c1(N0,N1)
  = d^3 * e^3 * (d-e)^8 * (d+e)^4 * Phi38(d,e)
```

where `Phi38` has total degree `38` and `280` terms.  Thus the old search is
controlled by one large residual plane curve plus four obvious degenerate
factors:

```text
d=0, e=0, d=e, d=-e.
```

Modulo `3`, the nondegenerate residual chart has no points:

```text
POINT_CONTACT_FINITE_TABLE
p  residual_solutions  smooth_open  good_pass_5  bad_or_boundary
3  0                   0            0            0
7  4                   2            2            2
11 2                   0            0            2
13 4                   2            2            2
17 10                  6            6            4
19 8                   6            6            2
23 28                  24           24           4
29 28                  22           22           6
31 26                  20           20           6
```

All `F_3` solutions of `N0=N1=0` are forced onto five degenerate centers:

```text
(d,e,c1) = (0,1,1)  tag D=0
(d,e,c1) = (2,0,2)  tag E=0
(d,e,c1) = (2,1,0)  tag D=-E
(d,e,c1) = (1,1,1)  tag D=E,r=1
(d,e,c1) = (2,2,2)  tag D=E,r=1
```

The first two have `b,c2` units and both map to the already-seen contact-7
boundary residue

```text
(a,b) = (2,2),  h(1)=0,
f = (x+2)^2 * (x^3 + 2*x^2 + 2*x + 2)  over F_3.
```

The last three are chart poles for the old parametrization:

```text
d=-e   means c2=0,
d=e    means b=0,
r=1    means the 5-contact point collides with the 7-contact point modulo 3.
```

These pole charts are exactly where dividing by `b` and `c2` loses structural
control.  They are better continuation targets than another `(a,b)` height
box.

A small 3-adic lift count also says the centers are not immediately dead:

```text
center       mod 9 lifts   mod 27 lifts   comment
(0,1,1)     3             9              maps to old (2,2), h(1)=0 branch
(2,0,2)     3             9              maps to old (2,2), h(1)=0 branch
(2,1,0)     9             27             c2=0 pole chart
(1,1,1)     27            243            b=0, r=1 pole chart
(2,2,2)     27            243            b=0, r=1 pole chart
```

## Interpretation

The point-contact route is **dead in the nondegenerate mod-3 chart**.  A
rational solution cannot reduce to an ordinary point of the residual
`(d,e,c1)` curve at `p=3`; it must either be non-integral at `3` or pass through
one of the five centers above.

It is not yet globally dead.  The useful next route is a local one:

1. Blow up the two `d=e, r=1` centers first.  These have the largest lift
   counts and are also the `b=0` pole of the previous formulas.
2. Then check the `d=-e` center, the `c2=0` pole.
3. Treat `d=0` and `e=0` as lower priority, because they land on the same
   `(a,b)=(2,2)` visible `h(1)=0` branch that already looked cold, though now
   with the 5-contact equations imposed.

Concretely, the next script should return to the original coefficient
equations in `(a,b,c0,c1,c2,r)` and use local charts such as `b=3*B` near
`d=e`, and `c2=3*C2` near `d=-e`, rather than solving by division through those
quantities.  After dividing by the correct powers of `3`, run a finite-field
viability table on the transformed equations.  If those transformed charts
have no viable `F_3` directions or no good-prime 5-divisibility survivors, the
point-contact route can be called genuinely dead.

I would only move to the degree-2 contact analogue after these five forced
point-contact centers are killed; otherwise the degree-2 setup adds variables
before we have used the obstruction already present in the smaller system.
