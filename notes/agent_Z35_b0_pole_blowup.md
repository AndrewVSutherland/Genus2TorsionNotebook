# Z/35 b=0 pole blow-up

Date: 2026-07-02.

This is the third-pass local check after `agent_Z35_next_route.md`.  I returned
to the original coefficient equations for

```text
q^2 - f = -(x-r)^5,       q=c0+c1*x+c2*x^2,
f=(h^2+(x-1)^7)/x^2,     h=1-(7/2)*x+a*x^2+b*x^3.
```

Using integral equations:

```text
G4 = c2^2 - b^2 - 5*r + 7
G3 = 2*c1*c2 - 2*a*b - 21 + 10*r^2
G2 = c1^2 + 2*c0*c2 - a^2 + 7*b + 35 - 10*r^3
G1 = 2*c0*c1 + 7*a - 2*b - 35 + 5*r^4
G0 = 4*c0^2 - 8*a + 35 - 4*r^5
```

No formula dividing by `b` or `c2` is used in the local chart.

Code:

```text
code/agent_Z35_b0_pole_blowup.m
```

Run:

```text
magma code/agent_Z35_b0_pole_blowup.m
```

## b=0, r=1 charts

The two finite original centers over `F_3` above the old `d=e,r=1` residual
centers are

```text
(a,b,c0,c1,c2,r) = (1,0,1,1,1,1),
(a,b,c0,c1,c2,r) = (1,0,2,2,2,1).
```

For `t=1,2` I used

```text
a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1, c2=t+3*C2, r=1+3*R.
```

In both charts every transformed equation has 3-adic content exactly `3`, so
each equation is divided by one power of `3`.  The first special-fiber equations
are linear.

For `t=1` they are:

```text
C2 + 2*R + 2
B + 2*C1 + 2*C2 + 2*R
A + B + 2*C0 + 2*C1 + 2*C2
A + B + 2*C0 + 2*C1 + 2*R + 2
A + 2*C0 + R
```

For `t=2` they are:

```text
C2 + R + 2
B + C1 + C2 + 2*R + 2
A + B + C0 + C1 + C2
A + B + C0 + C1 + 2*R + 1
A + C0 + R + 1
```

Finite direction/lift table:

```text
center t   F3 directions   Jacobian ranks   liftable mod 9   total mod 9 lifts   liftable mod 27   total mod 27 lifts
1          27              all 3            9                243                 9                 6561
2          27              all 3            9                243                 9                 6561
```

The lift histogram is the same for both centers:

```text
mod 9:  18 directions have 0 lifts, 9 directions have 27 lifts
mod 27: 18 directions have 0 lifts, 9 directions have 729 lifts
```

Thus the `b=0,r=1` pole chart definitely does **not** die at first or second
3-adic order.  The liftable directions look like honest 3-dimensional smooth
branches after the first obstruction: each surviving direction has `3^3` lifts
to mod `9` and `3^6` lifts to mod `27`.

The liftable `t=1` first directions are:

```text
<A,B,C0,C1,C2,R>
<0,2,0,1,1,0>
<0,2,1,2,2,1>
<0,2,2,0,0,2>
<1,1,0,2,0,2>
<1,1,1,0,1,0>
<1,1,2,1,2,1>
<2,0,0,0,2,1>
<2,0,1,1,0,2>
<2,0,2,2,1,0>
```

The liftable `t=2` first directions are:

```text
<A,B,C0,C1,C2,R>
<0,2,0,2,2,2>
<0,2,1,0,0,1>
<0,2,2,1,1,0>
<1,1,0,1,0,1>
<1,1,1,2,1,0>
<1,1,2,0,2,2>
<2,0,0,0,1,0>
<2,0,1,1,2,2>
<2,0,2,2,0,1>
```

## Good-prime table

I also enumerated the original point-contact equations over good primes, with
the `b=0` and `c2=0` poles handled by direct enumeration rather than division.
For `p != 3`, the translated `b=3*B` chart is just an affine reparameterization,
so this is a full original point-contact finite-field table, not a local
3-adic specialization.

```text
p   total point-contact solutions   smooth contact-7 curves   5 | #J(F_p)   b=0 points   c2=0 points
7   7                               2                         2             1            1
11  10                              6                         6             2            0
13  6                               2                         2             0            0
17  16                              10                        10            2            0
19  12                              8                         8             0            0
23  40                              32                        32            2            0
29  38                              30                        30            4            4
31  30                              22                        22            0            0
```

Every smooth finite-field point-contact solution in this table passed the
necessary `5 | #J(F_p)` test.  So there is no good-prime 5-divisibility
obstruction visible in this finite table; the point-contact condition is doing
what it should at smooth good reductions.

## c2=0 pole check

The old `d=-e` residual center is

```text
(d,e,c1) = (2,1,0),
```

which forces, in original finite variables,

```text
b=1, c1=0, c2=0, r=0.
```

Directly testing all finite `a,c0` over `F_3` in the original equations gives:

```text
finite_original_F3_lifts = 0.
```

Moreover, with `a,r` integral at this center, the equation `4*G0=0` contains
`4*c0^2` plus integral terms, so a simple `v_3(c0)<0` polar escape cannot
cancel the lowest valuation.  I did not pursue a larger nonintegral
multi-variable chart here.

## Verdict

The `b=0,r=1` pole chart is **viable**, not dead.  Both F3 centers have
first-order directions, and exactly one third of those directions lift cleanly
through mod `9` and mod `27`.  The good-prime finite table also has smooth
5-divisibility survivors at every tested good prime.

The finite original `c2=0` pole center, by contrast, dies immediately over
`F_3` in original variables.
