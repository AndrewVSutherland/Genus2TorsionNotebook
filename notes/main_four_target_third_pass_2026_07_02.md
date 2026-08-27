# Four target third pass

Date: 2026-07-02.

This continues `main_four_target_second_pass_2026_07_02.md`.  The goal was to
push all four missing-group lanes past the second-pass stopping points, using
separate agents where useful.

## Executive summary

The third pass did not produce a new rational genus-2 example, but it changed
the shape of the search:

```text
Z/5 x Z/5      full degree-2 Mumford norm is live over finite fields;
               the old B=1 slice was too small.

Z/2 x Z/24     the old A(2,12) fiber is bad, but small split fibers have
               degree-4 saturated halving branches.

Z/48           one new A16 gate survivor was certified, but its torsion is
               [16], not [48]; the remaining RTHeight=4 partitions are
               reproducible.

Z/35           the b=0,r=1 pole chart is genuinely alive 3-adically; the
               c2=0 pole dies in finite original variables over F_3.
```

## `Z/5 x Z/5`

Detailed files:

```text
code/agent_z5x5_full_mumford_norm.m
notes/agent_z5x5_full_mumford_norm.md
```

The full degree-2 Mumford norm system on the one-contact-5 family is:

```text
f = h^2 - K*x^5,
A(x)^2 - B(x)^2*f = Lambda*U(x)^5,
U = x^2 + s*x + t,      deg(A)<=5, deg(B)<=2.
```

On the open normalized chart

```text
h = 1 + h1*x + h2*x^2,
A = x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0,
B = b0 + b1*x + b2*x^2,
U = x^2 + s*x + t,
```

the coefficients `a4,...,a0` are forced by the degree `9,...,5` coefficients
of

```text
A^2 - B^2*(h^2 - K*x^5) - U^5.
```

The five residual equations in `(h1,h2,K,s,t,b0,b1,b2)` have sizes:

```text
E0: degree 30, terms 852
E1: degree 27, terms 560
E2: degree 24, terms 347
E3: degree 21, terms 208
E4: degree 18, terms 123
```

The height-2 normalized rational box was cold:

```text
rational_tested=300000
independent_hits=0
bad_U_boundary=128
```

All rational norm solutions in that box were the double-root `U=x^2`
boundary.  However, this is not evidence that the full chart is empty:
finite-field full-norm hits exist with `deg(B)=2`, outside the discarded
`B=1` slice.  One recorded hit over `F_11` is:

```text
h = 2*x^2 + 7*x + 1
K = 6
U = x^2 + 9*x + 2
V = 10*x
B = 4*x^2 + 2*x + 7
#J(F_11)=100
```

Magma verified the contact 5-torsion class, the degree-2 class, and their
independence over `F_5`.

Next move: keep the full norm system.  The useful attacks are:

```text
1. Saturate/eliminate on K != 0, disc(f) != 0, disc(U) != 0,
   gcd(B,U)=1, and (b1,b2) != (0,0).

2. Use the F_7/F_11 full-norm points as local lifting or modular-obstruction
   charts.

3. Treat b2=0 as a smaller side branch; its residual equations have degrees
   10,9,8,7,6 rather than 30,27,24,21,18.
```

## `Z/2 x Z/24`

Detailed files:

```text
code/agent_A2_24_branch_factor_scan.m
notes/agent_A2_24_branch_factor_scan.md
```

The previous split fiber

```text
p=-5/3, z=1, r=2/3
```

stays irreducible of saturated degree 16 after removing the `s4=0` boundary.
The new scan enumerated small rational split A(2,12) fibers and translated
the visible order-12 class by rational 2-torsion classes.

Height-4 summary:

```text
checked=11132
split_fibers=29
order12_split_fibers=29
translated_order12_rows=116
low_rows=40
errors=0
```

Every row had:

```text
raw resultant degree in M = 32
s4-boundary contribution = [<2,8>]
saturated affine degree = 16
gcd(E1,E0)=1
```

The useful change is factorization.  The best small fibers are:

```text
p=-1/3, z=-1, r= 4/3
p=-1/3, z= 1, r= 4/3
p= 1/3, z=-1, r=-4/3
p= 1/3, z= 1, r=-4/3
```

For each, the two extra 2-torsion translations have saturated factor degrees

```text
[<4,1>, <4,1>, <8,1>].
```

No rational factor appeared in the height-4 box, but degree 4 is a much better
target than the original irreducible degree-16 fiber.

Next move: extract these quartic branches explicitly, test their resolvents
and rational points, and only then decide whether a Chabauty/MW-sieve style
argument is appropriate.  The old `p=-5/3,z=1,r=2/3` branch should be
deprioritized.

## `Z/48`

Detailed files:

```text
code/agent_Z48_next_scan.m
notes/agent_Z48_next_scan.md
```

The new driver is a partitioned version of the sign-aware A(8)->A(16)
square-root scan.  It skips the completed RTHeight=3 box and runs a reproducible
slice of the RTHeight=4 boundary:

```text
ExcludeRTHeight:=3
SliceMod:=4
SliceClass:=0
```

Production result for `SliceClass:=0`:

```text
rawSlices=462 completedSkipped=182 eligibleSlices=280
partitionSkipped=210 runSlices=70
tested=1129030 commonRootPairs=15 rationalRoots=15
singular=9 nonsingular=6
pointGateReject=5 pointGatePass=1
exactTried=1 certified=1 z48Hits=0
```

The five nonsingular point-count rejects were first killed at:

```text
p=5 : 3
p=7 : 2
```

There was one new gate survivor:

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144,
N=5/8, z=125/96
gate primes: <13,192>, <17,336>, <19,384>
running gcd: 48
```

Exact A16 certification succeeded, but the integral square model had torsion
invariants:

```text
[16]
```

so this is an A16 false positive for the `48 | #J(F_l)` gate, not a Z/48
example.

A local check of the old completed RTHeight=3 box with the weaker 3-part-only
gate also killed all four nonsingular roots at `p=5`:

```text
tested=2935478 commonRootPairs=37 rationalRoots=37
singular=33 pointGateReject=4 pointGatePass=0 z48Hits=0
FIRST_POINT_GATE_KILLS: 5 : 4
```

Next move: run the remaining reproducible RTHeight=4 partitions:

```text
SliceClass:=1
SliceClass:=2
SliceClass:=3
```

with the same `RTHeight:=4 ExcludeRTHeight:=3 SearchBound:=10 PrimeBound:=43
MinGood:=3 SliceMod:=4` parameters.  If those are cold, the next genuinely new
route is to add the cubic-contact equations after the point-count gate rather
than continuing broad A8 prefilters indefinitely.

## `Z/35`

Detailed files:

```text
code/agent_Z35_b0_pole_blowup.m
notes/agent_Z35_b0_pole_blowup.md
```

The original point-contact equations are:

```text
q^2 - f = -(x-r)^5,       q=c0+c1*x+c2*x^2,
f=(h^2+(x-1)^7)/x^2,     h=1-(7/2)*x+a*x^2+b*x^3.
```

Using the original coefficient equations, with no division by `b` or `c2`,
the two `b=0,r=1` centers over `F_3` are:

```text
(a,b,c0,c1,c2,r) = (1,0,1,1,1,1)
(a,b,c0,c1,c2,r) = (1,0,2,2,2,1)
```

For both charts

```text
a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1, c2=t+3*C2, r=1+3*R,
```

every transformed equation has 3-adic content exactly `3`.  After dividing by
one power of `3`, the special-fiber equations are linear.  The lift table is:

```text
center t   F3 directions   liftable mod 9   total mod 9 lifts   liftable mod 27   total mod 27 lifts
1          27              9                243                 9                 6561
2          27              9                243                 9                 6561
```

Thus the `b=0,r=1` pole chart definitely does not die at first or second
3-adic order.  The liftable directions look like 3-dimensional smooth local
branches after the first obstruction.

Good-prime point-contact enumeration also found smooth 5-divisibility
survivors at every tested prime:

```text
p   total   smooth   pass 5 | #J(F_p)
7   7       2        2
11  10      6        6
13  6       2        2
17  16      10       10
19  12      8        8
23  40      32       32
29  38      30       30
31  30      22       22
```

The `c2=0` pole center from the old `d=-e` branch, by contrast, has no finite
original-variable `F_3` lifts.

Next move: parameterize the nine liftable first directions for each center,
push them one or two more 3-adic orders, and impose the genus-2 smoothness and
5-contact/Jacobian tests on those reduced branches.  This is now a better
Z/35 target than the nondegenerate residual chart.

## Smoke checks from the main thread

I also ran small local checks after closing the agents:

```text
Z/5 x Z/5: rational_height=1, finite_trials=100, cantor_trials=50
  ran cleanly; no tiny random hit.

Z/2 x Z/24: Height=1, MaxRows=1
  reproduced the known p=-5/3,z=1,r=2/3 row with saturated degree 16.

Z/48: RTHeight=2, ExcludeRTHeight=1, SearchBound=2, SliceClass=0
  ran cleanly; all three rational roots in this tiny partition were singular.

Z/35: full b=0 pole script
  reproduced the two viable F_3 centers, lift table, and good-prime table.
```

## Priority after this pass

The highest-yield next actions are:

```text
1. Z/2 x Z/24: extract and study the explicit quartic halving branches from
   the four best A(2,12) fibers.

2. Z/35: push the eighteen total liftable b=0,r=1 first directions beyond
   mod 27 and turn them into branch equations.

3. Z/48: finish SliceClass 1,2,3 for RTHeight=4 outside RTHeight=3.

4. Z/5 x Z/5: begin elimination/saturation on the full norm chart, starting
   with b2=0 as the smaller side branch.
```
