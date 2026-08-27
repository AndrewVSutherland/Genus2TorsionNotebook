# Orbit-12 CRT pullback and first-cover fibration

Date: 2026-07-11.

This note pulls the certified local orbit-12 branches at `11`, `19`, and
`23` back to the complete rational Clebsch--Klein chart, records the
first-radicand double cover, and runs a focused CRT-lattice search.  The
marked class is `{r1^2,r2^2}` and

```text
G0 = -(r1^2-r3^2)(r1^2-r4^2)(r1^2-r5^2).
```

## The chart and its first double cover

The complete labelled smooth chart is

```text
R1 = 1+t*(t+2)*m
R2 = t*m*(m-t-2)
R3 = -1+m+t*(t+1)*m^2
R4 = 1+t-m-t*m^2
R5 = -(1+t)*(1+t*m^2).
```

The exact pullback of `G0` factors into eleven reduced components:

```text
m, m-1, t, t+1, t-m+2, t*m+1, t*m+m-1, t*m+t+1,
t^2*m+t*m^2+2*t*m-t+m,
t^2*m^2-t^2*m+t*m^2-2*t*m+m-2,
t^2*m^2+t^2*m+t*m^2+2*t*m+t+2.
```

Thus `Y^2=G0(t,m)` has degree `12` in `m`; the direct `t`-fibers are
generically genus `5`, not elliptic.  A better genus-one direction is the
marked-root ratio

```text
q = r2/r1.
```

Normalize `r1=1`, put `r2=q`, `r3=x`, and call the other roots `u,v`.
The CK equations give

```text
u+v = -(1+q+x),
u*v = (q+x)*(1+q)*(1+x)/(1+q+x).
```

The condition that `u,v` be rational is the quartic genus-one equation

```text
V^2 = (q+x+1)*(
    q^3-q^2*x-q^2-q*x^2-2*q*x-q+x^3-x^2-x+1
).
```

The symmetric expression for the first radicand is

```text
G0 = -(1-x^2)*((1+u*v)^2-(u+v)^2).
```

It gives a sextic squareclass in `x` for generic `q`.  The exact symbolic
derivation and all factorization checks are in
`code/elkies22210_orbit12_first_cover_fibration.m`.

## How many local residues reach the chart

Projecting the certified boundary-lift states to the integral inverse-chart
locus (`r1=1` and `r1+r2` a unit) gives:

| prime power | all truncated cover states | distinct `(t,m)` | resolved states | distinct resolved `(t,m)` | resolved density |
|---:|---:|---:|---:|---:|---:|
| `11^3` | 182,289 | 109,247 | 29,040 | 22,086 | 1.2467% |
| `19^2` | 13,410 | 6,940 | 2,166 | 1,482 | 1.1372% |
| `23^2` | 19,128 | 10,614 | 3,174 | 3,174 | 1.1342% |

Here “all truncated” includes deep zero radicands, while “resolved” means
that all four valuations are already finite and even.  The reproducible
projection is

```text
python3 code/elkies22210_orbit12_chart_residue_counts.py
```

These counts are useful diagnostics, but they are not independent residue
conditions on `t` and `m`: the allowed objects are pairs `(t,m)`.

## The corrected selected CRT disk

Use the three certified primitive local seeds

```text
11^3: (1,242,959,9,120)
19^2: (1,248,98,3,11)
23^2: (1,392,118,8,10).
```

At `11` and `23`, their inverse chart scales are the units `485` and `420`.
At `19`, however, the inverse parameters are

```text
(t,m) = (41,347) mod 19^2 = (3,5) mod 19,
```

and all five polynomial chart coordinates vanish modulo `19`.  The chart
therefore meets a base point: using only `(t,m) mod 19^2` loses one digit of
the primitive CK tuple.

Fixing the first three coordinates of the seed, the CK Hensel equations give
the correction digits `(6,12)` in the last two coordinates.  Hence a lift is

```text
r = (1,248,98,2169,4343) mod 19^3,
(t,m) = (2929,347) mod 19^3.
```

Evaluating the chart modulo `19^3` and dividing its common factor `19`
gives

```text
(12,88,93,36,132) = 12*(1,248,98,3,11) mod 19^2,
```

which certifies that the lost digit has been restored.

The corrected parameter moduli and residues are

```text
moduli: (11^3,19^3,23^2) = (1331,6859,529)
(t,m):  (968,1202), (2929,347), (178,511).
```

Their product and CRT classes are

```text
M = 4,829,415,041,
t = 4,214,241,119 mod M,
m = 2,991,361,145 mod M.
```

The Gauss-reduced lattices `n-X*d=0 mod M` have bases

```text
t: (-44407,33027), (67155,58808)
m: (-52348,35765), (60221,51112).
```

For the genus-one fiber parameter, the same disk gives

```text
q mod (11^3,19^3,23^2) = (242,248,392),
q = 2,539,833,076 mod M,
q lattice basis = (38357,-22386), (-45012,-99637).
```

The shortest reduced fiber parameter in this disk is

```text
q = -38357/22386.
```

This is the natural first fiber for a subsequent Mordell--Weil or descent
calculation.

## Pullback of the first-cover square root

Choose the printed local `G0` roots.  The `19^2` root `97` lifts uniquely in
the chosen tangent disk to `3346 mod 19^3`.  After applying the chart scales,
the root data are

```text
(22,112,303),
```

where the middle entry denotes `Y/19^3 mod 19^2`, since the chart has a
common coordinate factor `19`.  Put globally `Z=Y/19^3`.  Then

```text
Z mod (1331,361,529) = (770,112,424).
```

The eight independent local sign choices give eight CRT classes modulo
`254,179,739`; the search code prints all of them.  This provides a weighted
`(t,m,Z)` lattice seed for work directly on the first double cover.

## Focused exact search

Run

```text
python3 code/elkies22210_orbit12_chart_crt_lattice.py --height 2000000
```

It enumerates reduced numerators and positive denominators bounded by
`2,000,000` in the selected `t` and `m` lattices.  The exact result is

```text
t fractions                         841
m fractions                         852
parameter pairs                 716,532
smooth CK tuples                716,532
G0 positive                     372,223
all four radicands positive     127,000
all four passing small-mod sieve  3,670
G0 exact rational squares             0
full orbit-12 cover hits               0
```

Every all-positive row has exact square mask `0`: none has even one exact
square radicand.  The detailed run record is
`data/elkies22210_orbit12_chart_crt_lattice_h2000000_summary.txt`.

This negative result is complete only for **one selected product of local
disks**, including one chosen tangent lift at `19`, and only to the stated
rational parameter height.  It does not eliminate the many other local
branches tabulated above.

## The special fiber `q=-2` is impossible

The only nonboundary rational specialization found by the branch-collision
calculation is `q=-2`.  On it the CK and first-cover square equations reduce
to elliptic quartics

```text
V^2 = (x-1)*(x^3+x^2-x-9),
U^2 = -x*(x-2)*(x-1)*(x+1).
```

However the **full** orbit-12 cover has a clean `Q_23` obstruction on this
fiber.  Direct enumeration modulo `23` shows that its only six reductions
are the permutations of

```text
(r3,r4,r5) = (0,2,-1).
```

Near the displayed reduction write

```text
r3=A,  r4=2+B,  r5=-1-A-B.
```

The CK equation has linear part `-3*A+9*B`, so at the first nonzero
`23`-adic digit an open lift must have `A/B=3`.  The radicand attached to
`r4` is

```text
G4 = (r4^2-4)*(1-r3^2)*(1-r5^2)
   = -32*B^2 + higher terms.
```

After removing its even valuation, its leading unit is

```text
-32 = -9 mod 23,
```

a nonsquare because `-1` is a nonsquare modulo `23`.  Permuting the three
unmarked roots gives the same argument at all six reductions.  Therefore
the `q=-2` fiber has no smooth `Q_23` point on the full orbit-12 cover and
cannot produce `[2,2,2,20]`.

