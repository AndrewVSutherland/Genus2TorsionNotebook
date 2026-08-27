# Deterministic marked identities and transverse tangent at the HLP `[60]` seed

## Fixed generator

Magma's image of `G.1` from `TorsionSubgroup` is not a deterministic choice
of generator: different Magma seeds can return different generators of the
cyclic group.  All calculations here therefore start with the fixed published
Mumford class

```text
P60 = (x^2 + 74/15*x + 128/25,
       -1081510/9*x - 1226624/3,
       2)
```

on

```text
f = -46250000*x^6 + 1761500625*x^4
    -22332312000*x^2 + 94277468160.
```

The deterministic primary multiples are

```text
D5 = 12*P60 = (x^2-4608/395, -(164280/79)*x, 2), order 5,
D3 = 20*P60 = (x^2-316/25,  -5476,             2), order 3,
D4 = 15*P60 = (x^2-1506/125,-24642/5,          2), order 4,
T2 = 30*P60 = (x^2-1728/125,0,                 2), order 2,
2*D4 = T2.
```

These equalities, the exact orders, `J(Q)_tors=[60]`, and all identities
below are checked by `code/hlp_z60_marked_identities_verify.m`.

## Exact low-degree identities

For the order-5 layer put

```text
q5 = x^2-4608/395,
B5 = x^2-1728/125,
A5 = (37/320)*x*(52705*x^4-1294080*x^2+7962624),
k5 = 341553260289/4096.
```

Then

```text
A5^2-f*B5^2 = k5*q5^5,
A5 = B5*v5 (mod q5),       v5=-(164280/79)*x.
```

Thus `A5-B5*y` selects `D5` and its norm is a fifth power.  The useful
structural feature is that `B5=q0`, where `q0=x^2-1728/125` is the support
of the order-2 class `T2`.

For the order-3 layer put

```text
q3 = x^2-316/25,
H3 = 29600-2775*x^2,
k3 = 46250000.
```

Then

```text
H3^2-f = k3*q3^3,
H3 = v3 (mod q3),           v3=-5476.
```

For the order-4 layer put

```text
q0   = x^2-1728/125,
u4   = x^2-1506/125,
ell4 = -2775*q0,
k4   = 46250000.
```

Then

```text
ell4^2-f = k4*q0*u4^2,
ell4 = -v4 (mod u4),         v4=-24642/5.
```

This is the compact halving identity for `2*D4=T2`.

## Normalized marked incidence and its differential

Use all seven coefficients of a degree-at-most-six sextic as base
coordinates.  The three normalized systems are

```text
A^2-f*B^2-k*q^5 = 0,        q and B monic;
H^2-f-k*q^3 = 0,            q monic;
(q0*L)^2-f-k*q0*u^2 = 0,    q0 and u monic.
```

Their auxiliary variable/equation counts are respectively `11/11`, `7/7`,
and `7/7`.  Exact rational differentiation at the HLP point gives

```text
auxiliary block ranks = 11, 7, 7,
combined rank         = 25 in 32 variables,
tangent dimension     = 7,
projection to df      = rank 7.
```

All three auxiliary determinants are nonzero.  Consequently every sextic
first-order deformation has a unique lift through these three normalized
marked identities.

## An exact transverse lift

Take the genuinely odd perturbation

```text
df = x.
```

The unique first derivatives in the monic-`B5` normalization are as follows;
all omitted derivatives vanish:

```text
order 5:
  dq5 = (2209/172868558400)*x,
  dA5 = 720396/6331625
        -(2246299/121567200)*x^2
        +(14177669/18672721920)*x^4,
  dB5 = -(43/1013060000)*x,
  dk5 = 0;

order 3:
  dq3 = -(3/1620896000)*x,
  dH3 = -(137/202612)*x +(75/1620896)*x^3,
  dk3 = 0;

order 4:
  dq0 = -(1/145880640)*x,
  du4 = (29/4668180480)*x,
  dL  = -(8125/175056768)*x,
  dk4 = 0.
```

Substitution in the 25 differentiated coefficient equations gives zero
exactly.

Magma computes `#Aut(C)=4`, so the reduced automorphism group has order two:
there is one non-hyperelliptic involution and hence one local Humbert-4
branch.  In the affine sextic coefficient chart, a primitive normal to its
tangent hyperplane is

```text
N = (0,81125,0,904800,0,9916416,0)
```

on `(df0,...,df6)`.  It annihilates the even-sextic directions and the two
odd infinitesimal `PGL(2)` directions.  Since

```text
N(x)=81125 != 0,
```

the lifted tangent is genuinely transverse to the geometrically split
locus, not merely a coordinate change.

## What this proves, and what remains arithmetic

This is a positive result: the split HLP point is not trapped in the split
locus by the simultaneous marked order-5, order-3, and order-4 equations.
There is no first-order deformation obstruction, and the formal implicit
function theorem gives a unique algebraic formal lift along `f+t*x` over
`Q[[t]]`.

It is not yet a new rational cyclic-`[60]` Jacobian.  A rational value of
`t` must make the finite order-5 and order-4 algebraic lifts rational
simultaneously.  Invertibility of the tangent blocks is expected from the
finite-etale nature of prime-to-characteristic torsion and does not by
itself supply such a specialization.  The next arithmetic task is therefore
to study the exact one-dimensional pullback above a simple transverse base
curve, not to enlarge blind coefficient boxes.

## Reproduction

```text
magma -b MemGB:=2 code/hlp_z60_marked_identities_verify.m
python3 code/hlp_z60_marked_tangent.py
```

The Magma verifier used about 83 MB before the automorphism check; both jobs
stay far below the requested 2 GB cap.
