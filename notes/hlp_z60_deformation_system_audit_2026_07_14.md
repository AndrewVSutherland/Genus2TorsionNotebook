# Audit of the simultaneous marked deformation at the HLP `[60]` seed

## Deterministic base point

The torsion generator returned as `G.1` by `TorsionSubgroup` is not a
deterministic choice.  The deformation system is therefore based on the
published class

```text
P60 = (x^2+74/15*x+128/25,
       -1081510/9*x-1226624/3,
       2)
```

on

```text
F = -46250000*x^6+1761500625*x^4
    -22332312000*x^2+94277468160.
```

Its primary multiples are

```text
D5 = 12*P60 = (x^2-4608/395, -(164280/79)*x, 2),
D3 = 20*P60 = (x^2-316/25, -5476, 2),
D4 = 15*P60 = (x^2-1506/125, -24642/5, 2),
D2 = 30*P60 = (x^2-1728/125, 0, 2).
```

Magma checks their exact orders `5,3,4,2`, respectively, and `2*D4=D2`.

## What carries over from the old `[6,12]` calculation

The order-3 contact block and the order-4 halving block from
`notes/m612_hlp_deformation.md` carry over unchanged in form:

```text
H^2-F=k3*q3^3,
(q2*L)^2-F=k4*u4^2*q2.
```

Only their seed coordinates change.  The new ingredient is the general
degree-two Mumford norm for `D5`:

```text
A^2-B^2*F=k5*q5^5.
```

The restricted identity `H^2-F=k*q^5` would miss this class.  In the
normalization `B(0)=1`, the exact data are

```text
q5 = x^2-4608/395,
B5 = 1-(125/1728)*x^2,
A5 = (48752125/110592)*x^5
     -(1558625/144)*x^3+66600*x,
k5 = 1778923230671875/4076863488.
```

Here `A5+B5*v5=0 mod q5`, for `v5=-(164280/79)*x`, so the selected
principal function corresponds to the fixed published `D5`.

For the other two layers one may take

```text
q3=x^2-316/25,       H3=29600-2775*x^2,    k3=46250000,
q2=x^2-1728/125,     u4=x^2-1506/125,
L4=2775,             k4=46250000.
```

All three polynomial identities and all Mumford remainders are verified
exactly.

## The normalized incidence system

Let all seven sextic coefficients vary.  With `B5(0)=1`, the order-5 block
has eleven auxiliary variables and eleven coefficient equations.  The
order-3 and order-4 blocks each have seven variables and seven equations.
Thus the simultaneous system has

```text
32 variables, 25 equations.
```

Both independent implementations give

```text
auxiliary block ranks = 11,7,7,
full Jacobian rank    = 25,
tangent dimension     = 7.
```

In particular, projection to the seven sextic coefficients has full tangent
rank.  Every rational first-order sextic perturbation has a unique rational
lift of all three marked layers.  After quotienting by the rank-four
`PGL2` plus equation-scaling directions, the marked moduli tangent has the
expected dimension three.

The exact coefficient slice

```text
F_t=F+t*(1+x)
```

adds six linear equations.  Its full Jacobian has rank `31` in `32`
variables, so its pullback to the marked incidence is a smooth algebraic
curve at `t=0`.  This is an algebraic curve with a rational formal branch;
it is not a rational parametrization of the marked curve.

## Transversality to the split locus

Magma gives `#Aut(C)=4`, so the reduced automorphism group is `C2` and there
is only one local Humbert-4 branch.  Its tangent hyperplane in sextic
coefficient space has primitive normal

```text
N=(0,81125,0,904800,0,9916416,0).
```

Consequently

```text
N(1+x)=81125 != 0.
```

The simultaneous marked incidence therefore has genuine rational tangent
directions leaving the geometrically split locus.  The HLP `[60]` point is
not infinitesimally trapped in Humbert `H_4`.

This does not yet give a new rational specialization.  Rational values of
the base parameter must split the finite order-5, order-3, and order-4
support covers simultaneously.  The tangent result removes the geometric
obstruction but leaves a difficult arithmetic rational-point problem.

## Reproduction

```text
magma -b MemGB:=2 code/hlp_z60_marked_torsion_extract.m
magma -b MemGB:=2 code/hlp_z60_marked_identities_verify.m
magma -b MemGB:=2 code/hlp_z60_deformation_system.m
python3 code/hlp_z60_deformation_tangent.py
magma -b MemGB:=2 code/hlp_z60_deformation_seed_geometry.m
```

All completed runs stayed far below the requested memory cap.
